#!/usr/bin/env python3
"""Machinery genes that the appendicularian profiles predict but Selenoprofiles drops at refiltering.

For each genome and family without a kept prediction in the machinery test (stage2/test_machinery_profiles.py),
Selenoprofiles is rerun with that family's profile only, keeping also the predictions dropped by refiltering
(-full_db, -state kept,refiltered). Each candidate protein is compared (blastp) with the proteins of the family found
in the other genomes (kept predictions and curation/machinery_proteins.faa) and with related proteins of other
families (curation/machinery_outgroups.faa, for example KTI12 for PSTK). A candidate whose best hit is a protein of
the family with an E-value of at most --evalue is reported as a fragment of the gene; otherwise as not supported.

Only families with related proteins in the outgroup file are analysed (currently PSTK, with KTI12): without them, a
discarded prediction cannot be told from a related protein (for example, SECp43 from other RNA-binding proteins).

Output: <out>/fragments.tsv (one row per candidate) and results/tables/stage2_machinery_fragments.tsv via make.
"""
import argparse, collections, gzip, os, re, shutil, subprocess, tempfile

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].rstrip("\n"), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def kept(gtf):
    """families with at least one kept prediction (any label) in a Selenoprofiles GTF"""
    fams = set()
    for l in open(gtf):
        m = re.search(r'gene_id "([^":]+)"', l)
        if m: fams.add(m.group(1).split(".")[0].removesuffix("_app"))
    return fams


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles-dir", default="work/stage2/machinery_profiles")
    ap.add_argument("--test-dir", default="work/stage2/machinery_test")
    ap.add_argument("--data", default="data/references/selenoprofiles_data")
    ap.add_argument("--genomes-dir", default="data/genomes")
    ap.add_argument("--curated", default="curation/machinery_proteins.faa")
    ap.add_argument("--outgroups", default="curation/machinery_outgroups.faa")
    ap.add_argument("--evalue", type=float, default=1e-10)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="work/stage2/machinery_fragments")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fams = sorted(os.path.basename(p)[:-3] for p in os.listdir(a.profiles_dir) if p.endswith(".fa"))
    home = tempfile.mkdtemp(); env = dict(os.environ, HOME=home)
    subprocess.run(["selenoprofiles", "-setup"], env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # reference proteins of each family: kept predictions of all genomes, curated proteins; outgroups
    refs = collections.defaultdict(list)
    for g in GENOMES:
        faa = f"{a.test_dir}/{g}.faa"
        if not os.path.exists(faa): continue
        for h, s in read_fasta(faa):
            pid = h.split()[0]
            refs[pid.split(".")[0].removesuffix("_app")].append((f"{g}|{pid}", s.replace("*", "")))
    for h, s in read_fasta(a.curated):
        f, g, _ = h.split("|", 2); refs[f].append((f"{g}|curated", s))
    for h, s in read_fasta(a.outgroups):
        f, name, _ = h.split("|", 2); refs[f"outgroup:{f}"].append((f"outgroup|{name}", s))
    rows = [("genome", "family", "candidate", "length", "contig", "best_hit", "best_hit_evalue", "status")]
    for g in GENOMES:
        gtf = f"{a.test_dir}/{g}.gtf"
        if not os.path.exists(gtf): continue
        missing = [f for f in fams if f not in kept(gtf) and refs.get(f"outgroup:{f}")]
        for f in missing:
            faa = f"{a.out}/{g}.{f}.faa"
            if not os.path.exists(faa):
                w = tempfile.mkdtemp(dir=a.out)
                with gzip.open(f"{a.genomes_dir}/{g}/genome.fna.gz", "rb") as i, open(f"{w}/genome.fna", "wb") as o:
                    shutil.copyfileobj(i, o)
                subprocess.run(["selenoprofiles", "-o", f"{w}/res", "-t", f"{w}/genome.fna", "-s", g,
                                "-p", os.path.abspath(f"{a.profiles_dir}/{f}.fa"), "-selenoprofiles_data_dir",
                                os.path.abspath(a.data), "-ncpus", str(a.cpus), "-temp", w, "-full_db",
                                "-state", "kept,refiltered", "-output_fasta_file", os.path.abspath(faa) + ".part"],
                               env=env, check=True, cwd=w, stdout=open(f"{a.out}/{g}.{f}.log", "w"),
                               stderr=subprocess.STDOUT)
                os.replace(faa + ".part", faa); shutil.rmtree(w)
            cands = [(h.split()[0], h.split()[1].split(":")[-1] if len(h.split()) > 1 else "", s.replace("*", ""))
                     for h, s in read_fasta(faa)]
            if not cands:
                rows.append((g, f, "-", 0, "-", "-", "-", "no candidate")); continue
            with tempfile.TemporaryDirectory() as d:
                with open(f"{d}/db.faa", "w") as o:
                    for n, s in refs[f] + refs[f"outgroup:{f}"]:
                        if not n.startswith(g + "|"): o.write(f">{n}\n{s}\n")
                with open(f"{d}/q.faa", "w") as o:
                    for n, _, s in cands: o.write(f">{n}\n{s}\n")
                subprocess.run(["formatdb", "-i", f"{d}/db.faa", "-p", "T"], check=True, cwd=d)
                res = subprocess.run(["blastall", "-p", "blastp", "-d", f"{d}/db.faa", "-i", f"{d}/q.faa",
                                      "-e", "10", "-m", "8", "-b", "1", "-v", "1"], capture_output=True, text=True,
                                     check=True).stdout
            best = {}
            for l in res.splitlines():
                q, sname, ev = l.split("\t")[0], l.split("\t")[1], float(l.split("\t")[10])
                if q not in best or ev < best[q][1]: best[q] = (sname, ev)
            for n, contig, s in cands:
                hit, ev = best.get(n, ("-", float("inf")))
                ok = hit != "-" and not hit.startswith("outgroup|") and ev <= a.evalue
                rows.append((g, f, n, len(s), contig, hit, f"{ev:.1e}" if hit != "-" else "-",
                             "fragment of the gene" if ok else "not supported"))
            print(g, f, len(cands), "candidates", flush=True)
    with open(f"{a.out}/fragments.tsv", "w") as o:
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    shutil.rmtree(home, ignore_errors=True)


if __name__ == "__main__":
    main()
