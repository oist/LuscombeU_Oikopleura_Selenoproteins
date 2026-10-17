#!/usr/bin/env python3
"""Is a machinery gene absent from a genome, or only too divergent for its profile?

References: the sequences of the gene's Selenoprofiles profile and the stage 1 predictions of the gene (labels
"homologue" or "selenocysteine") in all genomes. Outgroups: proteins of related families that a search can confuse
with the gene, taken from the Selenoprofiles UniRef50 database by a description pattern (--outgroup-pattern) and/or
from a FASTA file (--outgroup-faa; header "<gene>|<name>|<description>", only entries of the gene are used).

For each target (genome or transcriptome), the stage 1 predictions of the gene from the other genomes are searched with
tblastn. Hits within --locus-gap nucleotides on the same strand of a sequence form a locus; the translations of its
hits (genomes) or the open reading frame around its best hit (transcripts) give the candidate protein. Each candidate
is compared by blastp with the references (except those of its own genome) and the outgroups: its margin is the score
of its best reference hit minus that of its best outgroup hit. Positive controls are the stage 1 predictions of the
gene themselves, scored the same way against the references of the other genomes.

A candidate is "supported" when its margin reaches --min-margin bits, "ambiguous" (to be resolved by phylogeny) when it
reaches --ambiguous-margin bits, and a "related protein" otherwise. Divergent but genuine genes can fall in the
ambiguous range: for PSTK, the stage 1 prediction of O. vanhoeffeni scores a margin of about 11 bits.
Outputs in --out: <gene>.candidates.tsv (one row per candidate), <gene>.controls.tsv, <gene>.summary.tsv (per target).
"""
import argparse, collections, glob, gzip, os, re, shutil, subprocess, tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
APPENDICULARIANS = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
TRANSCRIPTOMES = {"OKI": "OKI_trinity", "OSKA": "OSKA_tx", "OdB3": "OdB3_est"}


def tr(s): return "".join(CODE.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].rstrip("\n"), []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def run(cmd, **kw): return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def blastp(query_faa, db_faa, d):
    run(["formatdb", "-i", db_faa, "-p", "T"], cwd=d)
    best = {}
    for l in run(["blastall", "-p", "blastp", "-d", db_faa, "-i", query_faa, "-e", "10", "-m", "8", "-b", "1000",
                  "-v", "1000"], cwd=d).splitlines():
        f = l.split("\t"); q, s, bits = f[0], f[1], float(f[11])
        if q.split("|")[1] == s.split("|")[1] and not s.startswith("OUT|"): continue  # same genome
        kind = "OUT" if s.startswith("OUT|") else "REF"
        if bits > best.get((q, kind), (0, ""))[0]: best[(q, kind)] = (bits, s)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True, help="Selenoprofiles family name, for example secp43 or pstk")
    ap.add_argument("--outgroup-pattern", default="", help="regular expression on UniRef50 descriptions")
    ap.add_argument("--outgroup-exclude", default="", help="regular expression of descriptions to exclude")
    ap.add_argument("--outgroup-faa", default="")
    ap.add_argument("--targets", nargs="+", default=APPENDICULARIANS)
    ap.add_argument("--data", default="data")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--min-margin", type=float, default=30.0)
    ap.add_argument("--ambiguous-margin", type=float, default=10.0)
    ap.add_argument("--locus-gap", type=int, default=20000)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="work/stage2/machinery_absence")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    d = tempfile.mkdtemp(dir=a.out)
    sp = f"{a.data}/references/selenoprofiles_data"
    # references and stage 1 predictions of the gene
    refs, preds = [], []
    for n, s in read_fasta(f"{sp}/selenoprotein_profiles/{a.gene}.fa"):
        refs.append((f"REF|profile|{n.split()[0]}", s.replace("-", "")))
    for faa in sorted(glob.glob(f"{a.stage1}/*.selenoprofiles.faa")):
        g = os.path.basename(faa).split(".")[0]
        for n, s in read_fasta(faa):
            pid = n.split()[0]
            if pid.split(".")[0] == a.gene and pid.rsplit(".", 1)[1] in ("homologue", "selenocysteine"):
                refs.append((f"REF|{g}|{pid}", s.replace("*", ""))); preds.append((f"CTRL|{g}|{pid}", s.replace("*", "")))
    outgroups = []
    if a.outgroup_pattern:
        for n, s in read_fasta(f"{sp}/uniref50_seleno_sim.fasta"):
            if re.search(a.outgroup_pattern, n) and not (a.outgroup_exclude and re.search(a.outgroup_exclude, n)):
                outgroups.append((f"OUT|uniref|{n.split()[0]}", s))
    if a.outgroup_faa:
        for n, s in read_fasta(a.outgroup_faa):
            f = n.split("|")
            if f[0] == a.gene: outgroups.append((f"OUT|file|{f[1]}", s))
    db = f"{d}/refs.faa"
    with open(db, "w") as o:
        for n, s in refs + outgroups: o.write(f">{n}\n{s}\n")
    # positive controls
    with open(f"{d}/ctrl.faa", "w") as o:
        for n, s in preds: o.write(f">{n}\n{s}\n")
    best = blastp("ctrl.faa", "refs.faa", d)
    with open(f"{a.out}/{a.gene}.controls.tsv", "w") as o:
        o.write("protein\tlength\tref_bits\toutgroup_bits\tmargin\n")
        for n, s in preds:
            r, x = best.get((n, "REF"), (0, ""))[0], best.get((n, "OUT"), (0, ""))[0]
            o.write(f"{n}\t{len(s)}\t{r}\t{x}\t{r - x:.1f}\n")
    print(f"{a.gene}: {len(refs)} references, {len(outgroups)} outgroups, {len(preds)} positive controls", flush=True)
    # targets
    rows, summary = [], []
    for g in a.targets:
        for kind, path in [("genome", f"{a.data}/genomes/{g}/genome.fna.gz")] + \
                          ([("transcriptome", f"{a.data}/transcripts/{TRANSCRIPTOMES[g]}.fa.gz")] if g in TRANSCRIPTOMES else []):
            if not os.path.exists(path): continue
            seqs = {n.split()[0]: s for n, s in read_fasta(path)}
            with open(f"{d}/t.fa", "w") as o:
                for n, s in seqs.items(): o.write(f">{n}\n{s}\n")
            run(["formatdb", "-i", "t.fa", "-p", "F"], cwd=d)
            with open(f"{d}/q.faa", "w") as o:
                for n, s in preds:
                    if n.split("|")[1] != g: o.write(f">{n}\n{s}\n")
            hits = collections.defaultdict(list)
            for l in run(["blastall", "-p", "tblastn", "-d", "t.fa", "-i", "q.faa", "-e", "1e-3", "-m", "8",
                          "-a", str(a.cpus)], cwd=d).splitlines():
                f = l.split("\t"); s0, s1 = int(f[8]), int(f[9])
                hits[(f[1], "+" if s1 > s0 else "-")].append((min(s0, s1), max(s0, s1), float(f[11])))
            cands = []
            for (sub, strand), hs in hits.items():
                hs.sort(); locus = [hs[0]]
                for h in hs[1:] + [None]:
                    if h is not None and h[0] - locus[-1][1] <= a.locus_gap:
                        locus.append(h); continue
                    seq = seqs[sub]
                    if kind == "transcriptome":
                        b = max(locus, key=lambda x: x[2]); nt = seq if strand == "+" else seq.translate(COMP)[::-1]
                        st = b[0] - 1 if strand == "+" else len(seq) - b[1]
                        fr = st % 3; prot = tr(nt[fr:]); p = (st - fr) // 3
                        i = prot.rfind("*", 0, p) + 1; j = prot.find("*", p); prot = prot[i:(j if j >= 0 else len(prot))]
                    else:  # translations of the non-overlapping hits, in gene order
                        parts, last = [], 0
                        for x in sorted(locus, key=lambda x: -x[2]):
                            if all(x[1] < y[0] or x[0] > y[1] for y in parts): parts.append(x)
                        parts.sort(key=lambda x: x[0] if strand == "+" else -x[0])
                        prot = "".join(tr(seq[x[0] - 1:x[1]] if strand == "+" else seq[x[0] - 1:x[1]].translate(COMP)[::-1])
                                       for x in parts)
                    cands.append((f"CAND|{g}|{kind}:{sub}:{locus[0][0]}-{locus[-1][1]}:{strand}", prot.replace("*", "")))
                    if h is not None: locus = [h]
            if cands:
                with open(f"{d}/c.faa", "w") as o:
                    for n, s in cands: o.write(f">{n}\n{s}\n")
                bc = blastp("c.faa", "refs.faa", d)
            top = None
            for n, s in cands:
                r, rs = bc.get((n, "REF"), (0, "-")); x, xs = bc.get((n, "OUT"), (0, "-"))
                m = r - x
                status = "supported" if m >= a.min_margin else "ambiguous" if m >= a.ambiguous_margin else "related protein"
                rows.append((g, kind, n.split("|", 2)[2], len(s), r, rs, x, xs, f"{m:.1f}", status))
                if top is None or m > top[0]: top = (m, n.split("|", 2)[2], status)
            summary.append((g, kind, len(cands), sum(1 for r in rows if r[0] == g and r[1] == kind and r[9] == "supported"),
                            sum(1 for r in rows if r[0] == g and r[1] == kind and r[9] == "ambiguous"),
                            f"{top[0]:.1f}" if top else "-", top[1] if top else "-"))
            print(g, kind, len(cands), "candidates", flush=True)
    with open(f"{a.out}/{a.gene}.candidates.tsv", "w") as o:
        o.write("target\ttype\tlocus\tlength\tref_bits\tbest_reference\toutgroup_bits\tbest_outgroup\tmargin\tstatus\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    with open(f"{a.out}/{a.gene}.summary.tsv", "w") as o:
        o.write("target\ttype\tcandidates\tsupported\tambiguous\tbest_margin\tbest_candidate\n")
        for r in summary: o.write("\t".join(map(str, r)) + "\n")
    shutil.rmtree(d)


if __name__ == "__main__":
    main()
