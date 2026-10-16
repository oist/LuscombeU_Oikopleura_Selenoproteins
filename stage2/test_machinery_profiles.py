#!/usr/bin/env python3
"""Search the genomes with the appendicularian machinery profiles and compare with stage 1.

Runs Selenoprofiles with the profiles built by build_machinery_profiles.py on each genome (resumable: a genome whose
GTF exists is skipped), then writes a table of predictions per genome and family, by Selenoprofiles label, next to
the stage 1 predictions of the standard profiles: <out>/comparison.tsv.
"""
import argparse, collections, glob, gzip, os, re, shutil, subprocess, tempfile

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]


def labels(gtf, suffix=""):
    c = collections.defaultdict(collections.Counter)
    seen = set()
    for l in open(gtf):
        m = re.search(r'gene_id "([^"]+)"', l)
        if not m or m.group(1) in seen or ":" in m.group(1): continue  # "Sec1:..." marks Sec positions
        seen.add(m.group(1))
        fam, _, lab = m.group(1).split(".")[0], None, m.group(1).rsplit(".", 1)[1]
        c[fam[:-len(suffix)] if suffix and fam.endswith(suffix) else fam][lab] += 1
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles-dir", default="work/stage2/machinery_profiles")
    ap.add_argument("--data", default="data/references/selenoprofiles_data")
    ap.add_argument("--genomes-dir", default="data/genomes")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--genomes", nargs="+", default=GENOMES)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="work/stage2/machinery_test")
    a = ap.parse_args()
    profiles = sorted(os.path.abspath(p) for p in glob.glob(f"{a.profiles_dir}/*_app.fa"))
    fams = [os.path.basename(p)[:-len("_app.fa")] for p in profiles]
    data = os.path.abspath(a.data)
    os.makedirs(a.out, exist_ok=True)
    home = tempfile.mkdtemp(); env = dict(os.environ, HOME=home)
    subprocess.run(["selenoprofiles", "-setup"], env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for g in a.genomes:
        gtf = os.path.abspath(f"{a.out}/{g}.gtf")
        if os.path.exists(gtf): continue
        w = tempfile.mkdtemp(dir=a.out)
        with gzip.open(f"{a.genomes_dir}/{g}/genome.fna.gz", "rb") as i, open(f"{w}/genome.fna", "wb") as o:
            shutil.copyfileobj(i, o)
        subprocess.run(["selenoprofiles", "-o", f"{w}/res", "-t", f"{w}/genome.fna", "-s", g, "-p", ",".join(profiles),
                        "-selenoprofiles_data_dir", data, "-ncpus", str(a.cpus), "-temp", w,
                        "-output_gtf_file", gtf + ".part", "-output_fasta_file", os.path.abspath(f"{a.out}/{g}.faa"),
                        "-output_three_prime_file", os.path.abspath(f"{a.out}/{g}.three_prime.fa")],
                       env=env, check=True, cwd=w, stdout=open(f"{a.out}/{g}.log", "w"), stderr=subprocess.STDOUT)
        os.replace(gtf + ".part", gtf)
        shutil.rmtree(w)
        print(g, "done", flush=True)
    with open(f"{a.out}/comparison.tsv", "w") as o:
        o.write("genome\tfamily\tstage1_standard_profile\tappendicularian_profile\n")
        for g in a.genomes:
            if not os.path.exists(f"{a.out}/{g}.gtf"): continue
            s1, s2 = labels(f"{a.stage1}/{g}.selenoprofiles.gtf"), labels(f"{a.out}/{g}.gtf", "_app")
            fmt = lambda c: ",".join(f"{n} {lab}" for lab, n in sorted(c.items())) or "-"
            for f in fams: o.write(f"{g}\t{f}\t{fmt(s1.get(f, {}))}\t{fmt(s2.get(f, {}))}\n")
    shutil.rmtree(home, ignore_errors=True)


if __name__ == "__main__":
    main()
