#!/usr/bin/env python3
"""Appendicularian Selenoprofiles profiles for the Sec machinery.

For each family (default: the Selenoprofiles "machinery" set), the standard profile alignment is extended with
appendicularian proteins:
  - stage 1 Selenoprofiles predictions of the family in the training genomes, labelled "homologue" or
    "selenocysteine", whose length is at least --min-fraction of the median length of the profile sequences;
  - curated proteins of the family in curation/machinery_proteins.faa (header: family|genome|source).
The added proteins are aligned with MAFFT --add --keeplength, which keeps the columns of the standard alignment, so
that column-based options of the standard configuration remain valid. The extended alignment is then built into a
profile with "selenoprofiles build", loading the standard configuration (filters, tags) and using all sequences as
BLAST queries. Profiles keep the standard names (SEPHS2, pstk, ...), so that they can replace the standard ones in a
copy of the Selenoprofiles data (see the pipeline option --selenoprofiles_overrides).

Held-out genomes (default: OdB3) contribute no sequence, so that they can test the profiles.
Outputs in --out: <family>.fa and its Selenoprofiles files, and members.tsv (every candidate sequence, with the
reason it was kept or not).
"""
import argparse, os, re, shutil, statistics, subprocess, tempfile

MACHINERY = ["SEPHS2", "pstk", "SecS", "eEFsec", "SBP2", "secp43"]
GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].rstrip("\n"), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", default="data/references/selenoprofiles_data/selenoprotein_profiles")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--curated", default="curation/machinery_proteins.faa")
    ap.add_argument("--families", nargs="+", default=MACHINERY)
    ap.add_argument("--held-out", nargs="*", default=["OdB3"])
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--out", default="work/stage2/machinery_profiles")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    train = [g for g in GENOMES if g not in a.held_out]
    curated = [(h.split("|"), s) for h, s in read_fasta(a.curated)]
    log = [("family", "source", "genome", "name", "length", "kept", "reason")]
    home = tempfile.mkdtemp()  # selenoprofiles reads its configuration from $HOME
    env = dict(os.environ, HOME=home)
    subprocess.run(["selenoprofiles", "-setup"], env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for fam in a.families:
        prof = list(read_fasta(f"{a.profiles}/{fam}.fa"))
        median = statistics.median(len(s.replace("-", "")) for _, s in prof)
        members = []
        for g in GENOMES:
            path = f"{a.stage1}/{g}.selenoprofiles.faa"
            if not os.path.exists(path): continue
            for h, s in read_fasta(path):
                pid = h.split()[0]
                if pid.split(".")[0] != fam: continue
                s = s.replace("*", "")
                label = pid.rsplit(".", 1)[1]
                if g in a.held_out: why = "held-out genome"
                elif label not in ("homologue", "selenocysteine"): why = f"label {label}"
                elif len(s) < a.min_fraction * median: why = f"shorter than {a.min_fraction} x median ({median:.0f})"
                else: why = ""
                log.append((fam, "stage1", g, pid, len(s), "yes" if not why else "no", why or "kept"))
                if not why: members.append((f"{g}_{pid}", s))
        for (f, g, src), s in curated:
            if f != fam: continue
            why = "held-out genome" if g in a.held_out else ""
            log.append((fam, "curated", g, src, len(s), "yes" if not why else "no", why or "kept"))
            if not why: members.append((f"{g}_curated_{len(members)}", s))
        name = fam
        with tempfile.TemporaryDirectory() as d:
            add = f"{d}/add.fa"
            with open(add, "w") as o:
                for n, s in members: o.write(f">{re.sub(r'[^A-Za-z0-9_.]', '_', n)}\n{s}\n")
            merged = f"{d}/{name}.fa"
            with open(merged, "w") as o:
                subprocess.run(["mafft", "--quiet", "--anysymbol", "--keeplength", "--add", add,
                                f"{a.profiles}/{fam}.fa"], stdout=o, check=True)
            out = os.path.join(os.path.abspath(a.out), f"{name}.fa")
            subprocess.run(["selenoprofiles", "build", "-i", merged, "-o", out, "-load", f"{a.profiles}/{fam}.fa.config",
                            "-queries", "all", "-name", name, "-y"], env=env, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        n_out = sum(1 for _ in read_fasta(os.path.join(a.out, f"{name}.fa")))
        print(f"{name}: {len(prof)} standard + {len(members)} appendicularian sequences -> {n_out} after redundancy removal",
              flush=True)
    with open(os.path.join(a.out, "members.tsv"), "w") as o:
        for row in log: o.write("\t".join(map(str, row)) + "\n")
    shutil.rmtree(home, ignore_errors=True)


if __name__ == "__main__":
    main()
