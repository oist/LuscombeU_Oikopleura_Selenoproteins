#!/usr/bin/env python3
"""Which filter of a Selenoprofiles profile removes a prediction: run Selenoprofiles with one profile on one genome
and report, for every candidate prediction before filtering, the measures and the value of each clause of the
profile's filters (stage2/refiltering_hook.py, called by a pre_filtering action of Selenoprofiles).

Jobs (--jobs): genome:family:build, where build is a directory of stage 2 profiles (work/stage2/profiles/<build>) or
"standard" (the Selenoprofiles data). Default: SELENOS in the three gold-standard and the three new O. dioica
assemblies with the standard, E0 and LOO_OKI profiles (the benchmark of 2026-10-28 found O. dioica SELENOS with the
standard profile in OSKA2016 in stage 1 and with LOO_OKI in OKI2018 and KUM, but with no other build), and GPX and
DIO in O. albicans (two stage 1 references not found with E0).
Selenoprofiles runs with the full Selenoprofiles data (Gene Ontology and tag database, used by some filters) and the
profile folder of the build. Outputs: --out/<genome>.<family>.<build>.tsv and .faa (the proteins of the candidate
predictions), --table (all jobs, with the reference gene, if any, that overlaps each prediction) and, with --faa, the
proteins of all jobs (named <genome>|<family>|<build>|<prediction>).
"""
import argparse, concurrent.futures, csv, gzip, os, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = [f"{g}:SELENOS:{b}" for g in ["OKI", "OSKA", "OdB3", "KUM", "AOM", "Bar2"]
           for b in ["standard", "E0", "LOO_OKI"]] + ["Oalb:GPX:standard", "Oalb:GPX:E0", "Oalb:DIO:standard",
                                                      "Oalb:DIO:E0"]


def job(spec, a):
    genome, family, build = spec.split(":")
    folder = os.path.abspath(f"{a.data}/selenoprotein_profiles" if build == "standard" else f"{a.profiles}/{build}")
    out = os.path.abspath(f"{a.out}/{genome}.{family}.{build}.tsv")
    if os.path.exists(out) and not a.force: return spec, out, "kept from an earlier run"
    for f in (out, out[:-4] + ".faa"):
        if os.path.exists(f): os.remove(f)
    work = tempfile.mkdtemp(dir=os.path.abspath(a.out))
    try:
        env = dict(os.environ, HOME=work, DIAG_CONFIG=f"{folder}/{family}.fa.config", DIAG_OUT=out,
                   DIAG_FAA=out[:-4] + ".faa")
        subprocess.run(["selenoprofiles", "-setup"], env=env, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        # Selenoprofiles 4.5.8 refuses actions on the command line (strict option parsing); its configuration file
        # accepts them
        action = f"import sys; sys.path.insert(0, {HERE!r}); import refiltering_hook; refiltering_hook.report(x)"
        with open(f"{work}/.selenoprofiles_config.txt", "a") as o:
            o.write(f"\nACTION.pre_filtering.zz_diagnose = {action}\n")
        target = f"{work}/{genome}.fa"
        with gzip.open(f"{a.genomes}/{genome}/genome.fna.gz", "rb") as i, open(target, "wb") as o:
            shutil.copyfileobj(i, o)
        p = subprocess.run(["selenoprofiles", "-o", f"{work}/results", "-t", target, "-s", genome, "-p", family,
                            "-selenoprofiles_data_dir", os.path.abspath(a.data), "-profiles_folder", folder,
                            "-temp", work, "-ncpus", "1"],
                           env=env, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode:
            return spec, out, f"selenoprofiles failed ({p.returncode}): {p.stderr[-1500:]}"
        if not os.path.exists(out):
            return spec, out, "no candidate prediction reported (action not run?)"
        return spec, out, "done"
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", nargs="+", default=DEFAULT)
    ap.add_argument("--profiles", default="work/stage2/profiles")
    ap.add_argument("--data", default="data/references/selenoprofiles_data")
    ap.add_argument("--genomes", default="data/genomes")
    ap.add_argument("--truth", default="results/tables/stage2_truth_sets.tsv")
    ap.add_argument("--out", default="work/stage2/diagnose_refiltering")
    ap.add_argument("--table", default="results/tables/stage2_refiltering_diagnosis.tsv")
    ap.add_argument("--faa", default="", help="proteins of the candidate predictions of all jobs")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--force", action="store_true", help="rerun jobs whose output exists")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(a.threads) as ex:
        for spec, out, status in ex.map(lambda s: job(s, a), a.jobs):
            print(f"{spec}: {status}", flush=True)
    truth = [r for r in csv.DictReader(open(a.truth), delimiter="\t") if r["cds"]]

    def reference(genome, contig, strand, lo, hi):
        for r in truth:
            if r["genome"] != genome or r["contig"] != contig or r["strand"] != strand: continue
            segs = [tuple(map(int, s.split("-"))) for s in r["cds"].split(";")]
            if lo <= segs[-1][1] and hi >= segs[0][0]: return r["gene"]
        return ""

    from refiltering_hook import MEASURES
    head = ["genome", "family", "build", "reference_gene", "prediction", "contig", "strand", "start", "end",
            "profile_sequences", "label", "coverage", "identity_with_profile", "awsi_z_score", "awsi_z_score_w",
            "sec_is_aligned", "p2g_filtering", "p2g_refiltering", "refiltering_clauses_false"]
    rows = []
    for spec in a.jobs:
        genome, family, build = spec.split(":")
        path = f"{a.out}/{genome}.{family}.{build}.tsv"
        if not os.path.exists(path) or not os.path.getsize(path): continue
        reader = list(csv.reader(open(path), delimiter="\t"))
        cols = reader[0]
        n = 5 + len(MEASURES)                    # prediction, position, then the measures
        for r in reader[1:]:
            d = dict(zip(cols, r))
            try: ref = reference(genome, d["contig"], d["strand"], int(d["start"]), int(d["end"]))
            except ValueError: ref = ""
            failing = [c for c, v in zip(cols[n + 2:], r[n + 2:]) if v != "True"]
            rows.append([genome, family, build, ref] + r[:n] + [r[n], r[n + 1], " | ".join(failing)])
    if a.faa:
        with open(a.faa, "w") as o:
            for spec in a.jobs:
                genome, family, build = spec.split(":")
                path = f"{a.out}/{genome}.{family}.{build}.faa"
                if not os.path.exists(path): continue
                name = None
                for l in open(path):
                    if l.startswith(">"): name = l[1:].strip()
                    elif name: o.write(f">{genome}|{family}|{build}|{name}\n{l.strip()}\n")
    if rows and a.table:
        with open(a.table, "w") as o:
            o.write("\t".join(head) + "\n")
            for r in rows: o.write("\t".join(r) + "\n")
        print(f"{len(rows)} predictions -> {a.table}")


if __name__ == "__main__":
    main()
