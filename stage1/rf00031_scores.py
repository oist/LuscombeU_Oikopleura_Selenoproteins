#!/usr/bin/env python3
"""Scores of the Rfam SECIS model RF00031 in the 3' regions of all stage 1 Selenoprofiles predictions.

Stage 1 searched these regions at the Rfam gathering threshold, which reports nothing below 35 bits. To see the
whole distribution, this script searches them again with no score threshold and with all filters off (cmsearch
--max -T 0, top strand only, since each region is written in the orientation of its gene), and keeps the best hit
of each region. Writes results/tables/stage1_rf00031_scores.tsv: genome, prediction, family, label (Selenoprofiles
label: selenocysteine, cysteine, ...), best score (bits), E-value, hit position in the region, hits >= 10 bits.
Needs Infernal (cmsearch)."""
import argparse, glob, os, subprocess, tempfile

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="work/stage1/pipeline")
    ap.add_argument("--cm", default="data/references/RF00031.cm")
    ap.add_argument("--out", default="results/tables/stage1_rf00031_scores.tsv")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    a = ap.parse_args()
    rows = []
    for fa in sorted(glob.glob(f"{a.outdir}/selenoprofiles/*.three_prime.fa")):
        genome = os.path.basename(fa)[: -len(".three_prime.fa")]
        regions = [l[1:].split()[0] for l in open(fa) if l.startswith(">")]
        best = {r: None for r in regions}
        n10 = {r: 0 for r in regions}
        with tempfile.NamedTemporaryFile(suffix=".tbl") as tbl:
            subprocess.run(["cmsearch", "--max", "-T", "0", "--toponly", "--cpu", str(a.cpus), "--tblout", tbl.name,
                            a.cm, fa], check=True, stdout=subprocess.DEVNULL)
            for l in open(tbl.name):
                if l.startswith("#"): continue
                f = l.split()
                score = float(f[14])
                if score >= 10: n10[f[0]] += 1
                if best[f[0]] is None or score > best[f[0]][0]:
                    best[f[0]] = (score, f[15], int(f[7]), int(f[8]))
        for r in regions:
            pred = r.split(":", 1)[1]
            family, label = pred.split(".")[0], pred.split(".")[-1]
            b = best[r]
            rows.append([genome, pred, family, label] + ([f"{b[0]:.1f}", b[1], str(b[2]), str(b[3])] if b else ["", "", "", ""]) + [str(n10[r])])
        print(f"{genome}: {len(regions)} regions", flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as out:
        out.write("genome\tprediction\tfamily\tlabel\tbest_bits\tevalue\thit_start\thit_end\thits_ge10\n")
        out.write("".join("\t".join(r) + "\n" for r in rows))

if __name__ == "__main__":
    main()
