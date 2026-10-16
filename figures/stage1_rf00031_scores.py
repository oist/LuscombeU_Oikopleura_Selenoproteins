#!/usr/bin/env python3
"""Figure: best RF00031 score in the 3' region of each stage 1 Selenoprofiles prediction, per genome.

Predictions labelled selenocysteine are compared with those labelled cysteine (homologs without Sec, whose 3'
regions are not expected to carry a SECIS) and with other labels (pseudo, unaligned, ...). Lines mark the Rfam
gathering threshold of RF00031 (35 bits) and 20 bits. Input: results/tables/stage1_rf00031_scores.tsv."""
import argparse, csv, random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENOMES = [("OKI", "O. dioica\nOKI2018"), ("OSKA", "O. dioica\nOSKA2016"), ("OdB3", "O. dioica\nOdB3"),
           ("Olon", "O. longicauda"), ("Oalb", "O. albicans"), ("Ovan", "O. vanhoeffeni"),
           ("Bsty", "B. stygius"), ("Mery", "M. erythrocephalus"), ("Fbor", "F. borealis")]
GROUPS = [("selenocysteine", "Sec", "#1b7837", -0.22), ("cysteine", "Cys homolog", "#7f7f7f", 0.0),
          ("other", "other label", "#c2a5cf", 0.22)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", default="results/tables/stage1_rf00031_scores.tsv")
    ap.add_argument("--out", default="results/figures/stage1_rf00031_scores")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table), delimiter="\t"))
    rnd = random.Random(1)
    fig, ax = plt.subplots(figsize=(10, 4.2))
    for i, (g, _) in enumerate(GENOMES):
        for label, name, colour, dx in GROUPS:
            sel = [r for r in rows if r["genome"] == g and
                   (r["label"] == label if label != "other" else r["label"] not in ("selenocysteine", "cysteine"))]
            ys = [float(r["best_bits"]) if r["best_bits"] else 0.0 for r in sel]
            xs = [i + dx + rnd.uniform(-0.07, 0.07) for _ in ys]
            ax.scatter(xs, ys, s=10, color=colour, alpha=0.8, linewidths=0, label=name if i == 0 else None)
    ax.axhline(35, color="black", lw=0.8, ls="--")
    ax.text(len(GENOMES) - 0.5, 35.5, "Rfam gathering threshold (35 bits)", ha="right", va="bottom", fontsize=8)
    ax.axhline(20, color="black", lw=0.6, ls=":")
    ax.text(len(GENOMES) - 0.5, 20.5, "20 bits", ha="right", va="bottom", fontsize=8)
    ax.set_xticks(range(len(GENOMES)))
    ax.set_xticklabels([n for _, n in GENOMES], fontsize=8, style="italic")
    ax.set_xlim(-0.6, len(GENOMES) - 0.4)
    ax.set_ylabel("best RF00031 score in the 3′ region (bits)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    import os
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=200)


if __name__ == "__main__":
    main()
