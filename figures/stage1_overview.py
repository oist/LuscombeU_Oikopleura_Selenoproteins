#!/usr/bin/env python3
"""Stage 1 overview: what the pipeline finds with its default settings in each genome.

Panels, one row per genome (appendicularians, other tunicates, vertebrates):
  A  selenocysteine predictions (Selenoprofiles) and how many have a SECIS element at the Rfam gathering threshold;
  B  fraction of the checked selenocysteine codons that fall in a correctly annotated coding sequence (genomes with an
     annotation only);
  C  machinery: tRNA-Sec loci supported by two or more tools, and the Selenoprofiles status of each machinery gene.
Inputs: results/tables/stage1_genes_summary.tsv, stage1_annotation.tsv, stage1_machinery.tsv.
"""
import argparse, collections, csv

GROUPS = [("Appendicularians", ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]),
          ("Other tunicates", ["Cint", "Crob", "Csav", "Hror", "Scla", "Tdem"]),
          ("Vertebrates", ["Hsap", "Mmus", "Drer"])]
NAMES = {"OKI": "O. dioica OKI2018", "OSKA": "O. dioica OSKA2016", "OdB3": "O. dioica OdB3", "Olon": "O. longicauda",
         "Oalb": "O. albicans", "Ovan": "O. vanhoeffeni", "Bsty": "B. stygius", "Mery": "M. erythrocephalus",
         "Fbor": "F. borealis", "Cint": "C. intestinalis", "Crob": "C. robusta", "Csav": "C. savignyi",
         "Hror": "H. roretzi", "Scla": "S. clava", "Tdem": "T. democratica", "Hsap": "H. sapiens", "Mmus": "M. musculus",
         "Drer": "D. rerio"}
MACHINERY = ["SEPHS2", "PSTK", "SEPSECS", "EEFSEC", "SECISBP2", "TRNAU1AP"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--figure", default="figures/stage1_overview")
    a = ap.parse_args()
    genes = {r["genome"]: r for r in csv.DictReader(open(f"{a.tables}/stage1_genes_summary.tsv"), delimiter="\t")}
    ann = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(f"{a.tables}/stage1_annotation.tsv"), delimiter="\t"): ann[r["genome"]][r["status"]] += 1
    mach = {}
    for l in open(f"{a.tables}/stage1_machinery.tsv"):
        if l.startswith("#") or l.startswith("Sample"): continue
        f = l.rstrip("\n").split("\t"); mach[f[0]] = f[1:]

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = [g for _, gs in GROUPS for g in gs]
    y = {g: len(order) - 1 - i for i, g in enumerate(order)}
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13, 6.5), sharey=True,
                                        gridspec_kw={"width_ratios": [3, 1.6, 3.2]})
    for g in order:
        r = genes[g]; n, s = int(r["selenocysteine"]), int(r["secis_supported"] or 0)
        axA.barh(y[g], n, color="#d5d8dc", edgecolor="none")
        axA.barh(y[g], s, color="#c0392b", edgecolor="none")
        axA.text(n + 0.5, y[g], f"{s}/{n}", va="center", fontsize=7)
        c = ann.get(g)
        if c:
            tot = sum(c.values()); ok = c.get("Well annotated", 0)
            axB.barh(y[g], ok / tot, color="#2c3e50")
            axB.text(ok / tot + 0.02, y[g], f"{ok}/{tot}", va="center", fontsize=7)
        else:
            axB.text(0.02, y[g], "no annotation", va="center", fontsize=7, color="#7f8c8d")
        m = mach.get(g, [""] * 7)
        trna = m[0]
        axC.text(0, y[g], trna, ha="center", va="center", fontsize=8)
        style = {"Sec": ("o", "#c0392b"), "found": ("o", "#2c3e50"), "partial": ("o", "#bdc3c7"), "-": ("x", "#7f8c8d")}
        for k, v in enumerate(m[1:], start=1):
            mk, col = style.get(v, ("x", "#7f8c8d"))
            axC.scatter(k, y[g], marker=mk, s=45 if mk == "o" else 25, color=col,
                        facecolors="none" if v == "partial" else None, linewidths=1.2)
    for ax in (axA, axB, axC):
        for _, gs in GROUPS[:-1]:
            ax.axhline(y[gs[-1]] - 0.5, color="black", lw=0.6)
    axA.set_yticks([y[g] for g in order]); axA.set_yticklabels([NAMES[g] for g in order], fontsize=8, style="italic")
    axA.set_xlabel("selenocysteine predictions\n(red: with a SECIS at the Rfam threshold)", fontsize=8)
    axB.set_xlim(0, 1.25); axB.set_xticks([0, 0.5, 1]); axB.set_xlabel("fraction in correctly\nannotated genes", fontsize=8)
    axC.set_xticks(range(len(MACHINERY) + 1)); axC.set_xticklabels(["tRNA-Sec\nloci"] + MACHINERY, fontsize=7, rotation=45, ha="right")
    axC.set_xlim(-0.6, len(MACHINERY) + 0.6)
    axC.set_xlabel("machinery (red: SEPHS2 with selenocysteine;\nopen: partial prediction; ×: not found)", fontsize=8)
    for (label, gs) in GROUPS:
        axC.text(len(MACHINERY) + 0.7, (y[gs[0]] + y[gs[-1]]) / 2, label, rotation=270, va="center", fontsize=8)
    for ax in (axA, axB, axC):
        ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="x", labelsize=7)
    fig.suptitle("Stage 1: default settings of the pipeline on appendicularian, other tunicate and vertebrate genomes",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.97, 0.96))
    for ext in ("png", "pdf"): fig.savefig(f"{a.figure}.{ext}", dpi=200)


if __name__ == "__main__":
    main()
