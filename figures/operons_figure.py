#!/usr/bin/env python3
"""Figure OPERONS (step 10): selenoprotein genes in the operons of the three O. dioica assemblies.

A. Schematic of an operon (after Ganot et al. 2004): a polycistronic pre-mRNA resolved into mRNAs by spliced-leader
   (SL) trans-splicing at the 5' end of each downstream mRNA and by cleavage and polyadenylation of the upstream one,
   23 to 30 nt before the next trans-splice site; a selenoprotein mRNA must hold its SECIS element in that 3' UTR.
B. Each selenoprotein gene followed by another gene on the same strand within --max-distance nt, in each assembly
   (results/tables/stage2_operons.tsv), aligned on its stop codon: the end of its coding sequence, its 3' UTR in its
   transcripts, its SECIS element (stem-loop), the window where its mRNA should end (23-30 nt before the first
   trans-splice site of the next gene), the trans-splice sites of the next gene (SL flags) and the next coding
   sequence; sorted by the distance to the next coding sequence.
C. Fraction of the loci in operons and, among them, fraction first in their operon: selenoprotein genes against all
   annotated loci (results/tables/stage2_operons_summary.tsv), with Fisher's exact test.
Outputs: <out>.pdf, <out>.svg and <out>.png.
"""
import argparse, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch

FS = 7
COL = {"OKI": "#1b7837", "OSKA": "#2166ac", "OdB3": "#b2182b"}
NAME = {"OKI": "OKI2018", "OSKA": "OSKA2016", "OdB3": "OdB3"}
GREY = "#9e9e9e"


def hairpin(ax, x, y, h, color):
    """A small stem-loop glyph standing on (x, y)."""
    ax.plot([x - 0.15, x - 0.15, x + 0.15, x + 0.15], [y, y + h, y + h, y], color=color, lw=1.2)
    ax.add_patch(plt.Circle((x, y + h + 0.17), 0.17, fill=False, color=color, lw=1.2))


def panel_a(ax):
    ax.set_xlim(0, 30); ax.set_ylim(0, 10); ax.axis("off")
    # genome / pre-mRNA
    y = 8.2
    ax.plot([0.5, 29.5], [y, y], color="black", lw=0.8)
    genes = [(2, 8, "gene 1", "#d9d9d9"), (11, 17, "selenoprotein", "#f4a582"), (21.5, 27.5, "gene 3", "#d9d9d9")]
    for a, b, lab, c in genes:
        ax.add_patch(Polygon([[a, y - 0.45], [b - 0.6, y - 0.45], [b, y], [b - 0.6, y + 0.45], [a, y + 0.45]],
                             closed=True, fc=c, ec="black", lw=0.6))
        ax.text((a + b) / 2, y + 0.75, lab, ha="center", va="bottom", fontsize=FS)
    ax.text(14.0, y - 0.05, "UGA", ha="center", va="center", fontsize=FS - 1.5)
    hairpin(ax, 18.7, y, 0.45, "#b2182b")
    ax.text(0.5, y + 1.9, "operon (pre-mRNA)", fontsize=FS, ha="left", fontweight="bold")
    ax.annotate("", xy=(15, 6.2), xytext=(15, 7.3), arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.text(15.4, 6.75, "SL trans-splicing,\ncleavage and polyadenylation", fontsize=FS - 1, va="center")
    # mRNAs
    for k, (a, b, lab, c) in enumerate(genes):
        yy = 4.6 - 1.7 * k
        ax.add_patch(Rectangle((a - 0.9, yy - 0.2), 0.6, 0.4, fc="black"))
        ax.plot([a - 0.3, a], [yy, yy], color="black", lw=0.8)
        ax.add_patch(Rectangle((a, yy - 0.3), b - a, 0.6, fc=c, ec="black", lw=0.6))
        utr = 2.6 if k == 1 else 1.4
        ax.plot([b, b + utr], [yy, yy], color="black", lw=0.8)
        ax.text(b + utr + 0.1, yy, "AAA", fontsize=FS - 1.5, va="center")
        if k == 1:
            ax.text(14.0, yy, "UGA", ha="center", va="center", fontsize=FS - 1.5)
            hairpin(ax, b + 1.5, yy, 0.45, "#b2182b")
            ax.text(b + 1.5, yy - 0.75, "SECIS", ha="center", fontsize=FS - 1, color="#b2182b")
    ax.text(0.2, 4.6 - 0.2, "SL", fontsize=FS - 1.5, va="top")
    ax.text(0.5, 0.2, "polyadenylation site 23-30 nt before the next trans-splice site (Ganot et al. 2004)",
            fontsize=FS - 1, ha="left")


def secis_glyph(ax, s0, s1, y, color, h=0.42):
    """SECIS element as a stem-loop drawn over [s0, s1]: two stems joined by a loop."""
    m, w = (s0 + s1) / 2, (s1 - s0)
    xs = [s0, s0 + 0.3 * w, s0 + 0.3 * w, m, s1 - 0.3 * w, s1 - 0.3 * w, s1]
    ys = [y, y, y + 0.55 * h, y + h, y + 0.55 * h, y, y]
    ax.plot(xs, ys, color=color, lw=1.4, solid_capstyle="round", solid_joinstyle="round")


def gene_arrow(ax, x0, x1, y, fc, ec="black", h=0.30):
    """A coding sequence drawn as an arrow from x0 to x1 (pointing right)."""
    tip = min(14, (x1 - x0) * 0.5)
    ax.add_patch(Polygon([[x0, y - h], [x1 - tip, y - h], [x1, y], [x1 - tip, y + h], [x0, y + h]], closed=True,
                         fc=fc, ec=ec, lw=0.5))


def panel_b(ax, rows, max_d):
    """One row per selenoprotein gene followed by a gene on the same strand: end of its coding sequence (arrow, left of
    0), the 3' UTR of its gold-standard model (black line; transcripts can run beyond the end of the mRNA), its SECIS element (stem-loop), the window where its mRNA
    should end (23-30 nt before the first trans-splice site of the next gene, shaded), the trans-splice sites of the
    next gene (SL flags: its mRNA starts there) and the next coding sequence (grey arrow)."""
    sel = []
    for r in rows:
        if r["next_gene_same_strand_after_stop"] == "" or not r["secis_after_stop"]: continue
        d = int(r["next_gene_same_strand_after_stop"])
        if d > max_d: continue
        s0, s1 = map(int, r["secis_after_stop"].split("-"))
        sl = [int(x.rstrip("*")) for x in r["next_gene_sl_sites_after_stop"].split(",") if x]
        utr = int(r["utr3_length"]) if r["utr3_length"] not in ("", "0") else None
        sel.append((d, r["gene"], r["genome"], s0, s1, sorted(set(sl)), utr))
    sel.sort(key=lambda x: (x[0], x[1]))
    right = max_d * 1.08
    for i, (d, gene, g, s0, s1, sl, utr) in enumerate(sel):
        y = len(sel) - i
        gene_arrow(ax, -70, 0, y, COL[g], h=0.26)                     # end of the selenoprotein coding sequence
        ax.plot([0, d], [y, y], color="#cfcfcf", lw=0.8)              # genome between the two coding sequences
        if utr: ax.plot([0, min(utr, d)], [y, y], color="black", lw=1.6, solid_capstyle="butt")
        if sl:
            f = min(sl)
            ax.add_patch(Rectangle((f - 30, y - 0.33), 7, 0.66, fc="#fde0b0", ec="none"))
            for x in sl:                                               # SL flags of the next mRNA
                ax.plot([x, x], [y, y + 0.45], color="black", lw=0.8)
                ax.add_patch(Rectangle((x, y + 0.30), 14, 0.15, fc="black"))
        if s1 <= d: secis_glyph(ax, s0, s1, y, COL[g])
        gene_arrow(ax, d, min(d + 60, right), y, "#d0d0d0")           # next coding sequence
        label = f"{gene} ({NAME[g]})" + (" *" if s1 > d else "")
        ax.text(-80, y, label, ha="right", va="center", fontsize=FS - 1)
    # key
    yk, xk, xt = len(sel) + 3.6, 360, 435
    gene_arrow(ax, xk, xk + 55, yk, COL["OKI"], h=0.26)
    ax.text(xt, yk, "selenoprotein coding sequence (colour: assembly)", fontsize=FS - 1, va="center")
    ax.plot([xk, xk + 55], [yk - 1, yk - 1], color="black", lw=1.6)
    ax.text(xt, yk - 1, "3' UTR of the gold-standard model", fontsize=FS - 1, va="center")
    secis_glyph(ax, xk + 5, xk + 55, yk - 2.2, COL["OKI"])
    ax.text(xt, yk - 2, "SECIS element", fontsize=FS - 1, va="center")
    ax.add_patch(Rectangle((xk + 20, yk - 3.33), 7, 0.66, fc="#fde0b0", ec="none"))
    ax.text(xt, yk - 3, "expected end of its mRNA (23-30 nt before the next SL)", fontsize=FS - 1,
            va="center")
    ax.plot([xk + 20, xk + 20], [yk - 4.2, yk - 3.75], color="black", lw=0.8)
    ax.add_patch(Rectangle((xk + 20, yk - 3.9), 14, 0.15, fc="black"))
    ax.text(xt, yk - 4, "trans-splice site of the next gene", fontsize=FS - 1, va="center")
    gene_arrow(ax, xk, xk + 55, yk - 5, "#d0d0d0")
    ax.text(xt, yk - 5, "next coding sequence, same strand", fontsize=FS - 1, va="center")
    ax.text(xk, yk - 6, "* SECIS only beyond the next coding sequence", fontsize=FS - 1, va="center")
    ax.set_xlim(-75, right); ax.set_ylim(0.3, len(sel) + 4.4)
    ax.set_yticks([]); ax.spines[["left", "top", "right"]].set_visible(False)
    ax.spines["bottom"].set_bounds(0, max_d)
    ax.set_xticks(range(0, max_d + 1, 100))
    ax.set_xlabel("nucleotides after the stop codon of the selenoprotein gene", fontsize=FS)
    ax.tick_params(labelsize=FS - 1)
    return len(sel)


def panel_c(ax, summary):
    gs = [r["genome"] for r in summary]
    x = range(len(gs))
    w = 0.2
    for k, r in enumerate(summary):
        n, op, first = int(r["loci"]), int(r["loci_in_operons"]), int(r["first_in_operon"])
        sn, sop, sfirst = int(r["selenoprotein_genes"]), int(r["selenoprotein_in_operons"]), int(r["selenoprotein_first"])
        vals = [(op / n, "all loci", "#cccccc"), (sop / sn, "selenoprotein genes", COL[r["genome"]])]
        vals2 = [(first / op, "#e8e8e8"), (sfirst / sop if sop else 0, COL[r["genome"]])]
        ax.bar(k - 1.5 * w, vals[0][0], w, color=vals[0][2], edgecolor="black", lw=0.4)
        ax.bar(k - 0.5 * w, vals[1][0], w, color=vals[1][2], edgecolor="black", lw=0.4)
        ax.bar(k + 0.5 * w, vals2[0][0], w, color=vals2[0][1], edgecolor="black", lw=0.4, hatch="///")
        ax.bar(k + 1.5 * w, vals2[1][0], w, color=vals2[1][1], edgecolor="black", lw=0.4, hatch="///", alpha=0.6)
        ax.text(k - 0.5 * w, sop / sn + 0.02, f"{sop}/{sn}", ha="center", fontsize=FS - 2)
        ax.text(k + 1.5 * w, (sfirst / sop if sop else 0) + 0.02, f"{sfirst}/{sop}", ha="center", fontsize=FS - 2)
        ax.text(k, -0.2, f"P = {float(r['fisher_p_operon']):.2g}; {float(r['fisher_p_first_among_operonic']):.2g}",
                ha="center", fontsize=FS - 2, transform=ax.get_xaxis_transform())
    ax.set_xticks(list(x)); ax.set_xticklabels([NAME[g] for g in gs], fontsize=FS)
    ax.set_ylim(0, 1); ax.set_ylabel("fraction", fontsize=FS); ax.tick_params(labelsize=FS - 1)
    ax.spines[["top", "right"]].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="#cccccc", ec="black", lw=0.4, label="in operons, all loci"),
                       Patch(fc="white", ec="black", lw=0.4, label="in operons, selenoprotein genes (colour)"),
                       Patch(fc="#e8e8e8", ec="black", lw=0.4, hatch="///", label="first in operon, among loci in operons")],
              fontsize=FS - 2, frameon=False, loc="upper left")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default="results/tables/stage2_operons.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_operons_summary.tsv")
    ap.add_argument("--max-distance", type=int, default=800)
    ap.add_argument("--out", default="figures/operons")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table), delimiter="\t"))
    summary = list(csv.DictReader(open(a.summary), delimiter="\t"))
    fig = plt.figure(figsize=(7.2, 7.4))
    axa = fig.add_axes([0.02, 0.70, 0.60, 0.28])
    axc = fig.add_axes([0.70, 0.74, 0.28, 0.22])
    axb = fig.add_axes([0.25, 0.06, 0.72, 0.60])
    panel_a(axa); panel_c(axc, summary); n = panel_b(axb, rows, a.max_distance)
    for ax, lab, (x, y) in ((axa, "A", (0.0, 0.985)), (axc, "C", (0.64, 0.985)), (axb, "B", (0.0, 0.66))):
        fig.text(x, y, lab, fontsize=FS + 3, fontweight="bold", va="top")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=200)
    print(f"{n} genes in panel B -> {a.out}.pdf")


if __name__ == "__main__":
    main()
