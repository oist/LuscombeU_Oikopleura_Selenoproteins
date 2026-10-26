#!/usr/bin/env python3
"""Figures of the Oikopleura SECIS model (step 3) and of its scores in genomes and controls.

figures/secis_model.{pdf,svg,png} (from stage2/secis_model_description.py):
  A. consensus sequence and structure: arcs for the 18 base pairs, coloured by the fraction of seeds with a canonical
     pair; letters shaded by the fraction of seeds with the most frequent residue; core motifs underlined.
  B. identity between seeds (consensus positions): all pairs of seeds, and those of orthologous genes of the three
     O. dioica assemblies.
  C. canonical pairs in the seeds against the best hits of the model in shuffled windows (decoys), aligned to the
     model in the same way, per base pair.
  D. changes between the orthologous SECIS of the three O. dioica assemblies: fraction of changes from a canonical
     pair that keep a canonical pair, observed and expected for random new bases; fraction of sites that differ.

figures/secis_scores.{pdf,svg,png} (from the step 4 tables):
  A. hits per Mb (both strands) at or above each score, in each genome and in its window-shuffled copy.
  B. hits at the genome threshold, by what lies upstream on the same strand: a selenoprotein gene (gold standard or
     stage 1 selenocysteine prediction), a copy of a known selenoprotein found by the comparative screen, KKSU (step
     5 loci), a stage 1 prediction without selenocysteine, a TGA pair shared with another species, or nothing found.
     The classes are also written to results/tables/stage2_secis_scan_classes.tsv.
  C. controls: best score per RefSeq mRNA with a frame through a TGA upstream of a hit of 20 bits or more, human and
     C. intestinalis, Oikopleura model and RF00031, for selenoprotein mRNAs and the others.
"""
import argparse, collections, csv, itertools, os, re, zlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

NAMES = {"OKI": "O. dioica OKI2018", "OSKA": "O. dioica OSKA2016", "OdB3": "O. dioica OdB3", "Olon": "O. longicauda",
         "Oalb": "O. albicans", "Ovan": "O. vanhoeffeni", "Bsty": "B. stygius", "Mery": "M. erythrocephalus",
         "Fbor": "F. borealis", "Cint": "C. intestinalis", "Csav": "C. savignyi", "Scla": "S. clava",
         "Tdem": "T. democratica", "Hsap": "H. sapiens"}
GROUPS = [("O. dioica", ["OKI", "OSKA", "OdB3"]), ("other Oikopleura", ["Olon", "Oalb", "Ovan"]),
          ("held-out appendicularians", ["Bsty", "Mery", "Fbor"]), ("non-target tunicates", ["Cint", "Csav", "Scla", "Tdem"])]
HELIX_COLOURS = {"1": "#4C72B0", "2": "#DD8452", "3": "#55A868"}
CLASSES = [("selenoprotein gene", "#C44E52"), ("copy of a known selenoprotein", "#DD8452"), ("KKSU", "#8172B3"),
           ("stage 1 prediction without selenocysteine", "#937860"), ("TGA pair shared with another species", "#64B5CD"),
           ("nothing found upstream", "#4D4D4D")]


def tsv(path): return list(csv.DictReader(open(path), delimiter="\t"))


def style():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"], "font.size": 7, "axes.titlesize": 7.5,
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False})


def save(fig, out):
    for ext in ("pdf", "svg", "png"): fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    print("written", out)


def model_figure(a):
    pos = tsv(f"{a.tables}/{a.prefix}_positions.tsv")
    pairs = tsv(f"{a.tables}/{a.prefix}_pairs.tsv")
    seeds = tsv(f"{a.tables}/{a.prefix}_seeds.tsv")
    summ = {r["measure"]: r["value"] for r in tsv(f"{a.tables}/{a.prefix}_summary.tsv")}
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.45, wspace=0.42)
    # A. structure
    ax = fig.add_subplot(gs[0, :])
    n = len(pos)
    sp = [r for r in pairs if r["set"] == "seeds"]
    cmap = plt.get_cmap("viridis")
    for r in sp:
        i, j, f = int(r["i"]), int(r["j"]), float(r["canonical_fraction"])
        c, w = (i + j) / 2, j - i
        ax.add_patch(Arc((c, 0.9), w, w * 0.32, theta1=0, theta2=180, color=cmap((f - 0.5) / 0.5), lw=1.4))
    for r in pos:
        k, frac = int(r["position"]), float(r["most_frequent_fraction"])
        shade = str(round(0.75 - 0.75 * max(0.0, (frac - 0.4) / 0.6), 2))
        ax.text(k, 0.25, r["consensus"].upper(), ha="center", va="center", fontsize=6.3, family="monospace",
                color=shade, fontweight="bold" if frac >= 0.9 else "normal")
    # helices and motifs
    hel = collections.defaultdict(list)
    for r in sp: hel[r["helix"]].append((int(r["i"]), int(r["j"])))
    for h, ps in hel.items():
        lo5, hi5 = min(p[0] for p in ps), max(p[0] for p in ps)
        lo3, hi3 = min(p[1] for p in ps), max(p[1] for p in ps)
        for x0, x1 in ((lo5, hi5), (lo3, hi3)):
            ax.plot([x0 - 0.4, x1 + 0.4], [-0.25, -0.25], color=HELIX_COLOURS[h], lw=2, solid_capstyle="butt")
        ax.text((lo5 + hi5) / 2, -0.75, f"helix {h}", ha="center", va="center", fontsize=6.3, color=HELIX_COLOURS[h])
    core = [(11, 13, "UGA"), (54, 55, "GA"), (26, 28, "AAA")]
    for x0, x1, lab in core:
        ax.plot([x0 - 0.4, x1 + 0.4], [-0.25, -0.25], color="black", lw=2, solid_capstyle="butt")
        ax.text((x0 + x1) / 2, -0.75, lab, ha="center", va="center", fontsize=6.3)
    ax.set_xlim(0, n + 1); ax.set_ylim(-1.1, n * 0.16 + 1.2)
    ax.set_aspect("auto"); ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 1.0))
    cb = fig.colorbar(sm, ax=ax, fraction=0.015, pad=0.01, aspect=12)
    cb.set_label("seeds with a canonical pair", fontsize=6.3); cb.ax.tick_params(labelsize=6)
    ax.set_title(f"A  consensus of the {summ['seeds']} seeds (66 positions; bold: residue in at least 90% of the seeds)",
                 loc="left", fontweight="bold")
    # B. identities
    ax = fig.add_subplot(gs[1, 0])
    ident = tsv(f"{a.tables}/{a.prefix}_identities.tsv")
    vals = [float(r["identity"]) for r in ident]
    orth = [float(r["identity"]) for r in ident if r["orthologous"] == "yes"]
    bins = [x / 40 for x in range(8, 41)]
    ax.hist(vals, bins=bins, color="#BBBBBB", label=f"all pairs of seeds ({len(vals):,})")
    ax.hist(orth, bins=bins, color="#C44E52", label=f"orthologs, O. dioica ({len(orth)})")
    ax.set_yscale("log")
    ax.set_xlabel("identity (consensus positions)"); ax.set_ylabel("pairs of seeds")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("B  diversity of the seeds", loc="left", fontweight="bold")
    # C. seeds against decoys
    ax = fig.add_subplot(gs[1, 1])
    dp = {(r["i"], r["j"]): float(r["canonical_fraction"]) for r in pairs if r["set"] == "decoys"}
    for r in sp:
        ax.scatter(dp[(r["i"], r["j"])], float(r["canonical_fraction"]), s=12, color=HELIX_COLOURS[r["helix"]], zorder=3)
    ax.plot([0.5, 1], [0.5, 1], color="#999999", lw=0.6, ls="--")
    for h, c in HELIX_COLOURS.items(): ax.scatter([], [], s=12, color=c, label=f"helix {h}")
    ax.set_xlim(0.5, 1.02); ax.set_ylim(0.5, 1.02)
    ax.set_xlabel(f"decoys ({summ['decoys (best hits in shuffled windows)']}; {summ['decoy scores (bits)']} bits)")
    ax.set_ylabel("seeds")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("C  canonical pairs per base pair", loc="left", fontweight="bold")
    # D. orthologs
    ax = fig.add_subplot(gs[1, 2])
    bars = []
    for kind, lab in (("single", "single\nchange"), ("double", "double\nchange")):
        k, nn = map(int, summ[f"orthologs: {kind} changes from a canonical pair, kept canonical"].split("/"))
        e = float(summ[f"orthologs: {kind} changes, expected fraction kept if random"])
        p = summ[f"orthologs: {kind} changes, binomial p (kept more than expected)"]
        bars.append((lab, k / nn, e, f"{k}/{nn}", p))
    x = range(len(bars))
    ax.bar([i - 0.18 for i in x], [b[1] for b in bars], width=0.36, color="#C44E52", label="observed")
    ax.bar([i + 0.18 for i in x], [b[2] for b in bars], width=0.36, color="#BBBBBB", label="random new base(s)")
    for i, b in zip(x, bars):
        ax.text(i - 0.18, b[1] + 0.02, b[3], ha="center", va="bottom", fontsize=6)
        ax.text(i, 1.04, f"p = {b[4]}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(list(x)); ax.set_xticklabels([b[0] for b in bars])
    ax.set_ylim(0, 1.45); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("changes from a canonical pair\nthat keep a canonical pair")
    ax.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02), handlelength=1, columnspacing=0.8)
    rates = []
    for cat, lab in (("paired", "paired"), ("core", "core motifs"), ("other unpaired", "other unpaired")):
        s, c = int(summ[f"orthologs: {cat} sites"]), int(summ[f"orthologs: {cat} sites with a change"])
        rates.append(f"{lab} {100 * c / s:.0f}%")
    ax.text(0.5, -0.36, "sites that differ: " + ", ".join(rates), transform=ax.transAxes, ha="center", fontsize=6)
    ax.set_title(f"D  orthologs ({summ['orthologs: genes compared (three assemblies)']} genes, 3 assemblies)",
                 loc="left", fontweight="bold")
    save(fig, f"{a.figures}/secis_model")


def classify(a):
    hits = tsv(f"{a.tables}/stage2_secis_scan_hits.tsv")
    scr = {r["id"]: r for r in tsv(f"{a.tables}/stage2_uga_screen.tsv")}
    loci = [r for r in tsv(f"{a.tables}/stage2_candidates_genomes.tsv") if r["candidate"].startswith("KKSU") and r["stop"]]
    out = []
    for h in hits:
        if h["at the genome threshold"] != "yes": continue
        up = h["known gene upstream (distance from its stop codon)"]
        cands = re.findall(r"cand\d+", h["candidates of the comparative screen"])
        st = [scr[c] for c in cands if c in scr]
        frm, to = int(h["from"]), int(h["to"])
        kksu = any(l["genome"] == h["genome"] and l["contig"] == h["contig"] and l["strand"] == h["strand"] and
                   0 <= ((frm - int(l["position"])) if h["strand"] == "+" else (int(l["position"]) - frm)) <= 3000 for l in loci)
        if up.startswith("gold standard KKSU") or kksu: c = "KKSU"
        elif up.startswith("gold standard") or "(selenocysteine" in up or any(r["status"] == "known" for r in st):
            c = "selenoprotein gene"
        elif any(r["status"] == "homolog of a known selenoprotein" for r in st): c = "copy of a known selenoprotein"
        elif up.startswith("stage 1"): c = "stage 1 prediction without selenocysteine"
        elif any(r["status"] == "candidate" and r["other_species"] for r in st): c = "TGA pair shared with another species"
        else: c = "nothing found upstream"
        out.append({"genome": h["genome"], "contig": h["contig"], "strand": h["strand"], "from": frm, "to": to,
                    "bits": h["bits"], "class": c, "upstream": up, "candidates": h["candidates of the comparative screen"]})
    with open(f"{a.tables}/stage2_secis_scan_classes.tsv", "w") as o:
        cols = list(out[0])
        o.write("\t".join(cols) + "\n")
        for r in out: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    return out


def scores_figure(a):
    summ = tsv(f"{a.tables}/stage2_secis_scan_summary.tsv")
    classes = classify(a)
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.5, 0.9], hspace=0.6, wspace=0.35)
    bins = [15, 20, 25, 30, 35, 40, 50]
    rows = {(r["genome"], r["set"]): r for r in summ if not r["genome"].startswith("#")}
    thr = {r["genome"]: r["genome threshold (bits)"] for r in summ if r["set"] == "genome"}
    colours = plt.get_cmap("tab10")
    for k, (title, gens) in enumerate(GROUPS):
        ax = fig.add_subplot(gs[0, k])
        for c, g in enumerate(gens):
            mb = float(rows[(g, "genome")]["Mb (both strands)"])
            for st, ls in (("genome", "-"), ("shuffled", "--")):
                y = [int(rows[(g, st)][f">= {b} bits"]) / mb for b in bins]
                xs = [b for b, v in zip(bins, y) if v > 0]
                ax.plot(xs, [v for v in y if v > 0], ls=ls, color=colours(c), lw=1,
                        label=NAMES[g] if st == "genome" else None)
        ax.set_yscale("log"); ax.set_ylim(5e-4, 50)
        ax.set_xlabel("score (bits)")
        if k == 0: ax.set_ylabel("hits per Mb at or above")
        ax.set_title(title, loc="left", fontsize=7)
        ax.legend(frameon=False, fontsize=5.5, loc="upper right", handlelength=1.2)
    fig.text(0.07, 0.925, "A  genome (solid) and window-shuffled copy (dashed)", fontweight="bold", fontsize=7.5)
    # B. classes
    ax = fig.add_subplot(gs[1, :])
    order = [g for _, gens in GROUPS for g in gens]
    col = dict(CLASSES)
    for y, g in enumerate(order):
        xs = [r for r in classes if r["genome"] == g]
        for r in xs:
            jitter = (zlib.crc32(f"{r['contig']}{r['from']}".encode()) % 7 - 3) * 0.06
            ax.scatter(float(r["bits"]), -y + jitter, s=9, color=col[r["class"]], zorder=3, lw=0)
        if thr.get(g) not in (None, "-"):
            ax.plot([float(thr[g])] * 2, [-y - 0.35, -y + 0.35], color="#999999", lw=0.8)
    ax.set_yticks([-y for y in range(len(order))]); ax.set_yticklabels([NAMES[g] for g in order], style="italic")
    ax.set_xlabel("score (bits)")
    for lab, c in CLASSES:
        nn = sum(r["class"] == lab for r in classes)
        if nn: ax.scatter([], [], s=12, color=c, label=f"{lab} ({nn})")
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    ax.set_title("B  hits at the genome threshold (grey tick), by what lies upstream on the same strand", loc="left",
                 fontweight="bold")
    # C. controls
    ax = fig.add_subplot(gs[2, :2])
    scr = tsv(f"{a.tables}/stage2_uga_screen.tsv")
    groups = []
    for g in ("Hsap", "Cint"):
        for model in ("oikopleura", "RF00031"):
            seqs = collections.defaultdict(lambda: [0.0, False])
            for r in scr:
                if r["genome"] == g and r["model"] == model:
                    s = seqs[r["sequence"]]; s[0] = max(s[0], float(r["bits"])); s[1] = s[1] or r["status"] == "known"
            groups.append((f"{NAMES[g]}, {'Oikopleura' if model == 'oikopleura' else 'RF00031'}", seqs))
    for y, (lab, seqs) in enumerate(groups):
        for b, kn in seqs.values():
            jitter = (zlib.crc32(f"{lab}{b}{kn}".encode()) % 7 - 3) * 0.05
            ax.scatter(b, -y + jitter, s=8, lw=0, color="#C44E52" if kn else "#4D4D4D", zorder=3)
    ax.set_yticks([-y for y in range(len(groups))]); ax.set_yticklabels([g[0] for g in groups])
    ax.scatter([], [], s=10, color="#C44E52", label="selenoprotein mRNAs"); ax.scatter([], [], s=10, color="#4D4D4D", label="other mRNAs")
    ax.legend(frameon=False, fontsize=6, loc="upper right", bbox_to_anchor=(1.45, 1.05))
    ax.set_xlabel("best score per mRNA (bits)")
    ax.set_title("C  controls: RefSeq mRNAs with a TGA frame upstream of a hit of 20 bits or more", loc="left", fontweight="bold")
    save(fig, f"{a.figures}/secis_scores")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--prefix", default="stage2_secis_v2_model")
    ap.add_argument("--figures", default="figures")
    a = ap.parse_args()
    style()
    model_figure(a)
    scores_figure(a)


if __name__ == "__main__":
    main()
