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
# class: colour (Okabe-Ito palette), marker, transparency (frequent classes only)
CLASSES = [("selenoprotein gene", "#D55E00", "o", 0.6), ("copy of a known selenoprotein", "#0072B2", "s", 1.0),
           ("KKSU", "#CC79A7", "D", 1.0), ("stage 1 prediction without selenocysteine", "#009E73", "^", 1.0),
           ("TGA pair shared with another species", "#E69F00", "v", 1.0), ("nothing found upstream", "#000000", "X", 0.6)]


def tsv(path): return list(csv.DictReader(open(path), delimiter="\t"))


def style():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"], "font.size": 7, "axes.titlesize": 7.5,
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False})


def save(fig, out):
    for ext in ("pdf", "svg", "png"): fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    print("written", out)


LOGO_COLOURS = {"A": "#009E73", "C": "#0072B2", "G": "#E69F00", "U": "#D55E00"}


def letter(ax, ch, x, y, w, h, colour):
    """Draw letter ch filling the box (x, y, w, h) in data coordinates."""
    from matplotlib.textpath import TextPath
    from matplotlib.patches import PathPatch
    from matplotlib.font_manager import FontProperties
    import matplotlib.transforms as mt
    tp = TextPath((0, 0), ch, size=1, prop=FontProperties(family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    tr = mt.Affine2D().translate(-bb.x0, -bb.y0).scale(w / bb.width, h / bb.height).translate(x, y)
    ax.add_patch(PathPatch(tr.transform_path(tp), facecolor=colour, edgecolor="none"))


def model_figure(a):
    pos = tsv(f"{a.tables}/{a.prefix}_positions.tsv")
    pairs = tsv(f"{a.tables}/{a.prefix}_pairs.tsv")
    summ = {r["measure"]: r["value"] for r in tsv(f"{a.tables}/{a.prefix}_summary.tsv")}
    n = len(pos)
    sp = [r for r in pairs if r["set"] == "seeds"]
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1], hspace=0.5, wspace=0.38)
    # A. structure, logo of the seeds, weight of each position in the model
    ga = gs[0, 0].subgridspec(4, 1, height_ratios=[1.5, 0.75, 0.55, 0.32], hspace=0.2)
    ax = fig.add_subplot(ga[0])
    cmap = plt.get_cmap("viridis")
    for r in sp:
        i, j, f = int(r["i"]), int(r["j"]), float(r["canonical_fraction"])
        ax.add_patch(Arc(((i + j) / 2, 0), j - i, (j - i) * 0.55, theta1=0, theta2=180, color=cmap((f - 0.5) / 0.5), lw=1.1))
    ax.set_xlim(0.5, n + 0.5); ax.set_ylim(0, (n - 1) * 0.55 / 2 + 1); ax.axis("off")
    cax = ax.inset_axes([0.0, 0.84, 0.16, 0.06])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 1.0)), cax=cax, orientation="horizontal")
    cb.set_ticks([0.5, 0.75, 1.0]); cb.ax.tick_params(labelsize=5, length=1.5, pad=1)
    cb.set_label("canonical pairs (seeds)", fontsize=5.5, labelpad=1); cb.ax.xaxis.set_label_position("top")
    ax.set_title(f"A  consensus structure of the {summ['seeds']} seeds", loc="left", fontweight="bold")
    # logo
    ax = fig.add_subplot(ga[1])
    for r in pos:
        k, info = int(r["position"]), float(r["information_bits"])
        cnt = {b: int(r[b]) for b in "ACGU"}
        tot = sum(cnt.values())
        y = 0.0
        for b in sorted("ACGU", key=lambda b: cnt[b]):
            h = info * cnt[b] / tot
            if h > 0.01: letter(ax, b, k - 0.43, y, 0.86, h, LOGO_COLOURS[b])
            y += h
    ax.set_xlim(0.5, n + 0.5); ax.set_ylim(0, 2)
    ax.set_yticks([0, 1, 2]); ax.set_ylabel("seeds\n(bits)", fontsize=6)
    ax.set_xticks([])
    for side in ("top", "right", "bottom"): ax.spines[side].set_visible(False)
    # weight in the model
    ax = fig.add_subplot(ga[2])
    hel = {}
    for r in sp:
        hel[int(r["i"])] = r["helix"]; hel[int(r["j"])] = r["helix"]
    core = {11: "UGA", 12: "UGA", 13: "UGA", 54: "GA", 55: "GA", 26: "AAA", 27: "AAA", 28: "AAA"}
    for r in pos:
        k = int(r["position"])
        c = HELIX_COLOURS[hel[k]] if k in hel else ("black" if k in core else "#AAAAAA")
        ax.bar(k, float(r["model_relative_entropy_bits"]), width=0.8, color=c)
    ax.set_xlim(0.5, n + 0.5); ax.set_ylim(0, 2)
    ax.set_yticks([0, 1, 2]); ax.set_ylabel("model\n(bits)", fontsize=6)
    ax.set_xticks([])
    for side in ("top", "right", "bottom"): ax.spines[side].set_visible(False)
    # helices and motifs
    ax = fig.add_subplot(ga[3])
    groups = collections.defaultdict(list)
    for r in sp: groups[r["helix"]].append((int(r["i"]), int(r["j"])))
    for h, ps in groups.items():
        for x0, x1 in ((min(p[0] for p in ps), max(p[0] for p in ps)), (min(p[1] for p in ps), max(p[1] for p in ps))):
            ax.plot([x0 - 0.4, x1 + 0.4], [0.75, 0.75], color=HELIX_COLOURS[h], lw=2, solid_capstyle="butt")
            ax.text((x0 + x1) / 2, 0.05, h, ha="center", va="center", fontsize=6, color=HELIX_COLOURS[h])
    for x0, x1, lab in ((11, 13, "UGA"), (54, 55, "GA"), (26, 28, "AAA")):
        ax.plot([x0 - 0.4, x1 + 0.4], [0.75, 0.75], color="black", lw=2, solid_capstyle="butt")
        ax.text((x0 + x1) / 2, 0.05, lab, ha="center", va="center", fontsize=6)
    ax.set_xlim(0.5, n + 0.5); ax.set_ylim(-0.4, 1); ax.axis("off")
    ax.text(0.5, -0.45, "consensus position 1-66 (helices 1, 2, 3 coloured; core motifs black)", ha="left", va="top",
            fontsize=5.5, color="#555555", transform=ax.transData)
    # D. orthologs
    ax = fig.add_subplot(gs[0, 1])
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
        s_, c = int(summ[f"orthologs: {cat} sites"]), int(summ[f"orthologs: {cat} sites with a change"])
        rates.append(f"{lab} {100 * c / s_:.0f}%")
    ax.text(0.5, -0.3, "positions that differ: " + ", ".join(rates), transform=ax.transAxes, ha="center", fontsize=6)
    ax.set_title(f"D  orthologous SECIS, {summ['orthologs: genes compared (three assemblies)']} genes, 3 assemblies",
                 loc="left", fontweight="bold")
    # B. identities
    ax = fig.add_subplot(gs[1, 0])
    ident = tsv(f"{a.tables}/{a.prefix}_identities.tsv")
    vals = [float(r["identity"]) for r in ident]
    orth = [r for r in ident if r["orthologous"] == "yes"]
    seeds_o = {r[k] for r in orth for k in ("seed_1", "seed_2")}
    genes_o = {x.split("|")[1].rsplit("/", 1)[0] for x in seeds_o}
    bins = [x / 40 for x in range(8, 41)]
    ax.hist(vals, bins=bins, color="#BBBBBB", label=f"all {summ['seeds']} seeds: {len(vals):,} pairs")
    ax.hist([float(r["identity"]) for r in orth], bins=bins, color="#C44E52",
            label=f"{len(seeds_o)} O. dioica seeds of {len(genes_o)} genes:\n{len(orth)} pairs of orthologs")
    ax.set_yscale("log")
    ax.set_xlabel("identity (consensus positions)"); ax.set_ylabel("pairs of seeds")
    ax.legend(frameon=False, loc="upper right", fontsize=6)
    ax.set_title("B  diversity of the seeds", loc="left", fontweight="bold")
    # C. seeds against decoys
    ax = fig.add_subplot(gs[1, 1])
    dp = {(r["i"], r["j"]): float(r["canonical_fraction"]) for r in pairs if r["set"] == "decoys"}
    for r in sp:
        ax.scatter(dp[(r["i"], r["j"])], float(r["canonical_fraction"]), s=12, color=HELIX_COLOURS[r["helix"]], zorder=3)
    ax.plot([0.5, 1], [0.5, 1], color="#999999", lw=0.6, ls="--")
    for h, c in HELIX_COLOURS.items(): ax.scatter([], [], s=12, color=c, label=f"helix {h}")
    ax.set_xlim(0.5, 1.02); ax.set_ylim(0.5, 1.02)
    ax.set_xlabel(f"decoys: best hits in shuffled windows\n({summ['decoys (best hits in shuffled windows)']} hits, "
                  f"{summ['decoy scores (bits)'].replace('-', ' to ')} bits)")
    ax.set_ylabel(f"seeds ({summ['seeds']})")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("C  canonical pairs at each base pair", loc="left", fontweight="bold")
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
    style_of = {lab: (c, m, al) for lab, c, m, al in CLASSES}
    for y, g in enumerate(order):
        if thr.get(g) not in (None, "-"):
            ax.plot([float(thr[g])] * 2, [-y - 0.32, -y + 0.32], color="#888888", lw=1.0, zorder=1)
        for lab, _, _, _ in CLASSES:   # frequent classes first, rare ones drawn on top
            for r in (r for r in classes if r["genome"] == g and r["class"] == lab):
                c, m, al = style_of[lab]
                jitter = (zlib.crc32(f"{r['contig']}{r['from']}".encode()) % 5 - 2) * 0.05
                ax.scatter(float(r["bits"]), -y + jitter, s=13 if m != "X" else 15, color=c, marker=m, alpha=al,
                           zorder=2 if al < 1 else 3, lw=0)
    ax.set_yticks([-y for y in range(len(order))]); ax.set_yticklabels([NAMES[g] for g in order], style="italic")
    ax.set_xlabel("score (bits)")
    for lab, c, m, al in CLASSES:
        nn = sum(r["class"] == lab for r in classes)
        if nn: ax.scatter([], [], s=16, color=c, marker=m, alpha=al, lw=0, label=f"{lab} ({nn})")
    ax.plot([], [], color="#888888", lw=1.0, marker="|", ls="none", markersize=7, label="genome threshold")
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    ax.set_title("B  hits at the genome threshold, by what lies upstream on the same strand", loc="left",
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
            jitter = (zlib.crc32(f"{lab}{b}{kn}".encode()) % 5 - 2) * 0.05
            ax.scatter(b, -y + jitter, s=13 if kn else 15, lw=0, marker="o" if kn else "X", alpha=0.7,
                       color="#D55E00" if kn else "#000000", zorder=3)
    ax.set_yticks([-y for y in range(len(groups))]); ax.set_yticklabels([g[0] for g in groups])
    ax.scatter([], [], s=13, lw=0, marker="o", alpha=0.7, color="#D55E00", label="selenoprotein mRNAs")
    ax.scatter([], [], s=15, lw=0, marker="X", alpha=0.7, color="#000000", label="other mRNAs")
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
