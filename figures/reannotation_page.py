#!/usr/bin/env python3
"""A4 portrait page of reannotated selenoprotein genes: three genes per row, as many rows as needed.

Each gene has two parts, drawn as in figures/reannotation_figure.py:
- gene models (top): in each O. dioica assembly, the public models (OKI2018 RefSeq, OSKA2016 ANISEED, OdB3 GenBank,
  grey) above the gold-standard model of step 2 (blue), in gene orientation and aligned on the selenocysteine TGA
  (x: nucleotides from the TGA, one scale per gene); introns coloured by splice site without changing their drawn
  length; SECIS elements: best hit of the current appendicularian model (cmsearch) after the stop codon. Public
  models are named, with the curator's note, in a caption above them (letters when several share a row);
  exons more than the margin outside the gold-standard transcript are counted, not drawn.
- proteins (bottom): reference proteins and the gold-standard proteins of the three assemblies on the coordinates
  of the first reference (global alignments, BLOSUM62, free end gaps), with its Pfam domains and the residues
  aligned to its selenocysteines.
The assembly names and the labels of the gold-standard proteins are written once per row of genes, on the left;
the key is at the bottom of the page. Sizes are physical (mm), so that rows of genes keep the same density.
Inputs: figures/inputs/reannotation_<gene>.{tsv,fa} (figures/reannotation_inputs.py) and the gold-standard outputs
of step 2. Output: <out>/<name>.{pdf,svg,png}.
"""
import argparse, collections, csv, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.textpath import TextToPath
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reannotation_figure as rf  # noqa: E402

plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 6, "hatch.linewidth": 0.4, "xtick.labelsize": 5, "xtick.major.size": 1.5,
                     "xtick.major.width": 0.4, "axes.linewidth": 0.4})
GOLD, PUBLIC, OTHER, SEC, SECIS, DOM = rf.GOLD, rf.PUBLIC, rf.OTHER, rf.SEC, rf.SECIS, rf.DOM
GENOMES, NAMES, SITE_COLOURS = rf.GENOMES, rf.NAMES, rf.SITE_COLOURS
MM = 1 / 25.4
PAGE_W, PAGE_H = 210, 297                  # mm
LEFT, RIGHT, COL_GAP, TOP, BOTTOM = 24, 5, 7, 8, 6
CDS_H, UTR_H, CHEV = 1.4, 0.55, 0.8        # gene models: coding and UTR exons, intron chevrons (half of the one-gene figure)
BAR_H = 1.4                                # proteins
CAPTION = 5                                # pt
LH = 2.1                                   # mm per caption line
_t2p = TextToPath()
BOX = dict(facecolor="white", edgecolor="none", pad=0.3)


def width_mm(s, size, weight="normal"):
    w, _, _ = _t2p.get_text_width_height_descent(s, FontProperties(family="DejaVu Sans", size=size, weight=weight), ismath=False)
    return w * 25.4 / 72


def wrap(text, width, size):
    """Lines of text no wider than width (mm) at the font size (pt), cut between words."""
    lines, cur = [], ""
    for w in text.split():
        t = f"{cur} {w}" if cur else w
        if cur and width_mm(t, size) > width: lines.append(cur); cur = w
        else: cur = t
    return lines + ([cur] if cur else [])


def short(name): return name.replace("Oidioi.CG.OSKA2016.", "")


def load(gene, a):
    rows = list(csv.DictReader(open(f"{a.inputs}/reannotation_{gene}.tsv"), delimiter="\t"))
    d = {"gene": gene, "seqs": rf.read_fasta(f"{a.inputs}/reannotation_{gene}.fa"), "public": collections.defaultdict(dict),
         "gold_notes": {}, "gold_sites": {}, "refs": collections.OrderedDict()}
    for r in rows:
        if r["kind"] == "public": d["public"][r["assembly"]].setdefault(r["name"], {})[r["field"]] = r["value"]
        elif r["kind"] == "gold": (d["gold_sites"] if r["field"] == "sites" else d["gold_notes"])[r["assembly"]] = r["value"]
        elif r["kind"] == "reference":
            x = d["refs"].setdefault(r["name"], {"domains": []})
            if r["field"] == "domain": x["domains"].append(r["value"])
            else: x[r["field"]] = r["value"]
    d["gold"] = {g: rf.gold_model(f"{a.gold}/{g}.gff3", gene) for g in GENOMES}
    d["gprot"] = {}
    for g in GENOMES:
        name = None
        for l in open(f"{a.gold}/{g}.proteins.faa"):
            if l.startswith(">"): name = l[1:].split("|")[0]
            elif name == gene: d["gprot"][g] = l.strip()
    d["secis"] = rf.secis_hits(a.cm, a.gold, gene)
    d["x"] = {g: ((lambda p, s0=sec[0]: p - s0) if st == "+" else (lambda p, s1=sec[1]: s1 - p))
              for g, (st, ex, cds, sec) in d["gold"].items() if ex}
    d["span"] = {g: (min(min(d["x"][g](p), d["x"][g](q)) for p, q in d["gold"][g][1]),
                     max(max(d["x"][g](p), d["x"][g](q)) + 1 for p, q in d["gold"][g][1])) for g in d["x"]}
    longest = max(b - a_ for a_, b in d["span"].values())
    d["margin"] = a.margin or int(min(400, max(150, 0.15 * longest)))
    d["X0"] = min(v[0] for v in d["span"].values()) - d["margin"]
    d["X1"] = max(v[1] for v in d["span"].values()) + d["margin"]
    for g in GENOMES:   # public models, 5' to 3', with the letters used when several share the row
        if g not in d["x"]: continue
        x = d["x"][g]
        ms = sorted(d["public"][g].items(), key=lambda kv: min(min(x(p), x(q)) for p, q in rf.blocks(kv[1]["exons"])))
        d["public"][g] = collections.OrderedDict(ms)
    return d


def caption(d, g, col_w):
    """Caption lines of the public models of an assembly, and whether letters are drawn on the models."""
    ms = d["public"].get(g, {})
    if not ms: return wrap("no annotated model", col_w, CAPTION), False
    parts = []
    for k, (name, m) in enumerate(ms.items()):
        lead = f"{'abcdefgh'[k]} " if len(ms) > 1 else ""
        parts.append(f"{lead}{short(name)}" + (f": {m['note']}" if m.get("note") else ""))
    return wrap("; ".join(parts), col_w, CAPTION), len(ms) > 1


def draw_model(ax, y, exons, cds, x, color, hatch=False, sites=None):
    """As reannotation_figure.draw_model, with the heights of this page (mm)."""
    ex = sorted(exons)
    for e1, e2 in [(e1, e2) for e1, e2 in zip(ex, ex[1:]) if e2[0] > e1[1] + 1]:
        r1, r2 = sorted([(min(x(e1[0]), x(e1[1])), max(x(e1[0]), x(e1[1])) + 1), (min(x(e2[0]), x(e2[1])), max(x(e2[0]), x(e2[1])) + 1)])
        c = rf.site_colour((sites or {}).get((e1[1], e2[0]), ""), color)
        ax.plot([r1[1], (r1[1] + r2[0]) / 2, r2[0]], [y, y + CHEV, y], color=c, lw=0.4 if c == color else 0.9,
                solid_joinstyle="miter", zorder=3 if c != color else 2)
    for a, b in ex:
        a0, a1 = min(x(a), x(b)), max(x(a), x(b)) + 1
        ax.add_patch(Rectangle((a0, y - UTR_H / 2), a1 - a0, UTR_H, facecolor=color, edgecolor="none"))
    for a, b in cds:
        a0, a1 = min(x(a), x(b)), max(x(a), x(b)) + 1
        ax.add_patch(Rectangle((a0, y - CDS_H / 2), a1 - a0, CDS_H, facecolor="white" if hatch else color,
                               edgecolor=color, lw=0.4, hatch="//////" if hatch else None))


def stem_loop(ax, x0, y, w, color):
    """SECIS glyph (w: width in data units) standing on a model drawn at y."""
    for xx in (x0, x0 + w): ax.plot([xx, xx], [y + 0.5, y + 1.9], color=color, lw=0.6)
    ax.add_patch(Ellipse((x0 + w / 2, y + 2.35), w * 1.7, 0.9, facecolor="none", edgecolor=color, lw=0.6))


def layout_models(genes, col_w):
    """y (mm, 0 at the top, downwards) of the public and gold-standard tracks of each assembly for one row of genes,
    with the caption lines of each gene, and the height of the part."""
    y, tracks = -3.6, {}
    for g in GENOMES:
        caps = {d["gene"]: caption(d, g, col_w) for d in genes}
        gnotes = {d["gene"]: wrap(d["gold_notes"].get(g, ""), col_w, CAPTION) for d in genes}
        n_cap = max(len(c[0]) for c in caps.values())
        letters = any(c[1] for c in caps.values())
        top = y
        y -= n_cap * LH + 1.6 + (1.9 if letters else 0)
        pub = y
        y -= 5.2
        gold = y
        y -= 1.4
        tracks[g] = {"caption_top": top, "pub": pub, "gold": gold, "gold_note_top": y, "caps": caps, "gnotes": gnotes}
        y -= max(len(v) for v in gnotes.values()) * LH + 2.0
    return tracks, -y


def draw_models(ax, d, tracks, col_w):
    X0, X1 = d["X0"], d["X1"]
    ntmm = (X1 - X0) / col_w                      # nt per mm
    last_gold = None
    for g in GENOMES:
        t = tracks[g]
        for k, line in enumerate(t["caps"][d["gene"]][0]):
            ax.text(X0, t["caption_top"] - k * LH, line, fontsize=CAPTION, color="0.35", va="top", ha="left", bbox=BOX, zorder=5)
        if g not in d["x"]: continue
        st, ex, cds, sec = d["gold"][g]
        x = d["x"][g]
        lo, hi = d["span"][g][0] - d["margin"], d["span"][g][1] + d["margin"]
        ms = d["public"][g]
        for k, (name, m) in enumerate(ms.items()):
            e, c = rf.blocks(m["exons"]), rf.blocks(m.get("cds", ""))
            inside = lambda b: lo <= min(x(b[0]), x(b[1])) and max(x(b[0]), x(b[1])) <= hi
            e_in, c_in = [b for b in e if inside(b)], [b for b in c if inside(b)]
            hatch = m.get("frame") == "other"
            col = OTHER if hatch else PUBLIC
            draw_model(ax, t["pub"], e_in, c_in, x, col, hatch=hatch, sites=rf.site_map(e, m.get("sites", "")))
            if not e_in: continue
            xm = min(min(x(p), x(q)) for p, q in e_in); xM = max(max(x(p), x(q)) + 1 for p, q in e_in)
            n_left = sum(max(x(p), x(q)) < lo for p, q in e); n_right = sum(min(x(p), x(q)) > hi for p, q in e)
            if n_left: ax.text(xm, t["pub"] - 1.3, f"◀ {n_left} more exons", fontsize=4.5, color=PUBLIC, va="top")
            if n_right: ax.text(xM, t["pub"] - 1.3, f"{n_right} more exons ▶", fontsize=4.5, color=PUBLIC, va="top", ha="right")
            if len(ms) > 1:
                ax.text((xm + xM) / 2, t["pub"] + 1.2, "abcdefgh"[k], fontsize=5, fontweight="bold", color=col, ha="center", va="bottom")
        draw_model(ax, t["gold"], ex, cds, x, GOLD, sites=rf.site_map(ex, d["gold_sites"].get(g, "")))
        for k, line in enumerate(t["gnotes"][d["gene"]]):
            ax.text(X0, t["gold_note_top"] - k * LH, line, fontsize=CAPTION, color=GOLD, va="top", ha="left", bbox=BOX, zorder=5)
        if g in d["secis"]:
            f0, f1, bits = d["secis"][g]
            stop = max(q for _, q in cds) if st == "+" else min(p for p, _ in cds)
            w = 1.1 * ntmm
            x0 = x(stop) + f0
            stem_loop(ax, x0, t["gold"], w, SECIS)
            label = f"{bits:.0f} bits"
            if x0 + w + 1.2 * ntmm + width_mm(label, 4.5) * ntmm < X1:
                ax.text(x0 + w + 1.2 * ntmm, t["gold"] + 2.2, label, fontsize=4.5, color=SECIS, va="center")
            else:
                ax.text(x0 - 1.2 * ntmm, t["gold"] + 2.2, label, fontsize=4.5, color=SECIS, va="center", ha="right")
        last_gold = t["gold"]
    ax.plot([1.5, 1.5], [last_gold - 1.5, -2.6], color=SEC, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
    ax.text(1.5, -2.4, "TGA", color=SEC, ha="center", fontsize=5, fontweight="bold", va="bottom")
    ax.set_xlim(X0, X1)
    step = next(s for s in (250, 500, 1000, 2000) if (X1 - X0) / s <= 6)
    ax.set_xticks(range(math.ceil(X0 / step) * step, int(X1) + 1, step))
    ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ("left", "right", "top")]
    ax.tick_params(axis="x", pad=1)


def layout_proteins(genes):
    """y (mm, downwards) of the label lines and bars of the references (as many rows as the gene with most) and of
    the gold-standard proteins, and the height of the part."""
    nref = max(len(d["refs"]) for d in genes)
    y = -2 * LH - 0.6                 # domain and selenocysteine labels (two lines when they overlap)
    refs = []
    for _ in range(nref):
        refs.append((y, y - LH - 0.9)); y -= LH + 0.9 + 2.6
    gold = []
    for _ in GENOMES:
        gold.append(y); y -= 2.3
    return refs, gold, -(y + 2.3) + 1.6


def draw_proteins(ax, d, refs_y, gold_y, col_w, fig):
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    al.end_insertion_score = al.end_deletion_score = 0
    axis_label = next(k for k, v in d["refs"].items() if v.get("axis") == "yes")
    axis = d["seqs"][axis_label]
    L = len(axis) * 1.13
    ammm = L / col_w
    secpos = [int(p) for p in d["refs"][axis_label].get("selenocysteine", "").split(",") if p]
    doms = []
    for s in d["refs"][axis_label]["domains"]:
        acc, span, name = s.split(":", 2)
        d0, d1 = map(int, span.split("-"))
        doms.append((d0, d1, name))
    rows = []
    for k, (lbl, v) in enumerate(d["refs"].items()):
        acc = v.get("accession", "").replace("UniProt ", "")
        text = f"{lbl} {acc}" + (f"; {v['note']}" if v.get("note") else "")
        rows.append((refs_y[k][1], d["seqs"][lbl], "black", (refs_y[k][0], text)))
    rows += [(gold_y[k], d["gprot"][g], GOLD, None) for k, g in enumerate(GENOMES) if g in d["gprot"]]
    for yy, p, col, label in rows:
        m = rf.protein_map(p, axis, al)
        h0, h1 = min(m.values()), max(m.values())
        ax.add_patch(Rectangle((h0, yy - BAR_H / 2), h1 - h0, BAR_H, facecolor="white", edgecolor=col, lw=0.5))
        for d0, d1, _ in doms:
            a0, a1 = max(d0, h0), min(d1, h1)
            if a1 > a0: ax.add_patch(Rectangle((a0, yy - BAR_H / 2), a1 - a0, BAR_H, facecolor=DOM, edgecolor=col, lw=0.5))
        for sp in secpos:
            site = next((p[i] for i, j in m.items() if j == sp), None)
            if site:
                ax.plot([sp], [yy], "o", ms=4.2, color=SEC if site == "U" else "white", mec=SEC if site == "U" else "black", mew=0.5)
                ax.text(sp, yy - 0.05, site, ha="center", va="center", fontsize=4, color="white" if site == "U" else "black", fontweight="bold")
        ax.text(h1 + 1.0 * ammm, yy, f"{len(p)} aa", va="center", fontsize=4.5, color="0.3")
        if label: ax.text(0, label[0], label[1], fontsize=CAPTION, va="top", ha="left", color="0.15")
    # domain and selenocysteine labels; a label that would overlap the previous one goes one line up
    items = sorted([((d0 + d1) / 2, name, "0.25", "normal") for d0, d1, name in doms] +
                   [(sp, f"U{sp}", SEC, "bold") for sp in secpos])
    placed = []
    for xc, text, col, wt in items:
        half = width_mm(text, CAPTION, wt) / 2 * ammm
        line = 0
        while any(l_ == line and not (xc + half < a_ or xc - half > b_) for a_, b_, l_ in placed): line += 1
        placed.append((xc - half - 0.8 * ammm, xc + half + 0.8 * ammm, line))
        ax.text(xc, -0.4 - (1 - min(line, 1)) * LH, text, ha="center", va="top", fontsize=CAPTION, color=col, fontweight=wt)
    ax.set_xlim(0, L)
    step = next(s for s in (50, 100, 200) if L / s <= 7)
    ax.set_xticks(range(0, len(axis) + 1, step))
    ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ("left", "right", "top")]
    ax.tick_params(axis="x", pad=1)
    ax.set_xlabel(f"residue of {axis_label} (global alignments)", fontsize=5, labelpad=1)


def key(fig, used_sites):
    """Key at the bottom of the page: items one after the other, a new line when the page width is reached; one
    group for the gene models, one for the proteins."""
    W = PAGE_W - LEFT - RIGHT
    G = 8.5                                       # width of a glyph, mm

    def exon(ax, cx, y, col="0.4"):
        ax.add_patch(Rectangle((cx, y - UTR_H / 2), 7, UTR_H, facecolor=col, edgecolor="none"))
        ax.add_patch(Rectangle((cx + 2.5, y - CDS_H / 2), 4.5, CDS_H, facecolor=col, edgecolor=col, lw=0.4))

    def chev(c, lw):
        return lambda ax, cx, y: ax.plot([cx + 1, cx + 4, cx + 7], [y - 0.4, y + 0.6, y - 0.4], color=c, lw=lw, solid_joinstyle="miter")

    def circle(letter, fill, edge, tc):
        def f(ax, cx, y):
            ax.plot([cx + 4], [y], "o", ms=4.2, color=fill, mec=edge, mew=0.5)
            ax.text(cx + 4, y - 0.05, letter, ha="center", va="center", fontsize=4, color=tc, fontweight="bold")
        return f
    models = [(exon, "exon: coding (thick), UTR (thin)", "black"),
              (lambda ax, cx, y: exon(ax, cx, y, PUBLIC), "public model", "black"),
              (lambda ax, cx, y: exon(ax, cx, y, GOLD), "gold-standard model", "black"),
              (lambda ax, cx, y: ax.add_patch(Rectangle((cx + 1, y - CDS_H / 2), 6, CDS_H, facecolor="white", edgecolor=OTHER,
                                                        lw=0.4, hatch="//////")), "coding sequence in another frame", "black"),
              (lambda ax, cx, y: ax.plot([cx + 4, cx + 4], [y - 1.5, y + 1.5], color=SEC, lw=0.5, ls=(0, (2, 1.5))),
               "selenocysteine TGA (x axis: nucleotides from it, gene orientation)", "black"),
              (lambda ax, cx, y: stem_loop(ax, cx + 3.4, y - 1.4, 1.1, SECIS), "SECIS element, bit score (current model)", "black"),
              (lambda ax, cx, y: ax.text(cx + 7, y, "◀ n", fontsize=4.5, color=PUBLIC, va="center", ha="right"),
               "exons beyond the drawn region", "black")]
    sites = ["GT-AG"] + [k for k in SITE_COLOURS if k in used_sites] + \
            (["other"] if any(v not in SITE_COLOURS and v not in ("GT-AG", "?", "") for v in used_sites) else [])
    introns = [(None, "intron splice sites (drawn length unchanged):", "black")]
    introns += [(chev("0.4" if k == "GT-AG" else SITE_COLOURS.get(k, "black"), 0.4 if k == "GT-AG" else 0.9),
                {"GT-AG": "GT-AG (colour of the model)", "other": "other site"}.get(k, k),
                "0.4" if k == "GT-AG" else SITE_COLOURS.get(k, "black")) for k in sites]
    proteins = [(lambda ax, cx, y: ax.add_patch(Rectangle((cx + 1, y - BAR_H / 2), 6, BAR_H, facecolor=DOM, edgecolor="black", lw=0.5)),
                 "Pfam domain of the first reference", "black"),
                (circle("U", SEC, SEC, "white"), "selenocysteine aligned to one of the first reference", "black"),
                (circle("C", "white", "black", "black"), "other residue aligned to it", "black")]
    lines, line, cx = [], [], 0
    for group in (models[:4], models[4:], introns, proteins):
        for f, text, col in group:
            w = (G if f else 0) + width_mm(text, 5.5) + 3.5
            if line and cx + w > W: lines.append(line); line, cx = [], 0
            line.append((cx, f, text, col)); cx += w
        lines.append(line); line, cx = [], 0
    h = 3.4 * len(lines)
    ax = fig.add_axes([LEFT / PAGE_W, BOTTOM / PAGE_H, W / PAGE_W, h / PAGE_H])
    ax.set_xlim(0, W); ax.set_ylim(-h, 0); ax.axis("off")
    for k, line in enumerate(lines):
        y = -1.7 - 3.4 * k
        for cx, f, text, col in line:
            if f: f(ax, cx, y)
            ax.text(cx + (G if f else 0), y, text, fontsize=5.5, va="center", color=col)
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", required=True, help="genes in reading order, three per row")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--inputs", default="figures/inputs")
    ap.add_argument("--cm", default="results/models/appendicularia_secis.cm")
    ap.add_argument("--margin", type=int, default=0, help="nt drawn beyond the gold-standard transcripts (0: 15%% of the longest, 150 to 400)")
    ap.add_argument("--name", default="reannotation_page")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    genes = [load(g, a) for g in a.genes]
    col_w = (PAGE_W - LEFT - RIGHT - 2 * COL_GAP) / 3
    fig = plt.figure(figsize=(PAGE_W * MM, PAGE_H * MM))
    fx = lambda mm: mm / PAGE_W
    fy = lambda mm: mm / PAGE_H                    # from the bottom of the page
    top = PAGE_H - TOP
    used = set()
    for d in genes:
        for g in GENOMES:
            used |= set(d["gold_sites"].get(g, "").split(";"))
            for m in d["public"].get(g, {}).values(): used |= set(m.get("sites", "").split(";"))
    for r in range(0, len(genes), 3):
        row = genes[r:r + 3]
        tracks, hA = layout_models(row, col_w)
        refs_y, gold_y, hB = layout_proteins(row)
        title_y = top
        axA_bottom = top - 3.5 - hA
        axB_bottom = axA_bottom - 7 - hB
        for k, d in enumerate(row):
            left = LEFT + k * (col_w + COL_GAP)
            fig.text(fx(left), fy(title_y), d["gene"], fontsize=8, fontweight="bold", va="top")
            axA = fig.add_axes([fx(left), fy(axA_bottom), fx(col_w), fy(hA)])
            axA.set_ylim(-hA, 0)
            draw_models(axA, d, tracks, col_w)
            axB = fig.add_axes([fx(left), fy(axB_bottom), fx(col_w), fy(hB)])
            axB.set_ylim(-hB, 0)
            draw_proteins(axB, d, refs_y, gold_y, col_w, fig)
            if k == 0:
                lab = blended_transform_factory(fig.transFigure, axA.transData)
                for g in GENOMES:
                    t = tracks[g]
                    db = rf.DB if hasattr(rf, "DB") else {"OKI": "RefSeq", "OSKA": "ANISEED", "OdB3": "GenBank"}
                    fig_x = fx(LEFT - 1.5)
                    axA.text(fig_x, t["pub"], db[g], transform=lab, fontsize=5.5, ha="right", va="center", color="0.25")
                    axA.text(fig_x, t["gold"], "gold standard", transform=lab, fontsize=5.5, ha="right", va="center", color=GOLD)
                    axA.text(fx(4), (t["pub"] + t["gold"]) / 2, NAMES[g], transform=lab, fontsize=7, fontweight="bold",
                             rotation=90, ha="center", va="center")
                axA.text(fx(LEFT - 1.5), -hA - 0.5, "nt from TGA", transform=lab, fontsize=5, ha="right", va="top", color="0.25")
                labB = blended_transform_factory(fig.transFigure, axB.transData)
                for g, yy in zip(GENOMES, gold_y):
                    axB.text(fx(LEFT - 1.5), yy, NAMES[g], transform=labB, fontsize=5.5, ha="right", va="center", color=GOLD)
                axB.text(fx(4), (gold_y[0] + gold_y[-1]) / 2, "gold standard", transform=labB, fontsize=6, rotation=90,
                         ha="center", va="center", color=GOLD)
        top = axB_bottom - 5.5
    key_h = key(fig, used)
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, a.name)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None)
    print("written", out, f"(lowest protein axis {top + 5.5:.0f} mm from the bottom of the page; key {BOTTOM + key_h:.0f} mm)")


if __name__ == "__main__":
    main()
