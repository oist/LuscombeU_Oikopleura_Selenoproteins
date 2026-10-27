#!/usr/bin/env python3
"""A4 portrait page of reannotated selenoprotein genes: three genes per row, as many rows as needed.

Each gene has three parts (thin rulers separate rows and columns of genes):
- gene models (top): in each O. dioica assembly, the public models (OKI2018 RefSeq, OSKA2016 ANISEED, OdB3 GenBank,
  grey) above the gold-standard model of step 2 (blue), in gene orientation and aligned on the selenocysteine TGA
  (x: nucleotides from the TGA, one scale per gene; stretches of more than 1,500 nt without any drawn exon, such as
  long introns, are drawn 200 nt long and marked by a break); introns coloured by splice site without changing their drawn
  length; SECIS elements: best hit of the Oikopleura SECIS model of step 3 (results/models/oikopleura_secis.cm,
  cmsearch) in the 2,000 nt after the stop codon, shortened at the next gold-standard coding sequence downstream,
  drawn as a green box with its bit score (open box below
  --secis-bits, 20 by default); letters when several public models share a row; public exons are cut at the margin
  beyond the gold-standard transcript, and those entirely beyond it are counted.
- paragraph (below the proteins): the public models of each assembly with the curator's notes, the notes on the
  gold-standard models, and the accessions of the reference proteins.
- proteins (bottom): reference proteins and the gold-standard proteins of the three assemblies on the coordinates
  of the first reference (local alignments, BLOSUM62; the gold-standard proteins of OSKA2016 and OdB3 through a
  global alignment to the first gold-standard protein; residues outside the aligned part drawn as a line as long as
  they are), with its Pfam domains and the residues aligned to its selenocysteines; right column: percentage of
  similar residues to the OdB3 protein (BLOSUM62 > 0, selenocysteine read as cysteine, over the aligned part of a
  global alignment with free end gaps, internal gaps counted as differences; none when OdB3 has no model).
The assembly names and the labels of the protein rows (one row per reference species, then the gold-standard
proteins) are written once per row of genes, on the left; the key is at the bottom of the page (with --fit, right
below the genes, and the page is cut below it). Sizes are physical (mm), so that rows of genes keep the same density.
Inputs: figures/inputs/reannotation_<gene>.{tsv,fa} (figures/reannotation_inputs.py) and the gold-standard outputs
of step 2; when figures/inputs/reannotation_<gene>.paralogs.{tsv,fa} exist (make txnl1-paralogs), the paralog of each
assembly is drawn below the gold-standard model (public model, aligned on its codon aligned to the selenocysteine)
and its protein below the gold-standard proteins. Output: <out>/<name>.{pdf,svg,png}.
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
PARALOG = "0.1"                            # paralog of the gene (public model and protein)
GENOMES, NAMES, SITE_COLOURS = rf.GENOMES, rf.NAMES, rf.SITE_COLOURS
MM = 1 / 25.4
PAGE_W, PAGE_H = 210, 297                  # mm
LEFT, RIGHT, COL_GAP, TOP, BOTTOM = 25.5, 5, 7, 8, 6
CDS_H, UTR_H, CHEV = 1.4, 0.55, 0.8        # gene models: coding and UTR exons, intron chevrons (half of the one-gene figure)
BAR_H = 1.4                                # proteins
CAPTION = 5                                # pt
NOTE, NOTE_LH = 5, 2.0                     # paragraph between gene models and proteins: font size (pt), line height (mm)
SECIS_H = 1.0                              # height of the SECIS box (mm)
SECIS_BITS = 20.0                          # SECIS hits below this score: open boxes (set by --secis-bits)
PAIR = 2.4                                 # distance between the public and the gold-standard track of an assembly (mm)
PAD = 120                                  # nt kept on each side of a shortened stretch
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


def short(name): return name.replace("Oidioi.CG.OSKA2016.", "").replace("rna-", "")


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
    d["par"], d["par_seqs"] = {}, {}
    if os.path.exists(f"{a.inputs}/reannotation_{gene}.paralogs.tsv"):
        for r in csv.DictReader(open(f"{a.inputs}/reannotation_{gene}.paralogs.tsv"), delimiter="\t"):
            d["par"].setdefault(r["assembly"], {"name": r["name"]})[r["field"]] = r["value"]
        d["par_seqs"] = rf.read_fasta(f"{a.inputs}/reannotation_{gene}.paralogs.fa")
    xp = {}
    for g, m in d["par"].items():   # nucleotides from the codon aligned to the selenocysteine
        c0, c1 = map(int, m["anchor"].split("-"))
        xp[g] = (lambda p, s0=c0: p - s0) if m["strand"] == "+" else (lambda p, s1=c1: s1 - p)
    pspan = {g: (min(xp[g](int(v)) for v in m["span"].split("-")), max(xp[g](int(v)) for v in m["span"].split("-")) + 1)
             for g, m in d["par"].items()}
    xr = {g: ((lambda p, s0=sec[0]: p - s0) if st == "+" else (lambda p, s1=sec[1]: s1 - p))
          for g, (st, ex, cds, sec) in d["gold"].items() if ex}       # nucleotides from the TGA
    span = {g: (min(min(xr[g](p), xr[g](q)) for p, q in d["gold"][g][1]),
                max(max(xr[g](p), xr[g](q)) + 1 for p, q in d["gold"][g][1])) for g in xr}
    longest = max(b - a_ for a_, b in span.values())
    d["margin"] = a.margin or int(min(400, max(150, 0.15 * longest)))
    win = {g: (v[0] - d["margin"], v[1] + d["margin"]) for g, v in span.items()}
    pwin = {g: (v[0] - d["margin"], v[1] + d["margin"]) for g, v in pspan.items()}
    X0, X1 = min(v[0] for v in list(win.values()) + list(pwin.values())), max(v[1] for v in list(win.values()) + list(pwin.values()))
    # stretches longer than --long nt without any drawn exon (a long intron) are drawn --short nt long, with a break
    feats = []
    for g in xr:
        feats += [(min(xr[g](p), xr[g](q)), max(xr[g](p), xr[g](q)) + 1) for p, q in d["gold"][g][1]]
        for m in d["public"][g].values():
            for p, q in rf.blocks(m["exons"]):
                a0, a1 = max(min(xr[g](p), xr[g](q)), win[g][0]), min(max(xr[g](p), xr[g](q)) + 1, win[g][1])
                if a1 > a0: feats.append((a0, a1))
    for g, m in d["par"].items():
        for p, q in rf.blocks(m["exons"]):
            a0, a1 = max(min(xp[g](p), xp[g](q)), pwin[g][0]), min(max(xp[g](p), xp[g](q)) + 1, pwin[g][1])
            if a1 > a0: feats.append((a0, a1))
    feats.sort()
    breaks, end = [], feats[0][1]
    for a0, a1 in feats[1:]:
        if a0 - end > a.long: breaks.append((end + PAD, a0 - PAD))
        end = max(end, a1)

    def fwd(v):
        out = v
        for b0, b1 in breaks:
            if v >= b1: out -= b1 - b0 - a.short
            elif v > b0: out -= (v - b0) * (1 - a.short / (b1 - b0))
        return out
    d["fwd"], d["raw_breaks"] = fwd, breaks
    d["breaks"] = [(fwd(b0), fwd(b1)) for b0, b1 in breaks]
    d["x"] = {g: (lambda p, f=xr[g]: fwd(f(p))) for g in xr}
    d["xp"] = {g: (lambda p, f=xp[g]: fwd(f(p))) for g in xp}
    d["pwin"] = {g: (fwd(v[0]), fwd(v[1])) for g, v in pwin.items()}
    d["xr"] = xr
    d["span"] = {g: (fwd(v[0]), fwd(v[1])) for g, v in span.items()}
    d["win"] = {g: (fwd(v[0]), fwd(v[1])) for g, v in win.items()}
    d["X0"], d["X1"], d["X0raw"], d["X1raw"] = fwd(X0), fwd(X1), X0, X1
    for g in GENOMES:   # public models, 5' to 3', with the letters used when several share the row
        if g not in d["x"]: continue
        x = d["x"][g]
        ms = sorted(d["public"][g].items(), key=lambda kv: min(min(x(p), x(q)) for p, q in rf.blocks(kv[1]["exons"])))
        d["public"][g] = collections.OrderedDict(ms)
    return d


DB = {"OKI": "RefSeq", "OSKA": "ANISEED", "OdB3": "GenBank"}


def notes(d):
    """Entries (label, text) of the paragraph of a gene: the public models of each assembly with the curator's notes
    (letters as on the models when several share a row), the gold-standard notes, and the reference proteins."""
    out = []
    for g in GENOMES:
        ms = d["public"].get(g, {})
        if not ms: text = f"no {DB[g]} model."
        elif len(ms) == 1:
            name, m = next(iter(ms.items()))
            text = f"{DB[g]} {short(name)}" + (f": {m['note']}." if m.get("note") else ".")
        else:
            text = f"{DB[g]} " + "; ".join(f"{'abcdefgh'[k]} {short(n)}" + (f": {m['note']}" if m.get("note") else "")
                                           for k, (n, m) in enumerate(ms.items())) + "."
        if g not in d["x"]:
            text = "no gold-standard model" + (f": {d['gold_notes'][g]}." if d["gold_notes"].get(g) else ".")
        elif d["gold_notes"].get(g): text += f" Gold standard: {d['gold_notes'][g]}."
        out.append((NAMES[g], text))
    if d["par"]:
        lbl = next(iter(d["par"].values())).get("label", "paralog")
        out.append((lbl, "; ".join(f"{NAMES[g]} {m.get('database', '')} {short(m['name'])}" + (f" ({m['note']})" if m.get("note") else "")
                                   for g, m in sorted(d["par"].items(), key=lambda kv: GENOMES.index(kv[0]))) + "."))
    out.append(("proteins", "; ".join(f"{lbl} {v.get('accession', '').replace('UniProt ', '')}" + (f" ({v['note']})" if v.get("note") else "")
                                      for lbl, v in d["refs"].items()) + "."))
    return out


def draw_model(ax, y, exons, cds, x, color, hatch=False, sites=None, clip=(-1e12, 1e12)):
    """As reannotation_figure.draw_model, with the heights of this page (mm); exons and coding sequences are cut at
    the window clip (plot coordinates), and introns drawn between the exons that are drawn."""
    span = lambda a, b: (max(min(x(a), x(b)), clip[0]), min(max(x(a), x(b)) + 1, clip[1]))
    ex = [e for e in sorted(exons) if span(*e)[1] > span(*e)[0]]
    drawn = []
    for e1, e2 in [(e1, e2) for e1, e2 in zip(ex, ex[1:]) if e2[0] > e1[1] + 1]:
        r1, r2 = sorted([span(*e1), span(*e2)])
        c = rf.site_colour((sites or {}).get((e1[1], e2[0]), ""), color)
        ax.plot([r1[1], (r1[1] + r2[0]) / 2, r2[0]], [y, y + CHEV, y], color=c, lw=0.4 if c == color else 0.9,
                solid_joinstyle="miter", zorder=3 if c != color else 2)
        drawn.append((r1[1], r2[0]))
    for a, b in ex:
        a0, a1 = span(a, b)
        ax.add_patch(Rectangle((a0, y - UTR_H / 2), a1 - a0, UTR_H, facecolor=color, edgecolor="none"))
    for a, b in cds:
        a0, a1 = span(a, b)
        if a1 <= a0: continue
        ax.add_patch(Rectangle((a0, y - CDS_H / 2), a1 - a0, CDS_H, facecolor="white" if hatch else color,
                               edgecolor=color, lw=0.4, hatch="//////" if hatch else None))
    return drawn


def break_mark(ax, xc, y, ntmm, color="0.3"):
    """Two slanted strokes on a white gap at (xc, y): a stretch drawn shorter than it is."""
    ax.add_patch(Rectangle((xc - 0.45 * ntmm, y - 0.9), 0.9 * ntmm, 1.8, facecolor="white", edgecolor="none", zorder=6, clip_on=False))
    for dx in (-0.3, 0.3):
        ax.plot([xc + (dx - 0.25) * ntmm, xc + (dx + 0.25) * ntmm], [y - 0.8, y + 0.8], color=color, lw=0.5, zorder=7, clip_on=False)


def mark_breaks(ax, d, y, introns, ntmm):
    """Break marks on the drawn introns (a, b: plot coordinates) that span a shortened stretch."""
    for b0, b1 in d["breaks"]:
        for a_, b_ in introns:
            if a_ <= b0 and b_ >= b1:
                xc, mid = (b0 + b1) / 2, (a_ + b_) / 2
                break_mark(ax, xc, y + CHEV * (1 - abs(xc - mid) / ((b_ - a_) / 2)), ntmm)


def secis_box(ax, x0, x1, y, strong=True):
    """SECIS element: a flat box on the model line, open when its score is below the threshold."""
    ax.add_patch(Rectangle((x0, y - SECIS_H / 2), x1 - x0, SECIS_H, facecolor=SECIS if strong else "white",
                           edgecolor=SECIS, lw=0 if strong else 0.5, zorder=4))


def layout_models(genes):
    """y (mm, 0 at the top, downwards) of the public, gold-standard and (when a gene of the row has one) paralog
    tracks of each assembly, the same for the genes of a row, and the height of the part."""
    y, tracks = -3.4, {}
    has_par = any(d["par"] for d in genes)
    for g in GENOMES:
        y -= 2.2                       # letters of the models and counts of exons beyond the window
        pub = y
        y -= PAIR                      # public and gold-standard models of one assembly: closer than two assemblies
        tracks[g] = {"pub": pub, "gold": y}
        if has_par:
            y -= PAIR + 0.6
            tracks[g]["par"] = y
        y -= 2.0
    return tracks, -y + 0.4


def draw_models(ax, d, tracks, col_w):
    X0, X1 = d["X0"], d["X1"]
    ntmm = (X1 - X0) / col_w                      # nt per mm
    last_gold = None
    for g in GENOMES:
        t = tracks[g]
        if g not in d["x"]:
            ax.text(d["X0"], t["gold"], "no gold-standard model", fontsize=4.5, color=GOLD, va="center")
            continue
        st, ex, cds, sec = d["gold"][g]
        x = d["x"][g]
        lo, hi = d["win"][g]
        ms = d["public"][g]
        for k, (name, m) in enumerate(ms.items()):
            e, c = rf.blocks(m["exons"]), rf.blocks(m.get("cds", ""))
            touches = lambda b: max(x(b[0]), x(b[1])) + 1 > lo and min(x(b[0]), x(b[1])) < hi   # exons cut at the window
            e_in, c_in = [b for b in e if touches(b)], [b for b in c if touches(b)]
            hatch = m.get("frame") == "other"
            col = OTHER if hatch else PUBLIC
            drawn = draw_model(ax, t["pub"], e_in, c_in, x, col, hatch=hatch, sites=rf.site_map(e, m.get("sites", "")), clip=(lo, hi))
            mark_breaks(ax, d, t["pub"], drawn, ntmm)
            if not e_in: continue
            xm = max(lo, min(min(x(p), x(q)) for p, q in e_in)); xM = min(hi, max(max(x(p), x(q)) + 1 for p, q in e_in))
            n_left = sum(max(x(p), x(q)) < lo for p, q in e); n_right = sum(min(x(p), x(q)) > hi for p, q in e)
            if n_left: ax.text(xm, t["pub"] + 1.0, f"◀ {n_left}", fontsize=4.5, color=PUBLIC, va="bottom")
            if n_right: ax.text(xM, t["pub"] + 1.0, f"{n_right} ▶", fontsize=4.5, color=PUBLIC, va="bottom", ha="right")
            if len(ms) > 1:   # letter over the model, after the count of exons beyond the window when there is one
                xl = (xm + xM) / 2
                if n_left: xl = max(xl, xm + (width_mm(f"◀ {n_left}", 4.5) + 1.2) * ntmm)
                if n_right: xl = min(xl, xM - (width_mm(f"{n_right} ▶", 4.5) + 1.2) * ntmm)
                ax.text(xl, t["pub"] + 1.0, "abcdefgh"[k], fontsize=5, fontweight="bold", color=col, ha="center", va="bottom")
        drawn = draw_model(ax, t["gold"], ex, cds, x, GOLD, sites=rf.site_map(ex, d["gold_sites"].get(g, "")))
        mark_breaks(ax, d, t["gold"], drawn, ntmm)
        end = d["span"][g][1]
        if g in d["secis"]:
            f0, f1, bits = d["secis"][g]
            stop = max(q for _, q in cds) if st == "+" else min(p for p, _ in cds)
            x0, x1 = d["fwd"](d["xr"][g](stop) + f0), d["fwd"](d["xr"][g](stop) + f1 + 1)
            secis_box(ax, x0, x1, t["gold"], bits >= SECIS_BITS)
            end = max(end, x1)
            ax.text(end + 0.8 * ntmm, t["gold"], f"{bits:.0f}", fontsize=4.5, color=SECIS, va="center")
        last_gold = t["gold"]
        if g in d["par"]:
            m, xpg = d["par"][g], d["xp"][g]
            plo, phi = d["pwin"][g]
            e, c = rf.blocks(m["exons"]), rf.blocks(m.get("cds", ""))
            touches = lambda b: max(xpg(b[0]), xpg(b[1])) + 1 > plo and min(xpg(b[0]), xpg(b[1])) < phi
            e_in, c_in = [b for b in e if touches(b)], [b for b in c if touches(b)]
            drawn = draw_model(ax, t["par"], e_in, c_in, xpg, PARALOG, sites=rf.site_map(e, m.get("sites", "")), clip=(plo, phi))
            mark_breaks(ax, d, t["par"], drawn, ntmm)
            if e_in:
                xm = max(plo, min(min(xpg(p), xpg(q)) for p, q in e_in)); xM = min(phi, max(max(xpg(p), xpg(q)) + 1 for p, q in e_in))
                n_left = sum(max(xpg(p), xpg(q)) < plo for p, q in e); n_right = sum(min(xpg(p), xpg(q)) > phi for p, q in e)
                if n_left: ax.text(xm, t["par"] + 1.0, f"◀ {n_left}", fontsize=4.5, color=PARALOG, va="bottom")
                if n_right: ax.text(xM, t["par"] + 1.0, f"{n_right} ▶", fontsize=4.5, color=PARALOG, va="bottom", ha="right")
            last_gold = t["par"]
    ax.plot([1.5, 1.5], [last_gold - 1.2, -2.4], color=SEC, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
    ax.text(1.5, -2.3, "TGA", color=SEC, ha="center", fontsize=5, fontweight="bold", va="bottom")
    ax.set_xlim(X0, X1)
    step = next(s for s in (250, 500, 1000, 2000) if (X1 - X0) / s <= 6)
    ticks = [v for v in range(math.ceil(d["X0raw"] / step) * step, int(d["X1raw"]) + 1, step)
             if not any(b0 - 0.6 * step < v < b1 + 0.6 * step for b0, b1 in d["raw_breaks"])
             and X0 + 0.04 * (X1 - X0) < d["fwd"](v) < X1 - 0.02 * (X1 - X0)]       # not at the very edges
    ax.set_xticks([d["fwd"](v) for v in ticks], [f"{v:g}".replace("-", "−") for v in ticks])
    bottom = ax.get_ylim()[0]
    for b0, b1 in d["breaks"]: break_mark(ax, (b0 + b1) / 2, bottom, ntmm, color="black")
    ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ("left", "right", "top")]
    ax.tick_params(axis="x", pad=1)


def layout_notes(genes, col_w):
    """Lines of the paragraph of each gene, with a hanging indent after the bold labels, and the height of the part."""
    indent = max(width_mm(lbl, NOTE, "bold") for lbl in list(NAMES.values()) + ["proteins"]) + 1.2
    out = {}
    for d in genes:
        out[d["gene"]] = [(lbl, wrap(text, col_w - indent, NOTE)) for lbl, text in notes(d)]
    n = max(sum(len(l) for _, l in v) for v in out.values())
    return out, indent, n * NOTE_LH + 0.6


def draw_notes(ax, entries, indent):
    y = -0.3
    for lbl, lines in entries:
        ax.text(0, y, lbl, fontsize=NOTE, fontweight="bold", va="top", color="0.15")
        for line in lines:
            ax.text(indent, y, line, fontsize=NOTE, va="top", color="0.15")
            y -= NOTE_LH
    ax.axis("off")


def layout_proteins(genes):
    """y (mm, downwards) of one row per reference species (the species of the genes of the row, in order) and of the
    gold-standard proteins, and the height of the part."""
    species = []
    for d in genes:
        species += [k for k in d["refs"] if k not in species]
    y = -2 * LH - 1.4                 # domain and selenocysteine labels (two lines when they overlap)
    refs = {}
    for sp in species:
        refs[sp] = y; y -= 2.3
    y -= 0.5
    gold = []
    for _ in GENOMES:
        gold.append(y); y -= 2.3
    par = []
    if any(d["par_seqs"] for d in genes):
        y -= 0.5
        for _ in GENOMES:
            par.append(y); y -= 2.3
    return refs, gold, par, -(y + 2.3) + 1.6


SIM_TO = "OdB3"   # gold-standard protein to which the similarity of the other proteins is given


def similarity(q, ref, aligner):
    """Percentage of similar residues of q to ref: residue pairs with a positive BLOSUM62 score (selenocysteine read
    as cysteine) over the aligned part of a global alignment with free end gaps (from the first to the last aligned
    pair; internal gaps count as differences)."""
    aln = aligner.align(q.replace("U", "C"), ref.replace("U", "C"))[0]
    bq, br = aln.aligned
    pairs = [(q[i], ref[j]) for (q0, q1), (r0, r1) in zip(bq, br) for i, j in zip(range(q0, q1), range(r0, r1))]
    columns = (bq[-1][1] - bq[0][0]) + (br[-1][1] - br[0][0]) - len(pairs)
    blosum = aligner.substitution_matrix
    similar = sum(a == b or blosum[a.replace("U", "C"), b.replace("U", "C")] > 0 for a, b in pairs)
    return 100 * similar / columns


def draw_proteins(ax, d, refs_y, gold_y, par_y, col_w, fig):
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    glob = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                 open_gap_score=-10, extend_gap_score=-0.5)
    glob.end_insertion_score = glob.end_deletion_score = 0
    axis_label = next((k for k, v in d["refs"].items() if v.get("axis") == "yes"), None)
    first = next(g for g in GENOMES if g in d["gprot"])
    if axis_label is None:   # no reference protein of another lineage: the first gold-standard protein is the axis
        axis = d["gprot"][first]
        secpos = [i + 1 for i, c in enumerate(axis) if c == "U"]
        ref_domains = []
        axis_name = f"{rf.NAMES[first]} gold-standard"
    else:
        axis = d["seqs"][axis_label]
        secpos = [int(p) for p in d["refs"][axis_label].get("selenocysteine", "").split(",") if p]
        ref_domains = d["refs"][axis_label]["domains"]
        axis_name = axis_label
    doms = []
    for s in ref_domains:
        acc, span, name = s.split(":", 2)
        d0, d1 = map(int, span.split("-"))
        doms.append((d0, d1, name))
    rows = []
    for lbl in d["refs"]:
        rows.append((refs_y[lbl], d["seqs"][lbl], "black", lbl, rf.protein_map(d["seqs"][lbl], axis, al)))
    # the gold-standard proteins of the other assemblies go through the first one (global alignment), so that the
    # three, nearly identical, are drawn alike where their similarity to the reference is low
    m_first = rf.protein_map(d["gprot"][first], axis, al) if axis_label else {i: i + 1 for i in range(len(axis))}
    # a selenocysteine a few residues beyond the end of the local alignment (C-terminal, as in TXNL1): the residues up to
    # it are placed after the last aligned one, so that it is drawn
    u = d["gprot"][first].find("U")
    if m_first and u > max(m_first) and u - max(m_first) <= 10:
        last = max(m_first)
        for i in range(last + 1, u + 1):
            if m_first[last] + i - last <= len(axis): m_first[i] = m_first[last] + i - last
    if not secpos:   # first reference without selenocysteine (a cysteine homolog): the residue aligned to the gold U
        secpos = [m_first[i] for i, c in enumerate(d["gprot"][first]) if c == "U" and i in m_first]
    for k, g in enumerate(GENOMES):
        if g not in d["gprot"]: continue
        m = m_first if g == first else {i: m_first[j - 1] for i, j in rf.protein_map(d["gprot"][g], d["gprot"][first], glob).items()
                                        if j - 1 in m_first}
        rows.append((gold_y[k], d["gprot"][g], GOLD, g, m))
    for k, g in enumerate(GENOMES):   # paralogs: on the axis through their own local alignment
        s = d["par_seqs"].get(f"paralog:{g}")
        if s and par_y: rows.append((par_y[k], s, PARALOG, f"paralog:{g}", rf.protein_map(s, axis, al)))
    ref = d["gprot"].get(SIM_TO)
    sim = {lbl: similarity(p, ref, glob) for _, p, _, lbl, _ in rows if ref and lbl != SIM_TO}
    ends = [(m[min(m)] - min(m), m[max(m)] + len(p) - 1 - max(m)) for _, p, _, _, m in rows]   # with unaligned residues
    lo = min(0, min(e[0] for e in ends))
    hi = max(len(axis), max(e[1] for e in ends))
    # right of the proteins: the length (after the end of each protein) and a column of similarities (mm)
    right = 1.0 + width_mm("0000 aa", 4.5) + 1.2 + width_mm("100", 4.5)
    L0, L1 = lo, hi + right * (hi - lo) / (col_w - right)
    ammm = (L1 - L0) / col_w
    for yy, p, col, label, m in rows:
        h0, h1 = min(m.values()), max(m.values())
        # residues outside the aligned part (free end gaps): a line as long as they are, at the scale of the axis
        t0, t1 = h0 - min(m), h1 + len(p) - 1 - max(m)
        ax.plot([t0, h0], [yy, yy], color=col, lw=0.6, solid_capstyle="butt")
        ax.plot([h1, t1], [yy, yy], color=col, lw=0.6, solid_capstyle="butt")
        ax.add_patch(Rectangle((h0, yy - BAR_H / 2), h1 - h0, BAR_H, facecolor="white", edgecolor=col, lw=0.5))
        for d0, d1, _ in doms:
            a0, a1 = max(d0, h0), min(d1, h1)
            if a1 > a0: ax.add_patch(Rectangle((a0, yy - BAR_H / 2), a1 - a0, BAR_H, facecolor=DOM, edgecolor=col, lw=0.5))
        for sp in secpos:
            site = next((p[i] for i, j in m.items() if j == sp), None)
            if site:
                ax.plot([sp], [yy], "o", ms=4.2, color=SEC if site == "U" else "white", mec=SEC if site == "U" else "black", mew=0.5)
                ax.text(sp, yy - 0.05, site, ha="center", va="center", fontsize=4, color="white" if site == "U" else "black", fontweight="bold")
        ax.text(t1 + 1.0 * ammm, yy, f"{len(p)} aa", va="center", fontsize=4.5, color="0.3")
        if label in sim: ax.text(L1, yy, f"{sim[label]:.0f}", va="center", ha="right", fontsize=4.5, color="black")
    # domain and selenocysteine labels; a label that would overlap the previous one goes one line up; the header of
    # the similarity column first, on the lower line
    items = sorted([((d0 + d1) / 2, name, "0.25", "normal") for d0, d1, name in doms] +
                   [(sp, f"{axis[sp - 1]}{sp}", SEC, "bold") for sp in secpos])
    placed = []
    if sim:   # header on two lines, so that it stays narrow
        for line, head in ((1, "% sim."), (0, f"to {rf.NAMES[SIM_TO]}")):
            ax.text(L1, -0.4 - (1 - line) * LH, head, ha="right", va="top", fontsize=4.5, color="black")
            placed.append((L1 - (width_mm(head, 4.5) + 0.8) * ammm, L1, line))
    for xc, text, col, wt in items:
        half = width_mm(text, CAPTION, wt) / 2 * ammm
        line = 0
        while any(l_ == line and not (xc + half < a_ or xc - half > b_) for a_, b_, l_ in placed): line += 1
        placed.append((xc - half - 0.8 * ammm, xc + half + 0.8 * ammm, line))
        ax.text(xc, -0.4 - (1 - min(line, 1)) * LH, text, ha="center", va="top", fontsize=CAPTION, color=col, fontweight=wt)
    ax.set_xlim(L0, L1)
    step = next(s for s in (50, 100, 200) if (L1 - L0) / s <= 7)
    ax.set_xticks(range(0, len(axis) + 1, step))
    ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ("left", "right", "top")]
    ax.tick_params(axis="x", pad=1)
    ax.set_xlabel(f"residue of the {axis_name} protein", fontsize=5, labelpad=1)


def key(fig, used_sites, sources, top=None, paralog=None):
    """Key at the bottom of the page (or below the genes, its top at top mm): items one after the other, a new line
    when the page width is reached; one group for the gene models, one for the proteins. Returns its height and its
    bottom (mm)."""
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
              (lambda ax, cx, y: (secis_box(ax, cx + 0.3, cx + 2.1, y), secis_box(ax, cx + 2.8, cx + 4.6, y, False),
                                  ax.text(cx + 5.1, y, "37", fontsize=4.5, color=SECIS, va="center")),
               "best SECIS hit and bit score (Oikopleura model; open box: below 20 bits)", "black"),
              (lambda ax, cx, y: ax.text(cx + 7, y, "◀ n", fontsize=4.5, color=PUBLIC, va="center", ha="right"),
               "number of exons beyond the drawn region", "black"),
              (lambda ax, cx, y: ax.text(cx + 6, y, "a b", fontsize=5, fontweight="bold", color=PUBLIC, va="center", ha="right"),
               "public models of the same row, as named in the paragraph", "black")]
    if paralog:
        models.insert(3, (lambda ax, cx, y: exon(ax, cx, y, PARALOG),
                          f"{paralog} (public model; aligned on its codon aligned to the selenocysteine)", "black"))
    sites = ["GT-AG"] + [k for k in SITE_COLOURS if k in used_sites] + \
            (["other"] if any(v not in SITE_COLOURS and v not in ("GT-AG", "?", "") for v in used_sites) else [])
    introns = [(None, "intron splice sites (drawn length unchanged):", "black")]
    introns += [(chev("0.4" if k == "GT-AG" else SITE_COLOURS.get(k, "black"), 0.4 if k == "GT-AG" else 0.9),
                {"GT-AG": "GT-AG (colour of the model)", "other": "other site"}.get(k, k),
                "0.4" if k == "GT-AG" else SITE_COLOURS.get(k, "black")) for k in sites]
    proteins = [(lambda ax, cx, y: ax.add_patch(Rectangle((cx + 1, y - BAR_H / 2), 6, BAR_H, facecolor=DOM, edgecolor="black", lw=0.5)),
                 "Pfam domain of the first reference", "black"),
                (circle("U", SEC, SEC, "white"), "selenocysteine aligned to one of the first reference", "black"),
                (circle("C", "white", "black", "black"), "other residue aligned to it", "black"),
                (lambda ax, cx, y: ax.plot([cx + 1, cx + 7], [y, y], color="black", lw=0.6), "residues not aligned to the first reference", "black"),
                (lambda ax, cx, y: ax.text(cx + 7, y, "92", fontsize=4.5, color="black", va="center", ha="right"),
                 f"% similar residues to the {rf.NAMES[SIM_TO]} protein (BLOSUM62 > 0; global alignment, aligned part)", "black"),
                (None, sources, "black")]
    lines, line, cx = [], [], 0
    split = 5 if paralog else 4
    for group in (models[:split], models[split:], introns, proteins):
        for f, text, col in group:
            w = (G if f else 0) + width_mm(text, 5.5) + 3.5
            if line and cx + w > W: lines.append(line); line, cx = [], 0
            line.append((cx, f, text, col)); cx += w
        lines.append(line); line, cx = [], 0
    h = 3.4 * len(lines)
    bottom = BOTTOM if top is None else top - h
    ax = fig.add_axes([LEFT / PAGE_W, bottom / PAGE_H, W / PAGE_W, h / PAGE_H])
    ax.set_xlim(0, W); ax.set_ylim(-h, 0); ax.axis("off")
    for k, line in enumerate(lines):
        y = -1.7 - 3.4 * k
        for cx, f, text, col in line:
            if f: f(ax, cx, y)
            ax.text(cx + (G if f else 0), y, text, fontsize=5.5, va="center", color=col)
    return h, bottom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", required=True, help="genes in reading order, three per row")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--inputs", default="figures/inputs")
    ap.add_argument("--cm", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--secis-bits", type=float, default=20.0, help="SECIS hits below this score are drawn as open boxes")
    ap.add_argument("--margin", type=int, default=0, help="nt drawn beyond the gold-standard transcripts (0: 15%% of the longest, 150 to 400)")
    ap.add_argument("--long", type=int, default=1500, help="stretches without exon longer than this (nt) are shortened")
    ap.add_argument("--short", type=int, default=200, help="drawn length (nt) of a shortened stretch")
    ap.add_argument("--name", default="reannotation_page")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--sources", default="human: UniProt; C. intestinalis: stage 1 predictions on its genome (accessions in the paragraphs)",
                    help="last line of the key: where the reference proteins come from")
    ap.add_argument("--labels", nargs="*", default=[], help="displayed names, as gene=label (for example TXNL1=TXNL1.b)")
    ap.add_argument("--fit", action="store_true", help="key right below the genes, and the page cut below the key")
    a = ap.parse_args()
    global SECIS_BITS
    SECIS_BITS = a.secis_bits
    genes = [load(g, a) for g in a.genes]
    labels = dict(x.split("=", 1) for x in a.labels)
    col_w = (PAGE_W - LEFT - RIGHT - 2 * COL_GAP) / 3
    fig = plt.figure(figsize=(PAGE_W * MM, PAGE_H * MM))
    fx = lambda mm: mm / PAGE_W
    fy = lambda mm: mm / PAGE_H                    # from the bottom of the page
    top = PAGE_H - TOP
    used, rulers = set(), []
    for d in genes:
        for g in GENOMES:
            used |= set(d["gold_sites"].get(g, "").split(";"))
            for m in d["public"].get(g, {}).values(): used |= set(m.get("sites", "").split(";"))
    for r in range(0, len(genes), 3):
        row = genes[r:r + 3]
        tracks, hA = layout_models(row)
        paras, indent, hN = layout_notes(row, col_w)
        refs_y, gold_y, par_y, hB = layout_proteins(row)
        axA_bottom = top - 3.5 - hA
        axB_bottom = axA_bottom - 6.0 - hB
        axN_bottom = axB_bottom - 7.0 - hN          # paragraph below the proteins (axis and its label above it)
        row_top, row_bottom = top + 1.5, axN_bottom - 1.0
        for k, d in enumerate(row):
            left = LEFT + k * (col_w + COL_GAP)
            fig.text(fx(left), fy(top), labels.get(d["gene"], d["gene"]), fontsize=8, fontweight="bold", va="top")
            axA = fig.add_axes([fx(left), fy(axA_bottom), fx(col_w), fy(hA)])
            axA.set_ylim(-hA, 0)
            draw_models(axA, d, tracks, col_w)
            axN = fig.add_axes([fx(left), fy(axN_bottom), fx(col_w), fy(hN)])
            axN.set_xlim(0, col_w); axN.set_ylim(-hN, 0)
            draw_notes(axN, paras[d["gene"]], indent)
            axB = fig.add_axes([fx(left), fy(axB_bottom), fx(col_w), fy(hB)])
            axB.set_ylim(-hB, 0)
            draw_proteins(axB, d, refs_y, gold_y, par_y, col_w, fig)
            if k == 0:
                lab = blended_transform_factory(fig.transFigure, axA.transData)
                for g in GENOMES:
                    t = tracks[g]
                    axA.text(fx(LEFT - 1.5), t["pub"], DB[g], transform=lab, fontsize=5, ha="right", va="center", color="0.25")
                    axA.text(fx(LEFT - 1.5), t["gold"], "gold standard", transform=lab, fontsize=5.5, ha="right", va="center", color=GOLD)
                    axA.text(fx(2), t["pub"], NAMES[g], transform=lab, fontsize=5.5, fontweight="bold", ha="left", va="center")
                    pm = next((x["par"][g] for x in row if g in x["par"]), None)
                    if "par" in t and pm:
                        axA.text(fx(LEFT - 1.5), t["par"], f"{pm.get('label', 'paralog')}, {pm.get('database', '')}", transform=lab,
                                 fontsize=5.5, ha="right", va="center", color=PARALOG)
                axA.text(fx(LEFT - 1.5), -hA - 0.5, "nt from TGA", transform=lab, fontsize=5, ha="right", va="top", color="0.25")
                labB = blended_transform_factory(fig.transFigure, axB.transData)
                for sp, yy in refs_y.items():
                    axB.text(fx(LEFT - 1.5), yy, sp, transform=labB, fontsize=5.5, ha="right", va="center",
                             style="normal" if sp == "human" else "italic")
                for g, yy in zip(GENOMES, gold_y):
                    axB.text(fx(LEFT - 1.5), yy, NAMES[g], transform=labB, fontsize=5.5, ha="right", va="center", color=GOLD)
                axB.text(fx(4), (gold_y[0] + gold_y[-1]) / 2, "gold standard", transform=labB, fontsize=5.5, rotation=90,
                         ha="center", va="center", color=GOLD)
                if par_y:
                    plabel = next((next(iter(x["par"].values())).get("label", "paralog") for x in row if x["par"]), "paralog")
                    for g, yy in zip(GENOMES, par_y):
                        axB.text(fx(LEFT - 1.5), yy, NAMES[g], transform=labB, fontsize=5.5, ha="right", va="center", color=PARALOG)
                    axB.text(fx(9), (par_y[0] + par_y[-1]) / 2, plabel, transform=labB, fontsize=5.5, rotation=90,
                             ha="center", va="center", color=PARALOG)
        rulers.append((row_top, row_bottom))
        top = axN_bottom - 6
    plabel = next((next(iter(d["par"].values())).get("label", "paralog") for d in genes if d["par"]), None)
    key_h, key_bottom = key(fig, used, a.sources, top=top - 1 if a.fit else None, paralog=plabel)
    # thin rulers: between rows of genes and above the key, and between the columns of each row
    from matplotlib.lines import Line2D
    rule = dict(color="0.8", lw=0.4, transform=fig.transFigure)
    W = PAGE_W - LEFT - RIGHT
    for k, (rt, rb) in enumerate(rulers):
        y = (rb + (rulers[k + 1][0] if k + 1 < len(rulers) else key_bottom + key_h)) / 2
        fig.add_artist(Line2D([fx(2), fx(PAGE_W - RIGHT)], [fy(y), fy(y)], **rule))
        top_ = (rt + (rulers[k - 1][1] if k else PAGE_H)) / 2 if k else rt
        for c in (1, 2):
            x = LEFT + c * col_w + (c - 0.5) * COL_GAP
            fig.add_artist(Line2D([fx(x), fx(x)], [fy(rt), fy(rb)], **rule))
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, a.name)
    from matplotlib.transforms import Bbox
    crop = Bbox.from_extents(0, (key_bottom - BOTTOM) * MM, PAGE_W * MM, PAGE_H * MM) if a.fit else None
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{out}.{ext}", dpi=200 if ext == "png" else None, bbox_inches=crop)
    print("written", out, f"(lowest paragraph ends {top + 6:.0f} mm from the bottom of the page; key {key_bottom + key_h:.0f} mm)"
          + (f"; page cut to {PAGE_H - key_bottom + BOTTOM:.0f} mm" if a.fit else ""))


if __name__ == "__main__":
    main()
