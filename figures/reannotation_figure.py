#!/usr/bin/env python3
"""Figure of the reannotation of one selenoprotein gene in the three O. dioica assemblies.

A. Gene models in gene orientation, aligned on the selenocysteine TGA in every assembly: the public models (OKI2018
   RefSeq, OSKA2016 ANISEED gene models, OdB3 GenBank) above the gold-standard model of step 2
   (curation/gold_standard/<genome>.gff3). Introns are coloured by splice site (GT-AG in the colour of the model),
   without changing their drawn length; the sites are those of the inputs table. SECIS elements are the best hits of the Oikopleura SECIS model of step 3
   (results/models/oikopleura_secis.cm, Infernal cmsearch) in the 2,000 nt after the stop codon
   (<genome>.after_stop.fa), shortened at the next gold-standard coding sequence downstream on the same strand, drawn
   when cmsearch is available. Exons more than --margin nt outside the
   gold-standard transcript of their assembly are counted, not drawn.
B. Gold-standard proteins and reference proteins drawn on the coordinates of the first reference protein through
   global alignments (BLOSUM62, free end gaps), with its Pfam domains (InterPro, short names) and the residues
   aligned to its selenocysteines.
Inputs: figures/inputs/reannotation_<gene>.tsv and .fa, written by figures/reannotation_inputs.py (public models,
reference proteins, curator notes, each value with its source). Outputs: figures/reannotation_<gene>.{pdf,svg,png}.
"""
import argparse, collections, csv, math, os, re, shutil, subprocess, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.transforms import blended_transform_factory

plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 7})
GOLD, PUBLIC, OTHER, SEC, SECIS, DOM = "#2b6cb0", "#8a8a8a", "#dd8a3c", "#c53030", "#2f855a", "#d6e4f0"
NAMES = {"OKI": "OKI2018", "OSKA": "OSKA2016", "OdB3": "OdB3"}
GENOMES = ["OKI", "OSKA", "OdB3"]
SITE_COLOURS = {"GA-AG": "#6b46c1", "GC-AG": "#d53f8c", "GG-AG": "#b7791f"}   # other sites (not an O. dioica class): black
H_CDS, H_UTR, CHEVRON = 0.2, 0.08, 0.07   # heights of coding and UTR exons, and of the intron chevrons
ROW, GAP = 0.5, 0.95                       # public models to gold standard; gold standard to the next assembly


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].strip(); out[name] = ""
        elif name is not None: out[name] += l.strip()
    return out


def blocks(s): return [tuple(map(int, x.split("-"))) for x in s.split(";")] if s else []


def site_map(exons, sites):
    """{(last base of an exon, first base of the next one): splice site} from the sites of the introns of a chain in
    genome order (exons that touch are not separated by an intron)."""
    ex = sorted(exons)
    introns = [(e1[1], e2[0]) for e1, e2 in zip(ex, ex[1:]) if e2[0] > e1[1] + 1]
    s = sites.split(";") if sites else []
    return dict(zip(introns, s)) if len(s) == len(introns) else {}


def gold_model(path, gene):
    """Strand, exons, CDS rows (1-based inclusive) and selenocysteine codon of one gold-standard model."""
    exons, cds, sec, strand = [], [], None, "+"
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "CDS"): continue
        if re.search(r"Parent=rna-([^;]+)", f[8]).group(1).split("_", 1)[1] != gene: continue
        strand = f[6]
        if f[2] == "exon": exons.append((int(f[3]), int(f[4])))
        else:
            cds.append((int(f[3]), int(f[4])))
            sec = tuple(map(int, re.search(r"pos:(?:complement\()?(\d+)\.\.(\d+)", f[8]).groups()))
    return strand, sorted(exons), sorted(cds), sec


def next_cds(path, gene):
    """Distance from the end of the stop codon of a gold-standard model to the first coding base of the next
    gold-standard model downstream on the same contig and strand (None when there is none)."""
    rows = collections.defaultdict(list)
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        rows[re.search(r"Parent=rna-([^;]+)", f[8]).group(1).split("_", 1)[1]].append((f[0], f[6], int(f[3]), int(f[4])))
    if gene not in rows: return None
    contig, strand = rows[gene][0][:2]
    if strand == "+":
        end = max(r[3] for r in rows[gene])
        d = [min(r[2] for r in v) - end - 1 for k, v in rows.items() if k != gene and v[0][:2] == (contig, strand)]
    else:
        end = min(r[2] for r in rows[gene])
        d = [end - max(r[3] for r in v) - 1 for k, v in rows.items() if k != gene and v[0][:2] == (contig, strand)]
    d = [x for x in d if x >= 0]
    return min(d) if d else None


def secis_hits(cm, gold, gene):
    """Best SECIS hit of a covariance model after the stop codon of each genome: genome -> (from, to, bits). The
    window (<genome>.after_stop.fa, 2,000 nt) is shortened at the next gold-standard coding sequence downstream on the
    same strand, as for the training windows of step 3, so that a gene is not given the SECIS of the next one."""
    if not shutil.which("cmsearch") or not os.path.exists(cm): return {}
    limit = {g: next_cds(f"{gold}/{g}.gff3", gene) for g in GENOMES if os.path.exists(f"{gold}/{g}.gff3")}
    with tempfile.TemporaryDirectory() as d:
        q = os.path.join(d, "q.fa")
        with open(q, "w") as o:
            for g in GENOMES:
                for k, v in read_fasta(f"{gold}/{g}.after_stop.fa").items():
                    if k.startswith(gene + "|"): o.write(f">{g}\n{v}\n")
        tbl = os.path.join(d, "t.tbl")
        subprocess.run(["cmsearch", "--toponly", "-T", "0", "--tblout", tbl, cm, q], check=True, stdout=subprocess.DEVNULL)
        best = {}
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split()
            if limit.get(f[0]) is not None and int(f[8]) > limit[f[0]]: continue
            if float(f[14]) > best.get(f[0], (0, 0, -1))[2]: best[f[0]] = (int(f[7]), int(f[8]), float(f[14]))
    return best


def site_colour(site, color):
    if site in ("GT-AG", "?", ""): return color
    return SITE_COLOURS.get(site, "black")


def draw_model(ax, y, exons, cds, x, color, hatch=False, sites=None):
    """Exons (UTR thin, coding thick) joined by intron chevrons; sites: {(last base of the exon before, first base of
    the exon after): splice site}, drawn as the colour of the chevron."""
    ex = sorted(exons)
    introns = [(e1, e2) for e1, e2 in zip(ex, ex[1:]) if e2[0] > e1[1] + 1]
    for e1, e2 in introns:
        r1, r2 = sorted([(min(x(e1[0]), x(e1[1])), max(x(e1[0]), x(e1[1])) + 1), (min(x(e2[0]), x(e2[1])), max(x(e2[0]), x(e2[1])) + 1)])
        a1, b0 = r1[1], r2[0]
        site = (sites or {}).get((e1[1], e2[0]), "")
        c = site_colour(site, color)
        ax.plot([a1, (a1 + b0) / 2, b0], [y, y + CHEVRON, y], color=c, lw=0.5 if c == color else 1.1,
                solid_joinstyle="miter", zorder=3 if c != color else 2)
    xs = [(min(x(a), x(b)), max(x(a), x(b)) + 1) for a, b in ex]
    for a0, a1 in xs:
        ax.add_patch(Rectangle((a0, y - H_UTR / 2), a1 - a0, H_UTR, facecolor=color, edgecolor="none"))
    for a, b in cds:
        a0, a1 = min(x(a), x(b)), max(x(a), x(b)) + 1
        ax.add_patch(Rectangle((a0, y - H_CDS / 2), a1 - a0, H_CDS, facecolor="white" if hatch else color,
                               edgecolor=color, lw=0.5, hatch="////" if hatch else None))


STEM = 22   # width of the SECIS glyph in nt (set from the panel width in main)


def stem_loop(ax, x0, y, color):
    """SECIS glyph standing on a model drawn at y."""
    ax.plot([x0, x0], [y + 0.04, y + 0.18], color=color, lw=0.8)
    ax.plot([x0 + STEM, x0 + STEM], [y + 0.04, y + 0.18], color=color, lw=0.8)
    ax.add_patch(Ellipse((x0 + STEM / 2, y + 0.23), STEM * 1.6, 0.1, facecolor="none", edgecolor=color, lw=0.8))


def protein_map(q, ref, aligner):
    """Coordinates (1-based) on ref of the aligned residues of q: {query index: ref position}."""
    bq, br = aligner.align(q.replace("U", "C"), ref.replace("U", "C"))[0].aligned
    return {i: j + 1 for (qa, qb), (ra, rb) in zip(bq, br) for i, j in zip(range(qa, qb), range(ra, rb))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--inputs", default="figures/inputs")
    ap.add_argument("--cm", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--margin", type=int, default=0, help="nt drawn beyond the gold-standard transcripts (0: 30%% of the longest, 250 to 600)")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    gene = a.gene
    rows = list(csv.DictReader(open(f"{a.inputs}/reannotation_{gene}.tsv"), delimiter="\t"))
    seqs = read_fasta(f"{a.inputs}/reannotation_{gene}.fa")
    public = collections.defaultdict(dict)
    gold_notes, gold_sites, refs = {}, {}, collections.OrderedDict()
    for r in rows:
        if r["kind"] == "public": public[r["assembly"]].setdefault(r["name"], {})[r["field"]] = r["value"]
        elif r["kind"] == "gold": (gold_sites if r["field"] == "sites" else gold_notes)[r["assembly"]] = r["value"]
        elif r["kind"] == "reference":
            d = refs.setdefault(r["name"], {"domains": []})
            if r["field"] == "domain": d["domains"].append(r["value"])
            else: d[r["field"]] = r["value"]
    gold = {g: gold_model(f"{a.gold}/{g}.gff3", gene) for g in GENOMES}
    gprot = {}
    for g in GENOMES:
        name = None
        for l in open(f"{a.gold}/{g}.proteins.faa"):
            if l.startswith(">"): name = l[1:].split("|")[0]
            elif name == gene: gprot[g] = l.strip()
    secis = secis_hits(a.cm, a.gold, gene)
    # x: nucleotides from the first base of the selenocysteine TGA, in gene orientation, in every assembly
    axes_x = {g: ((lambda p, s0=sec[0]: p - s0) if st == "+" else (lambda p, s1=sec[1]: s1 - p))
              for g, (st, ex, cds, sec) in gold.items()}
    span = {g: (min(min(axes_x[g](p), axes_x[g](q)) for p, q in gold[g][1]),
                max(max(axes_x[g](p), axes_x[g](q)) + 1 for p, q in gold[g][1])) for g in GENOMES}
    longest = max(b - a_ for a_, b in span.values())
    if not a.margin: a.margin = int(min(600, max(250, 0.3 * longest)))
    X0, X1 = min(v[0] for v in span.values()) - a.margin, max(v[1] for v in span.values()) + a.margin
    W = X1 - X0
    box = dict(facecolor="white", edgecolor="none", pad=0.4)
    global STEM
    STEM = 22 * W / 3100

    fig = plt.figure(figsize=(7.2, 5.0))
    axA = fig.add_axes([0.25, 0.40, 0.73, 0.57]); axB = fig.add_axes([0.25, 0.07, 0.73, 0.22])
    labA = blended_transform_factory(axA.transAxes, axA.transData)
    labB = blended_transform_factory(axB.transAxes, axB.transData)
    y, used_sites, notes = 0.0, set(), []
    for g in GENOMES:
        st, ex, cds, sec = gold[g]
        x = axes_x[g]
        lo, hi = span[g][0] - a.margin, span[g][1] + a.margin
        axA.text(-0.33, y - ROW / 2, NAMES[g], transform=labA, fontsize=7.5, fontweight="bold", va="center", ha="left")
        models = sorted(public[g].items(), key=lambda kv: min(min(x(p), x(q)) for p, q in blocks(kv[1]["exons"])))
        for i, (name, m) in enumerate(models):
            e, c = blocks(m["exons"]), blocks(m.get("cds", ""))
            inside = lambda b: lo <= min(x(b[0]), x(b[1])) and max(x(b[0]), x(b[1])) <= hi
            e_in, c_in = [b for b in e if inside(b)], [b for b in c if inside(b)]
            hatch = m.get("frame") == "other"
            sm = site_map(e, m.get("sites", ""))
            used_sites |= {v for k, v in sm.items() if any(b == k[0] for _, b in e_in) and any(a_ == k[1] for a_, _ in e_in)}
            draw_model(axA, y, e_in, c_in, x, OTHER if hatch else PUBLIC, hatch=hatch, sites=sm)
            if not e_in: continue
            xm = min(min(x(p), x(q)) for p, q in e_in); xM = max(max(x(p), x(q)) + 1 for p, q in e_in)
            n_left = sum(max(x(p), x(q)) < lo for p, q in e); n_right = sum(min(x(p), x(q)) > hi for p, q in e)
            if n_left: axA.text(xm, y - 0.15, f"◀ {n_left} more exons", fontsize=5.5, color=PUBLIC, va="top")
            if n_right: axA.text(xM, y - 0.15, f"{n_right} more exons ▶", fontsize=5.5, color=PUBLIC, va="top", ha="right")
            note = m.get("note", "")
            if len(models) == 1:
                if note: notes.append(axA.text(max(xm, span[g][0]), y + 0.15, note, fontsize=5.5, color=PUBLIC, va="bottom", bbox=box, zorder=5))
            else:   # several models on the row: each name and note above or below it, alternately
                txt = f"{name}: {note}" if note else name
                above = i % 2 == 0
                axA.text((xm + xM) / 2, y + (0.15 if above else -0.17), txt, fontsize=5.5, color=OTHER if hatch else PUBLIC,
                         ha="center", va="bottom" if above else "top", bbox=box, zorder=5)
        db = next(iter(public[g].values()))["database"] if public[g] else "public"
        axA.text(-0.02, y + 0.07, db, transform=labA, fontsize=6, ha="right", va="center")
        names = [n.replace("Oidioi.CG.OSKA2016.", "") for n, _ in models]
        axA.text(-0.02, y - 0.11, ", ".join(names) if len(names) < 3 else f"{len(names)} models", transform=labA,
                 fontsize=5, color="0.35", ha="right", va="center")
        y -= ROW
        sm = site_map(ex, gold_sites.get(g, ""))
        used_sites |= set(sm.values())
        draw_model(axA, y, ex, cds, x, GOLD, sites=sm)
        axA.text(-0.02, y, "gold standard", transform=labA, fontsize=6, color=GOLD, ha="right", va="center")
        if g in gold_notes:
            notes.append(axA.text(span[g][0], y - 0.17, gold_notes[g], fontsize=5.5, color=GOLD, va="top", bbox=box, zorder=5))
        if g in secis:
            f0, f1, bits = secis[g]
            stop = max(q for _, q in cds) if st == "+" else min(p for p, _ in cds)
            stem_loop(axA, x(stop) + f0, y, SECIS)
            axA.text(x(stop) + f0 + 1.8 * STEM, y + 0.2, f"SECIS {bits:.0f} bits", fontsize=5.5, color=SECIS, va="center", bbox=box, zorder=5)
        y_last = y
        y -= GAP
    axA.plot([1.5, 1.5], [y_last - 0.15, 0.3], color=SEC, lw=0.7, ls=(0, (2, 1.5)), zorder=1)
    axA.text(1.5, 0.32, "TGA", color=SEC, ha="center", fontsize=5.5, fontweight="bold", va="bottom")
    # key
    ky = y + GAP - 0.55
    unit = W / 2300   # key glyphs scaled to the width of the panel
    at = lambda f: X0 + f * W
    draw_model(axA, ky, [(1, int(120 * unit))], [(int(41 * unit), int(120 * unit))], lambda p: p - 1 + at(0.01), "0.4")
    axA.text(at(0.01) + 140 * unit, ky, "exons: coding (thick), UTR (thin)", fontsize=5.5, va="center")
    draw_model(axA, ky, [(1, int(80 * unit))], [(1, int(80 * unit))], lambda p: p - 1 + at(0.36), OTHER, hatch=True)
    axA.text(at(0.36) + 100 * unit, ky, "coding sequence in another frame", fontsize=5.5, va="center")
    stem_loop(axA, at(0.68), ky - 0.12, SECIS)
    axA.text(at(0.68) + 45 * unit, ky, "SECIS element (Oikopleura model)", fontsize=5.5, va="center")
    shown = [k for k in SITE_COLOURS if k in used_sites] + \
            (["other"] if any(v not in SITE_COLOURS and v not in ("GT-AG", "?", "") for v in used_sites) else [])
    bottom = ky
    if shown:
        bottom = ky - 0.27
        axA.set_xlim(X0, X1)   # the key items are placed one after the other, from the width of their labels
        renderer = fig.canvas.get_renderer()
        right = lambda t: axA.transData.inverted().transform(t.get_window_extent(renderer).get_points())[1][0]
        x0 = right(axA.text(at(0.01), bottom, "intron splice sites:", fontsize=5.5, va="center")) + 40 * unit
        for name in ["GT-AG"] + shown:
            c = "0.4" if name == "GT-AG" else SITE_COLOURS.get(name, "black")
            axA.plot([x0, x0 + 30 * unit, x0 + 60 * unit], [bottom - 0.04, bottom + 0.04, bottom - 0.04], color=c,
                     lw=0.5 if name == "GT-AG" else 1.1, solid_joinstyle="miter")
            label = {"GT-AG": "GT-AG (colour of the model)", "other": "other site"}.get(name, name)
            x0 = right(axA.text(x0 + 75 * unit, bottom, label, fontsize=5.5, va="center", color=c)) + 60 * unit
    axA.set_xlim(X0, X1); axA.set_ylim(bottom - 0.2, 0.55)
    renderer = fig.canvas.get_renderer()
    for t in notes:   # a note that would run beyond the right end of the panel ends there
        x1 = axA.transData.inverted().transform(t.get_window_extent(renderer).get_points())[1][0]
        if x1 > X1: t.set_x(t.get_position()[0] - (x1 - X1))
    axA.set_yticks([]); [axA.spines[s_].set_visible(False) for s_ in ("left", "right", "top")]
    step = 500 if W < 4500 else 1000
    axA.set_xticks(range(math.ceil(X0 / step) * step, int(X1) + 1, step))
    axA.set_xlabel("nucleotides from the selenocysteine TGA (gene orientation)", fontsize=6)
    axA.text(-0.33, 0.42, "A", transform=labA, fontsize=10, fontweight="bold")

    # ---- proteins on the coordinates of the first reference
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    al.end_insertion_score = al.end_deletion_score = 0
    axis_label = next(k for k, v in refs.items() if v.get("axis") == "yes")
    axis = seqs[axis_label]
    secpos = [int(p) for p in refs[axis_label].get("selenocysteine", "").split(",") if p]
    doms = []
    for d in refs[axis_label]["domains"]:
        acc, span, name = d.split(":", 2)
        d0, d1 = map(int, span.split("-"))
        doms.append((d0, d1, f"{name} ({acc})"))
    prots = [(k, v.get("accession", ""), seqs[k], "black") for k, v in refs.items()] + \
            [(f"{NAMES[g]} gold standard", "", gprot[g], GOLD) for g in GENOMES]
    for i, (lbl, acc, p, col) in enumerate(prots):
        yy = -i
        m = protein_map(p, axis, al)
        h0, h1 = min(m.values()), max(m.values())
        axB.add_patch(Rectangle((h0, yy - 0.15), h1 - h0, 0.3, facecolor="white", edgecolor=col, lw=0.7))
        for d0, d1, dn in doms:
            a0, a1 = max(d0, h0), min(d1, h1)
            if a1 > a0: axB.add_patch(Rectangle((a0, yy - 0.15), a1 - a0, 0.3, facecolor=DOM, edgecolor=col, lw=0.7))
        for sp in secpos:
            site = next((p[i_] for i_, j in m.items() if j == sp), None)
            if site:
                axB.plot([sp], [yy], "o", ms=6, color=SEC if site == "U" else "white", mec=SEC if site == "U" else "black", mew=0.7)
                axB.text(sp, yy - 0.01, site, ha="center", va="center", fontsize=5.5, color="white" if site == "U" else "black",
                         fontweight="bold")
        axB.text(-0.02, yy + (0.13 if acc else 0), lbl, transform=labB, ha="right", va="center", fontsize=6, color=col)
        note = refs.get(lbl, {}).get("note", "")
        if acc: axB.text(-0.02, yy - 0.2, f"{acc}; {note}" if note else acc, transform=labB, ha="right", va="center", fontsize=5, color="0.35")
        axB.text(h1 + len(axis) * 0.012, yy, f"{len(p)} aa", va="center", fontsize=5.5, color="0.3")
    for d0, d1, dn in doms:
        axB.text((d0 + d1) / 2, 0.3, dn, ha="center", fontsize=5.5, color="0.25", va="bottom")
    for sp in secpos:
        axB.text(sp, 0.62, f"U{sp}", ha="center", fontsize=5.5, color=SEC, fontweight="bold", va="bottom")
    axB.set_xlim(0, len(axis) * 1.07); axB.set_ylim(-len(prots) + 0.5, 0.95)
    axB.set_yticks([]); [axB.spines[s_].set_visible(False) for s_ in ("left", "right", "top")]
    axB.set_xlabel(f"position in {axis_label} (through global alignments)", fontsize=6)
    axB.text(-0.33, 0.75, "B", transform=labB, fontsize=10, fontweight="bold")
    out = os.path.join(a.out, f"reannotation_{gene}")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None)
    print("written", out)


if __name__ == "__main__":
    main()
