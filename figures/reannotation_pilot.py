#!/usr/bin/env python3
"""Pilot figure for the reannotation of one selenoprotein gene (SEPHS2): public annotations and gold-standard models
of the three O. dioica assemblies (A), and their proteins against human SEPHS2 and a Ciona SEPHS (B).

A. Gene models in gene orientation, aligned on the 5' end of the gold-standard transcript of each assembly: public
   models (OKI2018 RefSeq, OSKA2016 ANISEED gene models, OdB3 GenBank) above the gold-standard models of step 2
   (curation/gold_standard/<genome>.gff3; exon chains and atypical introns from <genome>.candidate_models.tsv).
   The selenocysteine TGA is drawn across each assembly; SECIS elements are the best hits of the current
   appendicularian model (results/models/appendicularia_secis.cm, Infernal cmsearch) in the 2,000 nt after the stop
   codon (<genome>.after_stop.fa), drawn when cmsearch is available.
B. Gold-standard proteins, with human SEPHS2 and a Ciona SEPHS, drawn on the coordinates of human SEPHS2 (UniProt
   Q99611) through global alignments (BLOSUM62, free end gaps): AIR synthase related protein domains (Pfam PF00586 and PF02769, InterPro coordinates on Q99611) and the
   residue at the selenocysteine position of human SEPHS2 (U60).
External inputs (public, with sources): figures/inputs/reannotation_pilot_sephs2.tsv and .fa. The Ciona protein comes
from the SEPHS2 profile of results/models/machinery_profiles/.
Outputs: figures/reannotation_pilot_SEPHS2.{pdf,svg,png}.
"""
import argparse, csv, os, re, shutil, subprocess, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "font.size": 7})
GOLD, PUBLIC, OTHER, SEC, SECIS, DOM = "#2b6cb0", "#8a8a8a", "#dd8a3c", "#c53030", "#2f855a", "#d6e4f0"
NAMES = {"OKI": "OKI2018", "OSKA": "OSKA2016", "OdB3": "OdB3"}


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].split()[0]; out[name] = ""
        elif name: out[name] += l.strip()
    return out


def parse_location(s):
    """GenBank location (join/complement) -> strand, [(start, end)] 1-based inclusive."""
    strand = "-" if s.startswith("complement") else "+"
    return strand, [tuple(map(int, x.split(".."))) for x in re.findall(r"\d+\.\.\d+", s)]


def blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(";")] if s else []


def gold_model(path, gene):
    """Strand, exons, CDS rows (1-based inclusive) and selenocysteine codon of one gold-standard model."""
    exons, cds, sec, strand = [], [], None, "+"
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "CDS"): continue
        parent = re.search(r"Parent=rna-([^;]+)", f[8]).group(1)
        if parent.split("_", 1)[1] != gene: continue
        strand = f[6]
        if f[2] == "exon": exons.append((int(f[3]), int(f[4])))
        else:
            cds.append((int(f[3]), int(f[4])))
            sec = tuple(map(int, re.search(r"pos:(?:complement\()?(\d+)\.\.(\d+)", f[8]).groups()))
    return strand, sorted(exons), sorted(cds), sec


def cds_from_residue(cds, strand, residue):
    """CDS segments of a model from codon `residue` (1-based) to its end."""
    pos = [p for a, b in sorted(cds) for p in range(a, b + 1)]
    if strand == "-": pos = pos[::-1]
    keep = sorted(pos[(residue - 1) * 3:])
    segs, start = [], keep[0]
    for x, y in zip(keep, keep[1:] + [None]):
        if y != x + 1: segs.append((start, x)); start = y
    return segs


def secis_hits(cm, fasta, genomes, gene):
    """Best SECIS hit of the current model in the after-stop sequence of each genome: genome -> (from, to, bits)."""
    if not shutil.which("cmsearch") or not os.path.exists(cm): return {}
    with tempfile.TemporaryDirectory() as d:
        q = os.path.join(d, "q.fa")
        with open(q, "w") as o:
            for g in genomes:
                seqs = read_fasta(fasta.format(g))
                k = next((k for k in seqs if k.startswith(gene + "|")), None)
                if k: o.write(f">{g}\n{seqs[k]}\n")
        tbl = os.path.join(d, "t.tbl")
        subprocess.run(["cmsearch", "--toponly", "-T", "0", "--tblout", tbl, cm, q], check=True, stdout=subprocess.DEVNULL)
        best = {}
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split()
            if float(f[14]) > best.get(f[0], (0, 0, -1))[2]: best[f[0]] = (int(f[7]), int(f[8]), float(f[14]))
    return best


H_CDS, H_UTR, CHEVRON = 0.2, 0.08, 0.07   # heights of coding and UTR exons, and of the intron chevrons
ROW, GAP = 0.5, 0.95                       # public model to gold standard; gold standard to the next assembly


def draw_model(ax, y, exons, cds, x, color, hatch=False, h_cds=H_CDS, h_utr=H_UTR):
    xs = sorted((min(x(a), x(b)), max(x(a), x(b)) + 1) for a, b in exons)
    for (a0, a1), (b0, b1) in zip(xs, xs[1:]):
        mid = (a1 + b0) / 2
        ax.plot([a1, mid, b0], [y, y + CHEVRON, y], color=color, lw=0.5, solid_joinstyle="miter")
    for a0, a1 in xs:
        ax.add_patch(Rectangle((a0, y - h_utr / 2), a1 - a0, h_utr, facecolor=color, edgecolor="none"))
    for a, b in cds:
        a0, a1 = min(x(a), x(b)), max(x(a), x(b)) + 1
        ax.add_patch(Rectangle((a0, y - h_cds / 2), a1 - a0, h_cds, facecolor="white" if hatch else color,
                               edgecolor=color, lw=0.5, hatch="////" if hatch else None))


def stem_loop(ax, x0, y, color):
    """SECIS glyph standing on a model drawn at y."""
    ax.plot([x0, x0], [y + 0.04, y + 0.18], color=color, lw=0.8)
    ax.plot([x0 + 22, x0 + 22], [y + 0.04, y + 0.18], color=color, lw=0.8)
    ax.add_patch(matplotlib.patches.Ellipse((x0 + 11, y + 0.23), 36, 0.1, facecolor="none", edgecolor=color, lw=0.8))


def protein_map(q, ref, aligner):
    """Human coordinates (1-based) of the aligned residues of q: {query index: human position}."""
    bq, br = aligner.align(q.replace("U", "C"), ref.replace("U", "C"))[0].aligned
    return {i: j + 1 for (qa, qb), (ra, rb) in zip(bq, br) for i, j in zip(range(qa, qb), range(ra, rb))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="SEPHS2")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--data", default="figures/inputs/reannotation_pilot_sephs2")
    ap.add_argument("--profile", default="results/models/machinery_profiles/SEPHS2.fa")
    ap.add_argument("--cm", default="results/models/appendicularia_secis.cm")
    ap.add_argument("--out", default="figures/reannotation_pilot_SEPHS2")
    a = ap.parse_args()
    gene, genomes = a.gene, ["OKI", "OSKA", "OdB3"]
    facts = list(csv.DictReader(open(a.data + ".tsv"), delimiter="\t"))
    seqs = read_fasta(a.data + ".fa")
    cand = {g: [r for r in csv.DictReader(open(f"{a.gold}/{g}.candidate_models.tsv"), delimiter="\t") if r["gene"] == gene]
            for g in genomes}
    secis = secis_hits(a.cm, f"{a.gold}/{{}}.after_stop.fa", genomes, gene)

    # ---- public models (exons, CDS, hatched, label) per genome
    gold = {g: gold_model(f"{a.gold}/{g}.gff3", gene) for g in genomes}
    public = {}
    st, ex, cds, sec = gold["OKI"]
    rs = next(r for r in cand["OKI"] if r["kind"] == "annotation" and r["source"].startswith("XP_"))
    res = int(next(f["value"] for f in facts if f["model"] == rs["source"] and f["feature"].startswith("protein start")))
    public["OKI"] = [(f"RefSeq {rs['source']}", blocks(rs["exons"]), cds_from_residue(cds, st, res), False)]
    gst, gex, _, _ = gold["OSKA"]
    lo, hi = min(a_ for a_, _ in gex), max(b for _, b in gex)
    osk = []
    for r in cand["OSKA"]:
        if r["item"].startswith("annotated") and r["splice_sites"].startswith("CDS"):
            b = blocks(r["exons"])
            if any(x <= hi and y >= lo for x, y in b): osk.append((r["source"], b))
    public["OSKA"] = []
    for src, b in osk:
        inside = [(x, y) for x, y in b if y >= lo - 600 and x <= hi + 600]
        public["OSKA"].append((f"ANISEED {src}" + (f" (and {len(b) - len(inside)} exons upstream)" if len(inside) < len(b) else ""),
                               inside, inside, False))
    od = []
    for model in ("CBY08599.1", "CBY08598.1"):
        loc = {f["feature"]: f["value"] for f in facts if f["model"] == model}
        od.append((f"GenBank {model}", parse_location(loc["mRNA"])[1], parse_location(loc["CDS"])[1], model == "CBY08599.1"))
    public["OdB3"] = od

    # ---- figure
    from matplotlib.transforms import blended_transform_factory
    fig = plt.figure(figsize=(7.2, 5.0))
    axA = fig.add_axes([0.25, 0.40, 0.73, 0.57]); axB = fig.add_axes([0.25, 0.07, 0.73, 0.22])
    labA = blended_transform_factory(axA.transAxes, axA.transData)
    labB = blended_transform_factory(axB.transAxes, axB.transData)
    gprot = {}
    for g in genomes:
        name = None
        for l in open(f"{a.gold}/{g}.proteins.faa"):
            if l.startswith(">"): name = l[1:].split("|")[0]
            elif name == gene: gprot[g] = l.strip()
    u = gprot["OKI"].index("U") + 1
    notes = {"OKI": f"coding sequence starts at codon {res}, {res - u} codons after the TGA (codon {u})",
             "OSKA": "TGA exon skipped; fused with the upstream gene"}
    y = 0.0
    for g in genomes:
        st, ex, cds, sec = gold[g]
        t0, t1 = min(x_ for x_, _ in ex), max(y_ for _, y_ in ex)
        x = (lambda p, t0=t0: p - t0) if st == "+" else (lambda p, t1=t1: t1 - p)
        xs = x(sec[0]) if st == "+" else x(sec[1])
        top = y
        axA.text(-0.33, y - ROW / 2, NAMES[g], transform=labA, fontsize=7.5, fontweight="bold", va="center", ha="left")
        # public models (one row; OdB3 has two models side by side)
        names = []
        for lbl, e, c, hatch in public[g]:
            draw_model(axA, y, e, c, x, OTHER if hatch else PUBLIC, hatch=hatch)
            names.append(lbl.split()[1])
            if g == "OdB3":
                xm = min(min(x(p_), x(q_)) for p_, q_ in e); xM = max(max(x(p_), x(q_)) for p_, q_ in e)
                txt = f"{lbl.split()[1]}: other frame" if hatch else f"{lbl.split()[1]}: starts at an internal ATG"
                axA.text((xm + xM) / 2 if not hatch else xs - 20, y + 0.15, txt, fontsize=5.5,
                         color=OTHER if hatch else PUBLIC, ha="center" if not hatch else "right", va="bottom")
            elif g == "OSKA" and "upstream" in lbl:
                xm = min(min(x(p_), x(q_)) for p_, q_ in e)
                n = re.search(r"and (\d+) exons", lbl).group(1)
                axA.annotate("", xy=(xm - 80, y), xytext=(xm - 10, y),
                             arrowprops=dict(arrowstyle="-|>", color=PUBLIC, lw=0.5, mutation_scale=5))
                axA.text(xm, y + 0.15, f"◀ {n} more exons", fontsize=5.5, color=PUBLIC, ha="left", va="bottom")
        src = {"OKI": "RefSeq", "OSKA": "ANISEED", "OdB3": "GenBank"}[g]
        axA.text(-0.02, y + 0.07, src, transform=labA, fontsize=6, ha="right", va="center")
        axA.text(-0.02, y - 0.11, ", ".join(n_.replace("Oidioi.CG.OSKA2016.", "") for n_ in names), transform=labA,
                 fontsize=5, color="0.35", ha="right", va="center")
        if g in notes:
            axA.text(xs + 60, y + 0.15, notes[g], fontsize=5.5, color=PUBLIC, va="bottom")
        y -= ROW
        # gold standard
        draw_model(axA, y, ex, cds, x, GOLD)
        axA.text(-0.02, y, "gold standard", transform=labA, fontsize=6, color=GOLD, ha="right", va="center")
        chosen = next((r for r in cand[g] if r["chosen"] == "yes"), None)
        if chosen:
            exs = blocks(chosen["exons"])
            for (a1, b1), (a2, b2), site in zip(exs, exs[1:], chosen["splice_sites"].split(";")):
                if site != "GT-AG":
                    axA.text((x(b1) + x(a2)) / 2, y + 0.09, site, ha="center", fontsize=5, color=GOLD, va="bottom")
        if g in secis:
            f0, f1, bits = secis[g]
            stop = max(b for _, b in cds) if st == "+" else min(a_ for a_, _ in cds)
            xs_ = x(stop) + f0
            stem_loop(axA, xs_, y, SECIS)
            axA.text(xs_ + 40, y + 0.2, f"SECIS {bits:.0f} bits", fontsize=5.5, color=SECIS, va="center")
        axA.plot([xs + 1.5, xs + 1.5], [y - 0.15, top + 0.33], color=SEC, lw=0.7, ls=(0, (2, 1.5)))
        axA.text(xs + 1.5, top + 0.35, "TGA", color=SEC, ha="center", fontsize=5.5, fontweight="bold", va="bottom")
        y -= GAP
    # key
    ky = y + GAP - 0.35
    draw_model(axA, ky, [(1, 120)], [(41, 120)], lambda p: p - 1 - 650, "0.4")
    axA.text(-510, ky, "exons: coding (thick), UTR (thin)", fontsize=5.5, va="center")
    draw_model(axA, ky, [(1, 80)], [(1, 80)], lambda p: p - 1 + 450, OTHER, hatch=True)
    axA.text(550, ky, "coding sequence in another frame", fontsize=5.5, va="center")
    stem_loop(axA, 1500, ky - 0.12, SECIS)
    axA.text(1545, ky, "SECIS element (current model)", fontsize=5.5, va="center")
    axA.set_xlim(-700, 2400); axA.set_ylim(ky - 0.2, 0.55)
    axA.set_yticks([]); [axA.spines[s_].set_visible(False) for s_ in ("left", "right", "top")]
    axA.set_xticks(range(0, 2001, 500))
    axA.set_xlabel("nucleotides from the 5′ end of the gold-standard transcript (gene orientation)", fontsize=6)
    axA.text(-0.33, 0.42, "A", transform=labA, fontsize=10, fontweight="bold")

    # ---- proteins on human coordinates
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    al.end_insertion_score = al.end_deletion_score = 0
    human = seqs["human_SEPHS2"]
    prof = read_fasta(a.profile)
    ciona = next(v for k, v in prof.items() if "198432475" in k).replace("-", "").upper()
    gp = {}
    for g in genomes:
        name = None
        for l in open(f"{a.gold}/{g}.proteins.faa"):
            if l.startswith(">"): name = l[1:].split("|")[0]
            elif name == gene: gp[g] = l.strip()
    prots = [("human SEPHS2", "UniProt Q99611", human, "black"),
             ("Ciona intestinalis SEPHS", "XP_002129526.1", ciona, "black"),
             ("OKI2018 gold standard", "", gp["OKI"], GOLD),
             ("OSKA2016 gold standard", "", gp["OSKA"], GOLD),
             ("OdB3 gold standard", "", gp["OdB3"], GOLD)]
    doms = [tuple(map(int, f["value"].split("-"))) + (f["feature"].split()[1],) for f in facts if f["feature"].startswith("domain")]
    for i, (lbl, acc, p, col) in enumerate(prots):
        yy = -i
        m = protein_map(p, human, al)
        h0, h1 = min(m.values()), max(m.values())
        axB.add_patch(Rectangle((h0, yy - 0.15), h1 - h0, 0.3, facecolor="white", edgecolor=col, lw=0.7))
        for d0, d1, dn in doms:
            a0, a1 = max(d0, h0), min(d1, h1)
            if a1 > a0: axB.add_patch(Rectangle((a0, yy - 0.15), a1 - a0, 0.3, facecolor=DOM, edgecolor=col, lw=0.7))
        site = next((p[i_] for i_, j in m.items() if j == 60), None)
        if site:
            axB.plot([60], [yy], "o", ms=6, color=SEC if site == "U" else "white", mec=SEC if site == "U" else "black", mew=0.7)
            axB.text(60, yy - 0.01, site, ha="center", va="center", fontsize=5.5, color="white" if site == "U" else "black",
                     fontweight="bold")
        axB.text(-0.02, yy + (0.13 if acc else 0), lbl, transform=labB, ha="right", va="center", fontsize=6, color=col)
        if acc: axB.text(-0.02, yy - 0.2, acc, transform=labB, ha="right", va="center", fontsize=5, color="0.35")
        axB.text(h1 + 6, yy, f"{len(p)} aa", va="center", fontsize=5.5, color="0.3")
    for d0, d1, dn in doms:
        axB.text((d0 + d1) / 2, 0.42, {"PF00586": "AIRS", "PF02769": "AIRS_C"}.get(dn, dn) + f" ({dn})", ha="center",
                 fontsize=5.5, color="0.25")
    axB.text(60, 0.42, "U60", ha="center", fontsize=5.5, color=SEC, fontweight="bold")
    axB.set_xlim(0, 480); axB.set_ylim(-len(prots) + 0.5, 0.75)
    axB.set_yticks([]); [axB.spines[s_].set_visible(False) for s_ in ("left", "right", "top")]
    axB.set_xlabel("position in human SEPHS2 (through global alignments)", fontsize=6)
    axB.text(-0.33, 0.6, "B", transform=labB, fontsize=10, fontweight="bold")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300 if ext == "png" else None)
    print("written", a.out)


if __name__ == "__main__":
    main()
