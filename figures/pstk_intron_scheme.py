#!/usr/bin/env python3
"""PSTK introns of chordates and appendicularians on the residues of human PSTK: a gene-structure scheme.

Each protein is placed on the residues of human PSTK (NP_001350460.1) through a multiple alignment (MAFFT L-INS-i):
a residue aligned with a human residue takes its number, residues in columns where human has a gap are collapsed onto
the boundary, and stretches of 10 or more human residues without a residue of the protein are drawn thin. Each gene is
drawn as its exons (alternating shades) joined by intron marks (a peak above the exon boundary; colour: splice-site
dinucleotides; number: phase). Vertical lines: the introns of human PSTK (grey, numbered 1-5) and of OdB3 PSTK
(purple, a-d). Two introns of different species are counted at the same position when they are at most one alignment
column apart and in the same phase.

Panel A, gene models: RefSeq models (intron positions of make pstk-introns, results/tables/pstk_intron_positions.tsv;
proteins from results/trees/pstk.faa), gold-standard models of the O. dioica assemblies, O. albicans and O. vanhoeffeni
(curation/gold_standard_machinery; splice sites read in the regions of make pstk-regions, introns with repeated end
bases shown at an O. dioica intron-class placement when one exists, as in figures/pstk_structures.py); one row per
group whose other members are compared with it and named on the right: mouse and the two zebrafish genes with human;
C. robusta, P. mammillata, S. clava and B. schlosseri with C. intestinalis (T. democratica: RefSeq transcript
XM_082946145.1, one exon, no protein here); the assembly of the other individual of each O. dioica species (KUM-M3-7f,
AOM-5-5f, Bar2_p4: their proteins of make pstk-regions aligned to their regions with exonerate).

Panel B, species without a gene model (O. longicauda, M. erythrocephalus, B. stygius; one row per contig): tblastn
HSPs of the O. vanhoeffeni protein in the regions of make pstk-regions (legacy BLAST, no filter, E at most --evalue;
searched at E 10 and filtered, as legacy BLAST drops HSPs with a low -e), with the stop codons in their frame (red
ticks). Junctions: an insertion of at least --min-insertion codons inside an HSP (gaps of the query), and each pair of
consecutive HSPs of a contig that are collinear on the gene strand. For each junction, every intron of at least 20 nt
with O. dioica class sites (GT-AG, GC-AG, GA-AG, GG-AG) whose donor lies at most --window codons upstream of the end
of the upstream block (or further, up to its first stop codon) and whose acceptor lies at most --window codons
downstream of the start of the downstream block (or further back, down to its last stop codon), and which keeps the
two frames, is spliced in silico; the spliced translation of the two blocks is aligned with the O. vanhoeffeni protein
(local, BLOSUM62, gaps -11/-1) and the placement with the highest score is drawn (dashed peak), with its length; ties
(equivalent placements or equal scores) are counted in the table. The phase and position come from the alignment of
the spliced protein. No transcript supports these placements.

Outputs: --out (PDF, PNG, SVG) and --table (TSV: panel, row, model, feature, position, phase, sites, length, x, note).
"""
import argparse, bisect, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "curation"))
from intron_placements import region, gold_exons, rc, translate  # noqa: E402
from intron_map import introns_of, exonerate_introns, read_fasta, CLASSES  # noqa: E402
from pstk_structures import shown  # noqa: E402
from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HUMAN = "NP_001350460.1"
ROWS = [  # group, label (drawn model), model key, [(other model, its name)] compared with it
    ("Vertebrates", "Homo sapiens", HUMAN,
     [("NP_001034623.1", "mouse"), ("XP_009296681.1", "zebrafish gene 1"), ("XP_009305284.1", "zebrafish gene 2")]),
    ("Amphioxus", "Branchiostoma lanceolatum", "XP_066278120.1", []),
    ("Ascidians, salp", "Ciona intestinalis", "XP_002120985.2",
     [("XP_002120985.1", "C. robusta"), ("XP_079949190.1", "P. mammillata"), ("XP_039273489.2", "S. clava"),
      ("XP_080181757.1", "B. schlosseri")]),
    ("O. dioica", "OKI2018 (Okinawa)", "OKI", [("KUM", "KUM-M3-7f")]),
    ("O. dioica", "OSKA2016 (North Pacific)", "OSKA", [("AOM", "AOM-5-5f")]),
    ("O. dioica", "OdB3 (North Atlantic)", "OdB3", [("Bar2", "Bar2_p4")]),
    ("Other Oikopleura", "O. albicans", "Oalb", []),
    ("Other Oikopleura", "O. vanhoeffeni", "Ovan", []),
]
FRAGMENTS = [("O. longicauda", "Olon"), ("M. erythrocephalus", "Mery"), ("B. stygius", "Bsty")]
GOLD = ["OKI", "OSKA", "OdB3", "Oalb", "Ovan"]
LIFTS = ["KUM", "AOM", "Bar2"]
SITE_COLOURS = {"GT-AG": "black", "class": "#7b3294", "other": "#d7191c", "unknown": "#808080"}
EXON_COLOURS = ("#3b75af", "#9ecae1")
HSP_FACE, HSP_EDGE = "#fdd0a2", "#a63603"


def site_cat(site):
    return "unknown" if not site else "GT-AG" if site == "GT-AG" else "class" if site in CLASSES else "other"


def blast_pairwise(protein, fasta, evalue):
    """tblastn HSPs with their alignment strings (legacy BLAST, -m 0), E at most evalue."""
    with tempfile.TemporaryDirectory() as d:
        recs = read_fasta(fasta)
        open(f"{d}/t.fa", "w").write("".join(f">{n}\n{s.upper()}\n" for n, s in recs.items()))
        open(f"{d}/q.fa", "w").write(f">q\n{protein}\n")
        subprocess.run(["formatdb", "-i", "t.fa", "-p", "F", "-o", "F"], cwd=d, check=True, stdout=subprocess.DEVNULL)
        out = subprocess.run(["blastall", "-p", "tblastn", "-d", f"{d}/t.fa", "-i", f"{d}/q.fa", "-e", "10", "-F", "F",
                              "-m", "0"], capture_output=True, text=True, check=True).stdout
    hsps = []
    for block in out.split("\n>")[1:]:
        name = block.split()[0]
        for h in block.split(" Score =")[1:]:
            ev = re.search(r"Expect(?:\(\d+\))? = (\S+)", h).group(1).rstrip(",")
            ev = float("1" + ev if ev.startswith("e") else ev)
            if ev > evalue: continue
            ident = int(re.search(r"Identities = \d+/\d+ \((\d+)%\)", h).group(1))
            q, s, qp, sp = [], [], [], []
            for l in h.splitlines():
                m = re.match(r"(Query|Sbjct)\s+(\d+)\s+(\S+)\s+(\d+)\s*$", l)
                if not m: continue
                (q if m.group(1) == "Query" else s).append(m.group(3))
                (qp if m.group(1) == "Query" else sp).append((int(m.group(2)), int(m.group(4))))
            hsps.append(dict(record=name, contig=name.split("|")[1].split(":")[0], evalue=ev, identity=ident,
                             qs=qp[0][0], qe=qp[-1][1], ss=sp[0][0], se=sp[-1][1], qaln="".join(q), saln="".join(s)))
    return hsps


def hsp_residues(h, L):
    """[(query index 0-based, gene-strand codon start)] of the aligned residue pairs, the stops of the HSP (query
    index of the column or of the residue before a gap) and the gaps of the query (query boundary, codons, stops)."""
    minus = h["ss"] > h["se"]
    ptr = (L - h["ss"]) if minus else (h["ss"] - 1)       # gene strand: reverse complement for minus-strand HSPs
    qi = h["qs"] - 1
    pairs, stops, gaps, run = [], [], [], None
    for qc, sc in zip(h["qaln"], h["saln"]):
        if qc == "-":
            if run is None: run = [qi, 0, 0, ptr]
            run[1] += 1; run[2] += sc == "*"
        else:
            if run is not None: gaps.append(tuple(run)); run = None
        if sc == "*": stops.append(qi - 0.5 if qc == "-" else qi)
        if qc != "-" and sc != "-": pairs.append((qi, ptr))
        if sc != "-": ptr += 3
        if qc != "-": qi += 1
    return pairs, stops, gaps


def place_junction(gs, up, down, query, aligner, window):
    """Best class-site intron between an upstream block and a downstream block of aligned residues
    [(query index, codon start)], scored by the local alignment of the spliced translation with the query."""
    u0, u_end = up[0][1], up[-1][1] + 3
    d0, d_end = down[0][1], down[-1][1] + 3
    f1, f2 = u0 % 3, d0 % 3
    s_after = next((p for p in range(u_end - 3 * window, len(gs) - 2, 3) if p >= u0 and gs[p:p + 3] in
                    ("TAA", "TAG", "TGA")), len(gs))
    s_before = max([p for p in range(d0 % 3, d0 + 3 * window, 3) if gs[p:p + 3] in ("TAA", "TAG", "TGA")] + [-3])
    donors = range(max(u0 + 3, u_end - 3 * window), s_after + 1)
    acceptors = range(max(s_before + 3, d0 - 3 * window if s_before < 0 else s_before + 3), min(d_end - 3, d0 + 3 * window) + 1)
    qlo, qhi = up[0][0], down[-1][0] + 1
    qseg = query[qlo:qhi]
    cands = []
    for dn in donors:
        for ac in acceptors:
            if ac - dn < 20 or (dn - f1) % 3 != (ac - f2) % 3: continue
            site = gs[dn:dn + 2] + "-" + gs[ac - 2:ac]
            if site not in CLASSES: continue
            spliced = gs[u0:dn] + gs[ac:d_end]
            prot = translate(spliced)
            if "*" in prot: continue
            cands.append((aligner.score(prot, qseg), dn, ac, site, prot))
    if not cands: return None
    best = max(cands, key=lambda c: c[0])
    score, dn, ac, site, prot = best
    ties = sum(1 for c in cands if c[0] >= score - 1e-9)
    al = aligner.align(prot, qseg)[0]
    jres, phase = (dn - u0) // 3, (dn - u0) % 3       # junction: after residue jres of the spliced protein (+ phase)
    qpos = None
    for (ps, pe), (qs, qe) in zip(*al.aligned):
        if ps <= jres <= pe: qpos = qs + (jres - ps); break
    if qpos is None:     # junction residue in a gap of the alignment: nearest aligned residue
        best_d = None
        for (ps, pe), (qs, qe) in zip(*al.aligned):
            for p, qq in ((ps, qs), (pe, qe)):
                if best_d is None or abs(p - jres) < best_d: best_d, qpos = abs(p - jres), qq + (jres - p)
    return dict(query_boundary=qlo + qpos, phase=phase, site=site, length=ac - dn, score=score, ties=ties,
                candidates=len(cands))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", default="results/tables/pstk_intron_positions.tsv")
    ap.add_argument("--proteins", default="results/trees/pstk.faa")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--regions", default="results/gene_regions/PSTK")
    ap.add_argument("--fragment-query", default="Ovan")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--min-insertion", type=int, default=5)
    ap.add_argument("--window", type=int, default=25)
    ap.add_argument("--out", default="figures/pstk_intron_scheme")
    ap.add_argument("--table", default="results/tables/pstk_intron_scheme.tsv")
    a = ap.parse_args()

    import csv
    refseq = {}
    for n, s in read_fasta(a.proteins).items():
        m = re.search(r"([NX]P_\d+\.\d+)", n)
        if m: refseq[m.group(1)] = s.replace("*", "").replace("-", "")
    gold = {g: next(s for n, s in read_fasta(f"{a.gold_dir}/{g}.proteins.faa").items()
                    if n.startswith("PSTK|")).rstrip("*") for g in GOLD}
    lifts = {n.split("|")[1]: s.rstrip("*") for n, s in read_fasta(f"{a.regions}/proteins.faa").items()}
    positions = list(csv.DictReader(open(a.positions), delimiter="\t"))
    keys = [k for _, _, k, others in ROWS for k in [k] + [o for o, _ in others]]
    seqs, introns = {}, {}
    for key in keys:
        if key in refseq:
            seqs[key] = refseq[key]
            introns[key] = [(int(float(r["residue"])), int(r["phase"]), None, int(r["intron_length"]))
                            for r in positions if r["id"] == key and r["residue"]]
        elif key in GOLD:
            seqs[key] = gold[key]
            r = region(f"{a.regions}/{key}.fa")
            introns[key] = [shown(res, ph, site, pl) + (length,) for res, ph, length, site, pl in
                            introns_of(gold_exons(f"{a.gold_dir}/{key}.gff3", "PSTK", r), r["gs"])]
        else:
            seqs[key] = lifts[key]
            r = region(f"{a.regions}/{key}.fa")
            introns[key] = [shown(res, ph, site, pl) + (length,) for res, ph, length, site, pl in
                            exonerate_introns(lifts[key], r, to_class=False)]
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/in.fa", "w") as o:
            for k, s in seqs.items(): o.write(f">{k}\n{s}\n")
        txt = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", f"{d}/in.fa"],
                             capture_output=True, text=True, check=True).stdout
    aln = {b.split("\n", 1)[0].strip(): b.split("\n", 1)[1].replace("\n", "") for b in txt.split(">")[1:]}
    cols = {k: [j for j, ch in enumerate(v) if ch != "-"] for k, v in aln.items()}
    hcols = cols[HUMAN]
    hset = set(hcols)
    nh = len(hcols)

    def hb(c):           # human residues in the columns before column c
        return bisect.bisect_left(hcols, c)

    def boundary_x(key, res, phase=0):   # intron after `res` residues of the protein, in a phase
        c = cols[key]
        if res >= len(c): return hb(c[-1]) + (1 if c[-1] in hset else 0)
        x = hb(c[res])
        return x + phase / 3 if (phase and c[res] in hset) else x

    def colpos(key, res, phase):         # alignment column of an intron (for the comparisons)
        c = cols[key]
        return c[-1] + 0.5 if res >= len(c) else c[res] - 0.5 + phase / 3

    def same_intron(k1, i1, k2, i2):
        return i1[1] == i2[1] and abs(colpos(k1, i1[0], i1[1]) - colpos(k2, i2[0], i2[1])) <= 1.0

    def matches(k1, k2):
        return sum(any(same_intron(k1, i, k2, j) for j in introns[k2]) for i in introns[k1])

    def covered(key):    # human residues aligned with a residue of the protein
        return {hb(c) for c in cols[key] if c in hset}

    table = []
    # fragment rows: HSPs and junctions
    aligner = Align.PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    q = a.fragment_query
    frag_rows = []
    for label, key in FRAGMENTS:
        recs = {n.split("|")[1].split(":")[0]: s.upper() for n, s in read_fasta(f"{a.regions}/{key}.fa").items()}
        by_contig = {}
        for h in blast_pairwise(gold[q], f"{a.regions}/{key}.fa", a.evalue):
            by_contig.setdefault(h["contig"], []).append(h)
        for contig, hs in sorted(by_contig.items(), key=lambda t: min(h["evalue"] for h in t[1])):
            seq = recs[contig]
            L = len(seq)
            minus = sum((h["ss"] > h["se"]) for h in hs) * 2 >= len(hs)
            hs = sorted([h for h in hs if (h["ss"] > h["se"]) == minus], key=lambda h: h["qs"])
            gs = rc(seq) if minus else seq
            blocks, stops, junctions = [], [], []
            for h in hs:
                pairs, st, gaps = hsp_residues(h, L)
                stops += st
                cuts = [g for g in gaps if g[1] >= a.min_insertion]
                cur = []
                for qi, p in pairs:
                    if cuts and qi >= cuts[0][0] and cur:
                        blocks.append(dict(res=cur, hsp=h, cut_after=cuts[0])); cur = []; cuts = cuts[1:]
                    cur.append((qi, p))
                blocks.append(dict(res=cur, hsp=h, cut_after=None))
            for b1, b2 in zip(blocks, blocks[1:]):
                if b2["res"][0][1] <= b1["res"][-1][1]: continue        # not collinear on the gene strand
                inside = b1["hsp"] is b2["hsp"]
                pl = place_junction(gs, b1["res"], b2["res"], gold[q], aligner, a.window)
                gap_nt = b2["res"][0][1] - (b1["res"][-1][1] + 3)
                junctions.append(dict(inside=inside, b1=b1, b2=b2, placement=pl, gap_nt=gap_nt,
                                      ins=b1["cut_after"] if inside else None))
            frag_rows.append((label, key, contig, hs, stops, junctions, L))

    # figure
    n_a, n_b = len(ROWS), len(frag_rows)
    fig_h = 0.55 * (n_a + n_b) + 2.8
    fig, ax = plt.subplots(figsize=(15, fig_h))
    y_a = [-(i * 1.0) for i in range(n_a)]
    y_b = [y_a[-1] - 1.9 - i * 1.0 for i in range(n_b)]
    H = 0.32
    top = 1.15
    bottom = y_b[-1] - 0.7
    human_x = [boundary_x(HUMAN, r, p) for r, p, _, _ in introns[HUMAN]]
    odb3_x = [boundary_x("OdB3", r, p) for r, p, _, _ in introns["OdB3"]]
    for i, x in enumerate(human_x):
        ax.plot([x, x], [bottom, top - 0.25], color="#9a9a9a", lw=0.9, ls="--", zorder=0)
        ax.text(x, top - 0.2, str(i + 1), color="#6a6a6a", ha="center", va="bottom", fontsize=8, weight="bold")
    for i, x in enumerate(odb3_x):
        ax.plot([x, x], [bottom, top - 0.25], color="#b07cc6", lw=0.9, ls=":", zorder=0)
        ax.text(x, top - 0.2, "abcd"[i], color="#7b3294", ha="center", va="bottom", fontsize=8, weight="bold")

    def draw_peak(x, y, colour, text, dashed=False, side="center"):
        ax.plot([x - 2.2, x, x + 2.2], [y + H / 2, y + H / 2 + 0.24, y + H / 2], color=colour, lw=1.3,
                ls="--" if dashed else "-", zorder=5, solid_joinstyle="miter")
        dx = {"center": 0, "left": -1.2, "right": 1.2}[side]
        ha = {"center": "center", "left": "right", "right": "left"}[side]
        ax.text(x + dx, y + H / 2 + 0.27, text, color=colour, ha=ha, va="bottom", fontsize=6.5)

    def sides(xs, near=12):      # labels of close peaks pushed apart
        out = []
        for i, x in enumerate(xs):
            prev_near = i > 0 and x - xs[i - 1] < near
            next_near = i + 1 < len(xs) and xs[i + 1] - x < near
            out.append("right" if prev_near else "left" if next_near else "center")
        return out

    def peak_text(phase, site, length=None):
        return ", ".join([str(phase)] + ([site] if site and site != "GT-AG" else []) +
                         ([f"{length} nt"] if length is not None else []))

    def draw_model(key, y):
        c = cols[key]
        x0, x1 = boundary_x(key, 0), boundary_x(key, len(c))
        cuts = [boundary_x(key, r, p) for r, p, _, _ in introns[key]]
        edges = [x0] + cuts + [x1]
        cov = covered(key)
        thin = []           # runs of >= 10 human residues without a residue of the protein, inside its span
        run = None
        for hres in range(int(x0), int(x1)):
            if hres not in cov:
                run = run or [hres, hres]; run[1] = hres + 1
            else:
                if run and run[1] - run[0] >= 10: thin.append(tuple(run))
                run = None
        if run and run[1] - run[0] >= 10: thin.append(tuple(run))
        for i, (e0, e1) in enumerate(zip(edges, edges[1:])):
            col = EXON_COLOURS[i % 2]
            segs, s = [], e0
            for t0, t1 in thin:
                if t1 <= e0 or t0 >= e1: continue
                if t0 > s: segs.append((s, t0, H))
                segs.append((max(t0, e0), min(t1, e1), 0.08)); s = min(t1, e1)
            if s < e1: segs.append((s, e1, H))
            for s0, s1, hh in segs:
                ax.add_patch(plt.Rectangle((s0, y - hh / 2), s1 - s0, hh, facecolor=col, edgecolor="none", zorder=2))
        for (r, p, site, length), x, sd in zip(introns[key], cuts, sides(cuts)):
            draw_peak(x, y, SITE_COLOURS[site_cat(site)], peak_text(p, site), side=sd)
            table.append(["A", key, key, "intron", r, p, site or "", length, f"{x:.2f}", ""])
        return x0, x1

    for (group, label, key, others), y in zip(ROWS, y_a):
        draw_model(key, y)
        ax.text(-6, y, label, ha="right", va="center", fontsize=9,
                style="italic" if key in refseq or key in ("Oalb", "Ovan") else "normal")
        n = len(introns[key])
        words = ["no intron" if n == 0 else f"{n} intron{'s' if n > 1 else ''}"]
        if key != HUMAN and n:
            at_h = matches(key, HUMAN)
            words.append("all at the human positions" if at_h == n else
                         "none at a human position" if at_h == 0 else f"{at_h} at human positions")
            if key == "Oalb":
                at_d = matches(key, "OdB3")
                words[-1] = words[-1].replace("none at a human position", "none at a human or O. dioica position") \
                    if at_d == 0 else words[-1] + f"; {at_d} at O. dioica positions"
        same = [nm for o, nm in others if len(introns[o]) == n and matches(o, key) == n]
        diff = [nm for o, nm in others if nm not in same]
        if same: words.append("same in " + ", ".join(same))
        if diff: words.append("differs in " + ", ".join(diff))
        if key == "XP_002120985.2": words.append("T. democratica: one exon (RefSeq transcript)")
        ax.text(nh + 10, y, "; ".join(words), va="center", fontsize=7.5, wrap=True)
    groups = {}
    for (group, label, key, others), y in zip(ROWS, y_a): groups.setdefault(group, []).append(y)

    # panel B
    for (label, key, contig, hs, stops, junctions, L), y in zip(frag_rows, y_b):
        ax.text(-6, y, f"{label}, {contig}", ha="right", va="center", fontsize=8.5, style="italic")
        ax.plot([boundary_x(q, 0), boundary_x(q, len(gold[q]))], [y, y], color="#d0d0d0", lw=0.8, zorder=1)
        for i_, h in enumerate(hs):
            x0, x1 = boundary_x(q, h["qs"] - 1), boundary_x(q, h["qe"])
            off = 0.06 if i_ % 2 else -0.06
            ax.add_patch(plt.Rectangle((x0, y + off - 0.11), x1 - x0, 0.22, facecolor=HSP_FACE, edgecolor=HSP_EDGE,
                                       lw=0.7, zorder=2, alpha=0.9))
            ax.text((x0 + x1) / 2, y - 0.2, f"{h['identity']}%, E {h['evalue']:.0e}", fontsize=6, ha="center",
                    va="top", color=HSP_EDGE)
            table.append(["B", contig, key, "hsp", f"{h['qs']}-{h['qe']}", "", "", f"{h['ss']}-{h['se']}",
                          f"{x0:.1f}-{x1:.1f}", f"identity {h['identity']}%, E {h['evalue']:.1e}"])
        for s in stops:
            x = boundary_x(q, int(s)) + (0.5 if s == int(s) else 0)
            ax.plot([x, x], [y - 0.13, y + 0.13], color="#d7191c", lw=1.4, zorder=3)
        notes = []
        for j in junctions:
            pl = j["placement"]
            kind = "insertion in an HSP" if j["inside"] else "between HSPs"
            q1, q2 = j["b1"]["res"][-1][0] + 1, j["b2"]["res"][0][0]
            if j["inside"]:
                gq, ncod, nst, _ = j["ins"]
                ev = f"insertion of {3 * ncod} nt after residue {gq}, {nst} stop codon{'s' if nst != 1 else ''}"
            elif q2 + 1 <= q1:
                ev = f"HSPs overlap on residues {q2 + 1}-{q1}, {j['gap_nt']} nt apart"
            else:
                ev = f"residues {q1 + 1}-{q2} not matched, HSPs {j['gap_nt']} nt apart"
            if pl:
                x = boundary_x(q, pl["query_boundary"], pl["phase"])
                draw_peak(x, y, SITE_COLOURS[site_cat(pl["site"])], peak_text(pl["phase"], pl["site"], pl["length"]),
                          dashed=True)
                table.append(["B", contig, key, "junction", pl["query_boundary"], pl["phase"], pl["site"], pl["length"],
                              f"{x:.2f}", f"{kind}: {ev}; best of {pl['candidates']} class-site placements, "
                              f"score {pl['score']:.0f}, {pl['ties']} with that score"])
            else:
                table.append(["B", contig, key, "junction", f"{q1}|{q2 + 1}", "", "", "", "",
                              f"{kind}: {ev}; no class-site placement"])
            notes.append(ev + ("" if pl else " (no class-site intron)"))
        ax.text(nh + 10, y, "\n".join(notes) if notes else "one HSP, no insertion", va="center", fontsize=7)
    for g, ys in groups.items():
        ax.text(-128, (max(ys) + min(ys)) / 2, g, ha="right", va="center", fontsize=9, weight="bold")
        ax.plot([-125, -125], [min(ys) - 0.3, max(ys) + 0.3], color="black", lw=1)
    yb = (y_b[0] + y_b[-1]) / 2
    ax.text(-128, yb, "No gene model", ha="right", va="center", fontsize=9, weight="bold")
    ax.plot([-125, -125], [y_b[-1] - 0.3, y_b[0] + 0.3], color="black", lw=1)
    ax.text(-200, y_a[0] + 0.85, "A  Gene models", fontsize=10.5, weight="bold", va="bottom")
    ax.text(-200, y_b[0] + 0.85, "B  Species without a gene model: tblastn HSPs of the O. vanhoeffeni protein, and the "
            "best class-site intron at each junction (no transcripts)", fontsize=10.5, weight="bold", va="bottom")
    ax.set_xticks(range(0, nh + 1, 50)); ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Residues of human PSTK (NP_001350460.1); other proteins placed by a MAFFT alignment, "
                  "residues absent from human collapsed", fontsize=9)
    ax.set_xlim(-205, nh + 150); ax.set_ylim(bottom - 0.1, top + 0.35); ax.set_yticks([])
    for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_bounds(0, nh)
    handles = [Patch(facecolor=EXON_COLOURS[0], label="exon"), Patch(facecolor=EXON_COLOURS[1], label="next exon"),
               Line2D([], [], color="black", lw=1.3, label="intron, GT-AG (first number: phase)"),
               Line2D([], [], color=SITE_COLOURS["class"], lw=1.3, label="intron, GC-AG, GA-AG or GG-AG"),
               Line2D([], [], color=SITE_COLOURS["other"], lw=1.3, label="intron, other sites"),
               Line2D([], [], color=SITE_COLOURS["unknown"], lw=1.3, label="intron, RefSeq (sites not read)"),
               Line2D([], [], color="#9a9a9a", ls="--", label="human introns 1-5"),
               Line2D([], [], color="#b07cc6", ls=":", label="OdB3 introns a-d"),
               Patch(facecolor=HSP_FACE, edgecolor=HSP_EDGE, label="tblastn HSP"),
               Line2D([], [], color="#d7191c", lw=1.4, marker="|", ls="none", ms=8, label="stop codon in the HSP frame"),
               Line2D([], [], color="black", lw=1.3, ls="--", label="intron placed from the genome only")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=4, fontsize=7.5, frameon=False)
    ax.set_title("PSTK: intron positions on the residues of human PSTK", fontsize=11.5)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{a.out}.{ext}", dpi=200 if ext == "png" else None)
    with open(a.table, "w") as o:
        o.write("panel\trow\tmodel\tfeature\tposition\tphase\tsites\tlength\tx\tnote\n")
        for r in table: o.write("\t".join(map(str, r)) + "\n")
    print(f"{a.out}.pdf/png/svg; {a.table}")


if __name__ == "__main__":
    main()
