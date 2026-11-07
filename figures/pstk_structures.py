#!/usr/bin/env python3
"""Intron-exon structures of PSTK in the species of the manuscript, drawn on a multiple alignment of the proteins.

Rows: RefSeq models of vertebrates, amphioxus and six other tunicates (intron positions of make pstk-introns,
results/tables/pstk_intron_positions.tsv; proteins from results/trees/pstk.faa), the gold-standard models of the three
O. dioica assemblies, O. albicans and O. vanhoeffeni (exon chains of curation/gold_standard_machinery, splice sites read
in the regions of make pstk-regions), the same genes in the assemblies of other individuals (KUM-M3-7f, AOM-5-5f,
Bar2_p4: their proteins from make pstk-regions aligned to their regions with exonerate, intron positions read from the
alignment), and, for O. longicauda, M. erythrocephalus and B. stygius, which have no model, the tblastn HSPs of the
O. vanhoeffeni protein in their regions (legacy BLAST, E at most --evalue): each HSP is drawn on the columns of the
O. vanhoeffeni residues it covers, and the genome distance between consecutive HSPs is written where they meet, as a
possible intron whose position and sites the HSPs do not fix. The proteins with a sequence are aligned with MAFFT
L-INS-i; each intron is drawn at the alignment column of the residue that follows it (phase 0) or that contains it
(phase 1 or 2, shifted by a third of a column per phase). Marker shape: phase; colour: splice sites (GT-AG, another
O. dioica intron class, outside the classes, or not read for RefSeq models). Vertical lines: the introns of human
PSTK (grey) and of OdB3 PSTK (purple). Introns whose placement is ambiguous (bases repeated at both ends) are drawn
at an intron-class placement when one exists, else at an AG acceptor. Between consecutive HSPs, a bar spans the
columns where the two HSPs end and start (an intron, if any, lies there). T. democratica (XM_082946145.1) has no protein here and is drawn as a label.
Outputs: --out (PDF, PNG, SVG) and --table (TSV of the drawn introns and HSPs).
"""
import argparse, csv, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "curation"))
from intron_placements import region, gold_exons  # noqa: E402
from intron_map import introns_of, exonerate_introns, read_fasta, CLASSES  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROWS = [  # group, label, kind, key
    ("Vertebrates", "Homo sapiens", "refseq", "NP_001350460.1"),
    ("Vertebrates", "Mus musculus", "refseq", "NP_001034623.1"),
    ("Vertebrates", "Danio rerio (1)", "refseq", "XP_009296681.1"),
    ("Vertebrates", "Danio rerio (2)", "refseq", "XP_009305284.1"),
    ("Cephalochordate", "Branchiostoma lanceolatum", "refseq", "XP_066278120.1"),
    ("Other tunicates", "Ciona intestinalis", "refseq", "XP_002120985.2"),
    ("Other tunicates", "Ciona robusta", "refseq", "XP_002120985.1"),
    ("Other tunicates", "Phallusia mammillata", "refseq", "XP_079949190.1"),
    ("Other tunicates", "Styela clava", "refseq", "XP_039273489.2"),
    ("Other tunicates", "Botryllus schlosseri", "refseq", "XP_080181757.1"),
    ("Other tunicates", "Thalia democratica", "nosequence", "XM_082946145.1"),
    ("O. dioica", "OKI2018 (Okinawa)", "gold", "OKI"),
    ("O. dioica", "KUM-M3-7f (Okinawa)", "lift", "KUM"),
    ("O. dioica", "OSKA2016 (North Pacific)", "gold", "OSKA"),
    ("O. dioica", "AOM-5-5f (North Pacific)", "lift", "AOM"),
    ("O. dioica", "OdB3 (North Atlantic)", "gold", "OdB3"),
    ("O. dioica", "Bar2_p4 (North Atlantic)", "lift", "Bar2"),
    ("Other Oikopleura", "O. albicans", "gold", "Oalb"),
    ("Other Oikopleura", "O. vanhoeffeni", "gold", "Ovan"),
    ("Other Oikopleura", "O. longicauda", "fragments", "Olon"),
    ("Giant larvaceans", "M. erythrocephalus", "fragments", "Mery"),
    ("Giant larvaceans", "B. stygius", "fragments", "Bsty"),
]
SITE_COLOURS = {"GT-AG": "black", "class": "#7b3294", "other": "#d7191c", "unknown": "#999999"}
PHASE_MARKERS = {0: "v", 1: "D", 2: "s"}


def shown(res, ph, site, pl):
    """The placement drawn for an intron among its equivalent placements (pl: shifts and sites): an intron class
    nearest to the given placement, else an AG acceptor nearest to it, else the given placement."""
    for ok in (lambda st: st in CLASSES, lambda st: st.endswith("-AG")):
        cand = [(k, st) for k, st in pl if ok(st)]
        if cand:
            k, st = min(cand, key=lambda t: (abs(t[0]), t[0]))
            c = 3 * res + ph + k
            return c // 3, c % 3, st
    return res, ph, site


def tblastn_hsps(protein, fasta, evalue):
    with tempfile.TemporaryDirectory() as d:
        recs = open(fasta).read()
        open(f"{d}/t.fa", "w").write(recs); open(f"{d}/q.fa", "w").write(f">q\n{protein}\n")
        subprocess.run(["formatdb", "-i", "t.fa", "-p", "F", "-o", "F"], cwd=d, check=True, stdout=subprocess.DEVNULL)
        out = subprocess.run(["blastall", "-p", "tblastn", "-d", f"{d}/t.fa", "-i", f"{d}/q.fa", "-e", "10",
                              "-F", "F", "-m", "8"], capture_output=True, text=True, check=True).stdout
    hsps = []
    for l in out.splitlines():
        f = l.split("\t")
        if float(f[10]) > evalue: continue      # legacy BLAST drops some HSPs below -e when it is set low
        hsps.append(dict(contig=f[1].split("|")[1].split(":")[0], record=f[1], identity=float(f[2]), qs=int(f[6]),
                         qe=int(f[7]), ss=int(f[8]), se=int(f[9]), evalue=float(f[10])))
    return hsps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", default="results/tables/pstk_intron_positions.tsv")
    ap.add_argument("--proteins", default="results/trees/pstk.faa")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--regions", default="results/gene_regions/PSTK")
    ap.add_argument("--fragment-query", default="Ovan")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--out", default="figures/pstk_structures")
    ap.add_argument("--table", default="results/tables/pstk_structures.tsv")
    a = ap.parse_args()

    refseq = {}
    for n, s in read_fasta(a.proteins).items():
        m = re.search(r"([NX]P_\d+\.\d+)", n)
        if m: refseq[m.group(1)] = s.replace("*", "").replace("-", "")
    gold = {}
    for g in ["OKI", "OSKA", "OdB3", "Oalb", "Ovan"]:
        gold[g] = next(s for n, s in read_fasta(f"{a.gold_dir}/{g}.proteins.faa").items() if n.startswith("PSTK|")).rstrip("*")
    lifts = {n.split("|")[1]: s.rstrip("*") for n, s in read_fasta(f"{a.regions}/proteins.faa").items()}
    positions = [r for r in csv.DictReader(open(a.positions), delimiter="\t")]

    seqs, introns, table = {}, {}, []      # row key -> protein; row key -> [(residue, phase, site, length)]
    for group, label, kind, key in ROWS:
        if kind == "refseq":
            seqs[key] = refseq[key]
            introns[key] = [(int(float(r["residue"])), int(r["phase"]), None, int(r["intron_length"]))
                            for r in positions if r["id"] == key and r["residue"]]
        elif kind == "gold":
            seqs[key] = gold[key]
            r = region(f"{a.regions}/{key}.fa")
            introns[key] = [shown(res, ph, site, pl) + (length,) for res, ph, length, site, pl in
                            introns_of(gold_exons(f"{a.gold_dir}/{key}.gff3", "PSTK", r), r["gs"])]
        elif kind == "lift":
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
    ncol = len(next(iter(aln.values())))
    cols = {k: [j for j, ch in enumerate(v) if ch != "-"] for k, v in aln.items()}   # residue index -> column

    def xpos(key, res, phase):          # intron after residue res (0-based count of residues before it), phase
        c = cols[key]
        if res >= len(c): return c[-1] + 0.5
        return c[res] - 0.5 + phase / 3 if phase else c[res] - 0.5

    def identity(k1, k2):
        pairs = [(x, y) for x, y in zip(aln[k1], aln[k2]) if x != "-" and y != "-"]
        return 100 * sum(x == y for x, y in pairs) / len(pairs)

    rows = []      # one row per contig for the species without a model
    for group, label, kind, key in ROWS:
        if kind != "fragments": rows.append((group, label, kind, key, None)); continue
        by_contig = {}
        for h in tblastn_hsps(gold[a.fragment_query], f"{a.regions}/{key}.fa", a.evalue):
            by_contig.setdefault(h["contig"], []).append(h)
        for contig, hs in sorted(by_contig.items(), key=lambda t: min(h["evalue"] for h in t[1])):
            rows.append((group, f"{label}, {contig}", kind, key, hs))
    fig_h = 0.42 * len(rows) + 2.2
    fig, ax = plt.subplots(figsize=(13, fig_h))
    yof = {}
    y = 0
    group_y = {}
    for group, label, kind, key, hs in rows:
        yof[(label, key)] = y; group_y.setdefault(group, []).append(y); y -= 1
    # reference lines
    for res, ph, site, length in introns["NP_001350460.1"]:
        ax.axvline(xpos("NP_001350460.1", res, ph), color="#bbbbbb", lw=1, ls="--", zorder=0)
    for res, ph, site, length in introns["OdB3"]:
        ax.axvline(xpos("OdB3", res, ph), color="#c2a5cf", lw=1, ls=":", zorder=0)
    for group, label, kind, key, hs in rows:
        yy = yof[(label, key)]
        ax.text(-8, yy, label, ha="right", va="center", fontsize=8.5, style="italic" if kind == "refseq" else "normal")
        if kind in ("refseq", "gold", "lift"):
            c = cols[key]
            ax.plot([c[0], c[-1]], [yy, yy], color="#555555", lw=1, zorder=1)            # span, gaps thin
            runs, start = [], c[0]
            for p, q in zip(c, c[1:] + [None]):
                if q is None or q != p + 1: runs.append((start, p)); start = q
            for s0, s1 in runs:
                ax.add_patch(plt.Rectangle((s0 - 0.5, yy - 0.18), s1 - s0 + 1, 0.36, color="#9ecae1", lw=0, zorder=2))
            last_x = None
            for res, ph, site, length in introns[key]:
                x = xpos(key, res, ph)
                cat = "unknown" if site is None else "GT-AG" if site == "GT-AG" else "class" if site in CLASSES else "other"
                ax.plot(x, yy + 0.32, marker=PHASE_MARKERS[ph], color=SITE_COLOURS[cat], ms=6, zorder=4, ls="none")
                if site and site != "GT-AG":
                    near = last_x is not None and x - last_x < 25
                    ax.text(x + (1.5 if near else -1.5), yy + 0.36, site, fontsize=6.5, color=SITE_COLOURS[cat],
                            va="bottom", ha="left" if near else "right")
                last_x = x
                table.append([label, key, "intron", res, ph, site or "", length, f"{x:.2f}", ""])
            ident = f"{len(seqs[key])} aa; {identity(key, 'NP_001350460.1'):.0f}% / {identity(key, 'OdB3'):.0f}%"
            ax.text(ncol + 6, yy, ident, va="center", fontsize=7.5)
        elif kind == "nosequence":
            ax.text(cols["XP_002120985.2"][0], yy, "RefSeq XM_082946145.1: 312 codons, one exon (protein not available here, "
                    "not aligned)", va="center", fontsize=7.5, color="#555555")
        else:   # tblastn HSPs of the fragment query
            q = a.fragment_query
            contig = hs[0]["contig"]
            note = []
            for dy in (0.0,):
                strand = 1 if sum((h["se"] > h["ss"]) - (h["se"] < h["ss"]) for h in hs) >= 0 else -1
                hs = sorted([h for h in hs if ((h["se"] > h["ss"]) == (strand == 1))], key=lambda h: h["qs"])
                for i_, h in enumerate(hs):
                    x0, x1 = cols[q][h["qs"] - 1], cols[q][h["qe"] - 1]
                    off = 0.05 if i_ % 2 else -0.05
                    ax.add_patch(plt.Rectangle((x0 - 0.5, yy + dy + off - 0.08), x1 - x0 + 1, 0.16, facecolor="#fdae6b",
                                               edgecolor="#a63603", lw=0.6, zorder=2))
                    ax.text((x0 + x1) / 2, yy + dy - 0.16, f"{h['identity']:.0f}%, E {h['evalue']:.0e}", fontsize=5.5,
                            ha="center", va="top", color="#a63603")
                    table.append([label, key, "hsp", f"{h['qs']}-{h['qe']}", "", contig, f"{h['ss']}-{h['se']}",
                                  f"{x0}-{x1}", f"{h['identity']:.1f} {h['evalue']:.1e}"])
                for h1, h2 in zip(hs, hs[1:]):
                    gap = (h2["ss"] - h1["se"]) * strand - 1
                    xa, xb = sorted((cols[q][min(h1["qe"], len(cols[q])) - 1], cols[q][h2["qs"] - 1]))
                    x = (xa + xb) / 2
                    lab = f"{gap} nt between HSPs" if gap >= 0 else f"HSPs overlap {-gap} nt"
                    ax.plot([xa, xb] if xb - xa > 1 else [x - 0.8, x + 0.8], [yy + dy + 0.2] * 2, color="#a63603", lw=2,
                            zorder=4, solid_capstyle="butt")
                    ax.text(x, yy + dy + 0.24, lab, fontsize=5.5, ha="center", va="bottom", color="#a63603")
                    table.append([label, key, "between HSPs", f"{h1['qe']}|{h2['qs']}", "", contig, gap, f"{x:.1f}", ""])
            ax.text(ncol + 6, yy, f"tblastn HSPs of the O. vanhoeffeni\nprotein (best E {min(h['evalue'] for h in hs):.0e})",
                    va="center", fontsize=6.5)
    for group, ys in group_y.items():
        ax.text(-168, (max(ys) + min(ys)) / 2, group, ha="right", va="center", fontsize=9, weight="bold")
        ax.plot([-164, -164], [min(ys) - 0.3, max(ys) + 0.3], color="black", lw=1)
    # axis: human residues
    hcols = cols["NP_001350460.1"]
    ticks = [hcols[i - 1] for i in range(1, len(hcols) + 1) if i % 50 == 0]
    ax.set_xticks(ticks); ax.set_xticklabels([str(i) for i in range(1, len(hcols) + 1) if i % 50 == 0], fontsize=8)
    ax.set_xlabel("Alignment columns, labelled with the residues of human PSTK (NP_001350460.1)", fontsize=9)
    ax.set_xlim(-172, ncol + 75); ax.set_ylim(y + 0.3, 1.2); ax.set_yticks([])
    for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
    ax.text(ncol + 6, 0.9, "length; identity to human / OdB3", fontsize=7.5, weight="bold")
    handles = [Line2D([], [], marker=PHASE_MARKERS[p], ls="none", color="black", label=f"phase {p}") for p in (0, 1, 2)]
    handles += [Line2D([], [], marker="o", ls="none", color=c, label=l) for l, c in
                [("GT-AG", SITE_COLOURS["GT-AG"]), ("GC-AG, GA-AG, GG-AG", SITE_COLOURS["class"]),
                 ("outside the intron classes", SITE_COLOURS["other"]), ("RefSeq (sites not read)", SITE_COLOURS["unknown"])]]
    handles += [Line2D([], [], color="#bbbbbb", ls="--", label="human introns"), Line2D([], [], color="#c2a5cf", ls=":", label="OdB3 introns"),
                plt.Rectangle((0, 0), 1, 1, color="#9ecae1", label="aligned residues"),
                plt.Rectangle((0, 0), 1, 1, color="#fdae6b", label="tblastn HSP (identity to O. vanhoeffeni)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=5, fontsize=7.5, frameon=False)
    ax.set_title("PSTK: intron positions on the protein alignment", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{a.out}.{ext}", dpi=200 if ext == "png" else None)
    with open(a.table, "w") as o:
        o.write("row\tkey\tfeature\tresidue_or_query\tphase\tsites_or_contig\tlength_or_genome\tcolumn\tidentity_evalue\n")
        for r in table: o.write("\t".join(map(str, r)) + "\n")
    print(f"{a.out}.pdf/png/svg; {a.table}")


if __name__ == "__main__":
    main()
