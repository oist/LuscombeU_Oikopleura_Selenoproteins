#!/usr/bin/env python3
"""How fast PSTK evolves: PSTK against its paralog KTI12 in the same chordates, and against other selenocysteine genes.

1. PSTK and KTI12 (--out-pairs). The PSTK tree holds both families (KTI12 is its outgroup), so that the distances
   between two species in one family and in the other come from the same alignment and the same model. For each pair
   of chordates with both proteins (PAIRS), the table gives the identity of the two PSTK and of the two KTI12 (global
   alignment of the whole proteins, BLOSUM62, gap opening 10, extension 0.5, end gaps free; identical residues over
   aligned columns) and their distance in the trees (substitutions per site, sum of the branch lengths between the two
   leaves): the whole alignment (results/trees/pstk.treefile, selected model; results/trees/lba/pstk.pmsf.treefile,
   PMSF) and its two parts (results/trees/domains/pstk.{kinase,C_terminal}[.pmsf].treefile; identity there over the
   gap-free columns of the trimmed alignment), with the ratio of the PSTK distance to the KTI12 distance.
   Zebrafish: the first of its two PSTK genes (XP_009296681.1).
2. Other genes (--out-genes). For each gene with human and Ciona reference proteins in figures/inputs/reannotation_*.fa
   (one row per human reference) and an OKI2018 gold-standard protein (curation/gold_standard*/OKI.proteins.faa): the
   identity between the human and Ciona proteins, the human and O. dioica proteins, and the Ciona and O. dioica
   proteins, aligned as above (selenocysteine read as cysteine).
"""
import argparse, glob, os, statistics
from Bio import Phylo, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

PSTK = {"human": "Human_PSTK_NP_001350460.1", "mouse": "Mmus_RefSeq_PSTK_NP_001034623.1",
        "zebrafish": "Drer_RefSeq_PSTK_XP_009296681.1", "amphioxus": "Blan_RefSeq_PSTK_XP_066278120.1",
        "Ciona": "Cint_RefSeq_PSTK_XP_002120985.2", "Phallusia": "Pmam_RefSeq_PSTK_XP_079949190.1",
        "Styela": "Scla_RefSeq_PSTK_XP_039273489.2", "Botryllus": "Bsch_RefSeq_PSTK_XP_080181757.1",
        "O. dioica": "OKI_gold_PSTK"}
KTI12 = {"human": "Human_KTI12_NP_612426.1", "mouse": "Mmus_RefSeq_KTI12_NP_083847.1",
         "zebrafish": "Drer_RefSeq_KTI12_NP_001119890.1", "amphioxus": "Blan_RefSeq_KTI12_XP_066270687.1",
         "Ciona": "Cint_RefSeq_KTI12_XP_026695236.2", "Phallusia": "Pmam_RefSeq_KTI12_XP_079953735.1",
         "Styela": "Scla_RefSeq_KTI12_XP_039256908.2", "Botryllus": "Bsch_RefSeq_KTI12_XP_080176166.1",
         "O. dioica": "KTI12_KTI12_Oikopleura_dioica_XP_081641927.1"}
PAIRS = [("human", "mouse"), ("human", "zebrafish"), ("human", "amphioxus"), ("human", "Ciona"), ("human", "O. dioica"),
         ("Ciona", "Phallusia"), ("Ciona", "Styela"), ("Ciona", "Botryllus"), ("Ciona", "O. dioica")]

ALIGNER = PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.mode, ALIGNER.open_gap_score, ALIGNER.extend_gap_score, ALIGNER.end_gap_score = "global", -10, -0.5, 0


def identity(a, b):
    """Percent identical residues over the aligned columns of a global alignment."""
    clean = lambda s: s.replace("*", "").replace("-", "").replace("U", "C")
    aln = ALIGNER.align(clean(a), clean(b))[0]
    cols = [(x, y) for x, y in zip(aln[0], aln[1]) if x != "-" and y != "-"]
    return 100 * sum(x == y for x, y in cols) / len(cols)


def column_identity(a, b):
    """Percent identity over the gap-free columns of two rows of an alignment."""
    cols = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    return 100 * sum(x == y for x, y in cols) / len(cols) if cols else float("nan")


def fasta(path):
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def pairs_table(a):
    seqs = fasta(a.faa)
    trees = [("whole", "selected model", a.tree, None), ("whole", "PMSF", a.pmsf_tree, None)]
    for part in ("kinase", "C_terminal"):
        trees += [(part, "selected model", f"{a.domains}/pstk.{part}.treefile", f"{a.domains}/pstk.{part}.fa"),
                  (part, "PMSF", f"{a.domains}/pstk.{part}.pmsf.treefile", f"{a.domains}/pstk.{part}.fa")]
    rows = []
    for part, model, tpath, apath in trees:
        if not os.path.exists(tpath): continue
        t = Phylo.read(tpath, "newick")
        leaf = {c.name: c for c in t.get_terminals()}
        aln = fasta(apath) if apath else None
        for s1, s2 in PAIRS:
            if apath:
                ip, ik = (column_identity(aln[F[s1]], aln[F[s2]]) for F in (PSTK, KTI12))
            else:
                ip, ik = (identity(seqs[F[s1]], seqs[F[s2]]) for F in (PSTK, KTI12))
            dp = t.distance(leaf[PSTK[s1]], leaf[PSTK[s2]]); dk = t.distance(leaf[KTI12[s1]], leaf[KTI12[s2]])
            rows.append([part, model, s1, s2, f"{ip:.1f}", f"{ik:.1f}", f"{dp:.3f}", f"{dk:.3f}", f"{dp / dk:.2f}"])
    with open(a.out_pairs, "w") as h:
        h.write("region\tmodel\tspecies1\tspecies2\tPSTK_identity\tKTI12_identity\tPSTK_distance\tKTI12_distance\tratio\n")
        for r in rows: h.write("\t".join(r) + "\n")
    whole = [float(r[8]) for r in rows if r[0] == "whole" and r[1] == "selected model" and r[2:4] != ["human", "mouse"]]
    print(f"{a.out_pairs}: PSTK/KTI12 distance ratio, whole alignment, selected model, all pairs but human-mouse: "
          f"{min(whole):.2f} to {max(whole):.2f}")


def genes_table(a):
    gold = {}
    for f in a.gold:
        for name, s in fasta(f).items(): gold.setdefault(name.split("|")[0], s)
    rows, seen = [], set()
    for f in sorted(glob.glob(f"{a.references}/reannotation_*.fa")):
        gene = os.path.basename(f)[len("reannotation_"):-len(".fa")]
        ref = fasta(f)
        ciona = next((v for k, v in ref.items() if k.startswith("C.")), None)
        if "human" not in ref or ciona is None or gene not in gold or ref["human"] in seen: continue
        seen.add(ref["human"])
        rows.append([gene, identity(ref["human"], ciona), identity(ref["human"], gold[gene]), identity(ciona, gold[gene])])
    rows.sort(key=lambda r: r[1])
    with open(a.out_genes, "w") as h:
        h.write("gene\thuman_Ciona\thuman_Odioica\tCiona_Odioica\n")
        for r in rows: h.write(r[0] + "".join(f"\t{x:.1f}" for x in r[1:]) + "\n")
    hc = [r[1] for r in rows]
    print(f"{a.out_genes}: {len(rows)} genes; human-Ciona identity median {statistics.median(hc):.1f}, "
          f"lowest {rows[0][0]} {rows[0][1]:.1f}, highest {rows[-1][0]} {rows[-1][1]:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", default="results/trees/pstk.faa")
    ap.add_argument("--tree", default="results/trees/pstk.treefile")
    ap.add_argument("--pmsf-tree", default="results/trees/lba/pstk.pmsf.treefile")
    ap.add_argument("--domains", default="results/trees/domains")
    ap.add_argument("--references", default="figures/inputs")
    ap.add_argument("--gold", nargs="+", default=["curation/gold_standard/OKI.proteins.faa",
                                                   "curation/gold_standard_machinery/OKI.proteins.faa"])
    ap.add_argument("--out-pairs", default="results/tables/pstk_kti12_distances.tsv")
    ap.add_argument("--out-genes", default="results/tables/sec_gene_identities.tsv")
    a = ap.parse_args()
    pairs_table(a)
    genes_table(a)


if __name__ == "__main__":
    main()
