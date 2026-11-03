#!/usr/bin/env python3
"""Intron positions of one gene across species, on the coordinates of a reference protein.

For each model (annotated protein, annotated transcript or gold-standard model), the coding rows are read from the
annotation (samplesheet column "annotation"; RefSeq GFF3 with protein_id and Parent attributes) or from the
gold-standard GFF3, spliced on the genome (samplesheet column "fasta") and translated. Each intron is placed in the
protein (residue after which it lies, and phase: 0 between codons, 1 or 2 inside a codon), and the protein is aligned
to the reference protein (global alignment, BLOSUM62, free end gaps, Biopython), so that each intron gets a position
on the reference. Introns of different species at the same reference position and phase are candidates for a shared
(ancestral) intron; a position shifted by a few residues with another phase is not the same intron.
Models: --reference NAME=ID (the reference, an annotated protein), --models NAME=ID ... (ID: a protein_id, or a
transcript ID such as XM_..., whose coding rows are used), --gold NAME[:GENE] ... (gold-standard model of GENE in
<gold-dir>/<NAME>.gff3; GENE defaults to --gene).
Output (--out, TSV): model, id, source, exons, protein_length, stop (translation ends with a stop codon),
internal_stops, identity (to the reference, over the aligned columns), intron (number in gene order), residue (in the
model protein, x.0 = after that residue, x.33 and x.67 inside the next codon), phase, reference_residue (aligned
position on the reference, with the same convention; empty when the residue is in a gap), intron_length.
"""
import argparse, csv, gzip, os, re

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_contigs(path, names):
    out, name, buf = {}, None, []
    with opener(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in names: out[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            elif name in names: buf.append(l.strip())
    if name in names: out[name] = "".join(buf).upper()
    return out


def annotation_cds(path, ident):
    """Coding rows (contig, start, end, strand) of a protein_id, or of the transcript whose ID or Name is ident."""
    rows, by_parent, tx = [], {}, set()
    attr = lambda f, k: (re.search(rf"(?:^|;){k}=([^;]+)", f[8]) or [None, None])[1]
    with opener(path) as fh:
        for l in fh:
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) != 9: continue
            if f[2] == "CDS":
                if attr(f, "protein_id") == ident: rows.append(f)
                by_parent.setdefault(attr(f, "Parent"), []).append(f)
            elif ident in (attr(f, "ID") or "", attr(f, "Name") or "", attr(f, "transcript_id") or "") or \
                    (attr(f, "ID") or "").endswith("-" + ident):
                tx.add(attr(f, "ID"))
    if not rows:
        for t in tx: rows += by_parent.get(t, [])
    return [(f[0], int(f[3]), int(f[4]), f[6]) for f in rows]


def gold_cds(path, gene):
    rows = []
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "CDS" and re.search(rf"Parent=rna-[^;]*_{re.escape(gene)}(;|$)", f[8]):
            rows.append((f[0], int(f[3]), int(f[4]), f[6]))
    return rows


def model(rows, genome_path):
    """Coding sequence (gene strand), exon lengths in gene order, intron lengths."""
    contig, strand = rows[0][0], rows[0][3]
    seq = read_contigs(genome_path, {contig})[contig]
    ex = sorted((a, b) for _, a, b, _ in rows)
    introns = [b2[0] - b1[1] - 1 for b1, b2 in zip(ex, ex[1:])]
    parts = [seq[a - 1:b] for a, b in ex]
    if strand == "-":
        parts = [p.translate(COMP)[::-1] for p in parts[::-1]]; ex = ex[::-1]; introns = introns[::-1]
    return "".join(parts), [len(p) for p in parts], introns


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="PSTK")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--models", nargs="*", default=[])
    ap.add_argument("--gold", nargs="*", default=[])
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--samplesheet", nargs="+", default=["stage1/samplesheet.csv"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sheet = {}
    for s in a.samplesheet:
        for r in csv.DictReader(open(s)): sheet.setdefault(r["sample"], r)
    specs = [("reference",) + tuple(a.reference.split("=", 1))] + [("annotation",) + tuple(m.split("=", 1)) for m in a.models]
    specs += [("gold", g.split(":")[0], g.split(":")[1] if ":" in g else a.gene) for g in a.gold]
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.mode = "global"; al.end_insertion_score = al.end_deletion_score = 0
    rows_out, ref = [], None
    for kind, name, ident in specs:
        if kind == "gold": rows = gold_cds(f"{a.gold_dir}/{name}.gff3", ident)
        else: rows = annotation_cds(sheet[name]["annotation"], ident)
        if not rows:
            print(f"{name} {ident}: no coding rows"); continue
        cds, exl, intl = model(rows, sheet[name]["fasta"])
        prot = str(Seq(cds[:len(cds) // 3 * 3]).translate())
        stop = prot.endswith("*"); prot = prot.rstrip("*")
        if ref is None: ref = prot
        aln = al.align(prot, ref)[0]
        A, B = aln[0], aln[1]
        m2r, ia, ib = {}, 0, 0             # model residue index (0-based) -> reference residue index or None
        for x, y in zip(A, B):
            if x != "-": m2r[ia] = ib if y != "-" else None
            if x != "-": ia += 1
            if y != "-": ib += 1
        cols = [(x, y) for x, y in zip(A, B) if x != "-" and y != "-"]
        ident_frac = sum(x == y for x, y in cols) / max(len(cols), 1)
        base = dict(model=name, id=ident, source=kind, exons=len(exl), protein_length=len(prot), stop="yes" if stop else "no",
                    internal_stops=prot.count("*"), identity=round(ident_frac, 3))
        if len(exl) == 1: rows_out.append(dict(base, intron="", residue="", phase="", reference_residue="", intron_length=""))
        cum = 0
        for i, (e, il) in enumerate(zip(exl, intl)):
            cum += e
            res, ph = cum // 3, cum % 3     # intron after `res` complete codons; phase = bases of the next codon before it
            r = m2r.get(res if ph else res - 1)
            if r is None: refpos = ""
            else: refpos = round((r + 1 if not ph else r) + ph / 3, 2)
            rows_out.append(dict(base, intron=i + 1, residue=round(res + ph / 3, 2), phase=ph, reference_residue=refpos,
                                 intron_length=il))
    cols = ["model", "id", "source", "exons", "protein_length", "stop", "internal_stops", "identity", "intron", "residue",
            "phase", "reference_residue", "intron_length"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows_out: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows_out)} rows written to {a.out}")


if __name__ == "__main__":
    main()
