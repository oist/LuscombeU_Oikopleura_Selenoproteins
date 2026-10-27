#!/usr/bin/env python3
"""Intron positions of TXNL1.b (gold standard) and TXNL1.a (public model) in the three O. dioica assemblies, on the
coordinates of the TXNL1.a protein.

For each coding intron, the residue of the protein whose codon holds the first coding base after the intron (counted
from 1) and the intron phase (0: between codons; 1 or 2: after the first or the second base of a codon). The TXNL1.b
positions are carried onto TXNL1.a through a local alignment of the two proteins (BLOSUM62, as in txnl1_paralogs.py);
an intron of one gene is shared when the other has an intron of the same phase at the same aligned residue.
Inputs (repository files only): curation/gold_standard/<genome>.{gff3,proteins.faa} and
figures/inputs/reannotation_TXNL1.b.paralogs.{tsv,fa} (make txnl1-paralogs; the TXNL1.a protein there is the part
aligned to human TXNL1, residues given by the note for a merged model). Output: results/tables/stage2_txnl1_introns.tsv.
"""
import argparse, csv, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from txnl1_paralogs import aligner, pair, mapped   # noqa: E402


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].split()[0]; out[name] = []
        elif name: out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def introns(cds, strand):
    """(residue, phase) of the coding introns of a model: coding segments (1-based inclusive, genome order)."""
    segs = sorted(cds) if strand == "+" else sorted(cds)[::-1]
    out, off = [], 0
    for x, y in segs[:-1]:
        off += y - x + 1
        out.append((off // 3 + 1, off % 3))
    return out


def gold_cds(path, gene):
    cds, strand = [], "+"
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS" and re.search(rf"Parent=rna-[^_;]+_{re.escape(gene)}(;|$)", f[8]):
            cds.append((int(f[3]), int(f[4]))); strand = f[6]
    return cds, strand


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gene", default="TXNL1.b")
    ap.add_argument("--paralogs", default="figures/inputs/reannotation_TXNL1.b.paralogs")
    ap.add_argument("--out", default="results/tables/stage2_txnl1_introns.tsv")
    a = ap.parse_args()
    par = {}
    for r in csv.DictReader(open(f"{a.paralogs}.tsv"), delimiter="\t"):
        par.setdefault(r["assembly"], {"name": r["name"]})[r["field"]] = r["value"]
    pfa = read_fasta(f"{a.paralogs}.fa")
    al = aligner()
    rows = []
    for g in a.genomes:
        gold = next(s for k, s in read_fasta(f"{a.gold}/{g}.proteins.faa").items() if k.split("|")[0] == a.gene)
        cb, sb = gold_cds(f"{a.gold}/{g}.gff3", a.gene)
        p = par[g]
        ca = [tuple(map(int, x.split("-"))) for x in p["cds"].split(";")]
        m = re.search(r"drawn: residues (\d+)-", p.get("note", ""))
        shift = int(m.group(1)) - 1 if m else 0          # merged model: residues before the part drawn
        A = pfa[f"paralog:{g}"]
        blocks = pair(al, gold, A)
        ib = introns(cb, sb)
        ia = [(r - shift, ph) for r, ph in introns(ca, p["strand"]) if r - shift >= 1]
        on_a = {}
        for r, ph in ib:
            pos, _ = mapped(blocks, r - 1, 0, len(A))
            on_a[(r, ph)] = (pos + 1) if pos is not None else None
        shared_a = {(on_a[k], k[1]) for k in on_a if on_a[k] is not None}
        for r, ph in ib:
            rows.append([g, a.gene, r, ph, on_a[(r, ph)] or "", "yes" if (on_a[(r, ph)], ph) in set(ia) else "no"])
        for r, ph in ia:
            rows.append([g, p.get("label", "TXNL1.a"), r, ph, r, "yes" if (r, ph) in shared_a else "no"])
    with open(a.out, "w") as o:
        o.write("genome\tgene\tresidue\tphase\tresidue_on_TXNL1.a\tshared\n")
        for r in rows: o.write("\t".join(str(x) for x in r) + "\n")
    for g in a.genomes:
        sb = [r for r in rows if r[0] == g and r[1] == a.gene]
        sa = [r for r in rows if r[0] == g and r[1] != a.gene]
        print(f"{g}: {a.gene} {sum(r[5] == 'yes' for r in sb)} of {len(sb)} introns shared; "
              f"TXNL1.a {len(sa)} introns ({sum(r[5] == 'yes' for r in sa)} shared)")


if __name__ == "__main__":
    main()
