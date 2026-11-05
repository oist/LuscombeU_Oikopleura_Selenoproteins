#!/usr/bin/env python3
"""Splice junctions of an exon chain rescanned against reference proteins.

For a gene and an assembly, the chain is read from --overrides (the curator's row, 1-based start-end of the coding
exons on the contig) or given with --chain, and placed on the region of make machinery-regions
(results/gene_regions/<gene>/<NAME>.fa, first record; gene strand). For each intron, every donor and acceptor within
--window nt of its ends that keeps the downstream reading frame is tried; the whole chain is translated (assembly gaps
read as X) and scored by the sum of its global alignment scores to the reference proteins (BLOSUM62, gap opening 11,
extension 1, free end gaps; X scored as A); pairs whose protein has an internal stop are dropped. The output lists,
per intron, the best pairs (equivalent placements of one intron, which give the same protein, score alike), and the
best pair whose sites are of an O. dioica intron class (GT-AG, GC-AG, GA-AG, GG-AG). It also reports, per intron,
the residue after which it lies and its phase, and the protein around the junction (stdout; --out, TSV: gene,
assembly, contig, strand, intron, residue, phase, sites, length, junction, chain_score, rank, score, candidate_sites,
candidate_length, donor_shift, acceptor_shift, candidate_intron, intron_class; shifts in nt along the gene, the
intron in contig coordinates). This is the evidence of the
curator's choice where no transcript tests the junctions: a site outside the classes is kept when the class sites
nearby add or remove residues against all references, and a class site is taken when the scores do not separate them.
References: the gold-standard proteins of --references (assemblies; curation/gold_standard_machinery).
"""
import argparse, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from intron_placements import region, translate   # noqa: E402
from Bio import Align                              # noqa: E402
from Bio.Align import substitution_matrices        # noqa: E402

CLASSES = {"GT-AG", "GC-AG", "GA-AG", "GG-AG"}


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].strip(); out[name] = []
        elif name: out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--assemblies", nargs="+", required=True)
    ap.add_argument("--overrides", default="curation/machinery_overrides.tsv")
    ap.add_argument("--chain", help="1-based start-end of the coding exons on the contig, separated by ';' "
                                    "(default: the curator's row of --overrides)")
    ap.add_argument("--regions", default="results/gene_regions")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--references", nargs="+", default=["OSKA", "OdB3"])
    ap.add_argument("--window", type=int, default=30)
    ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--out", help="TSV of the candidates shown")
    a = ap.parse_args()
    rows = []
    for asm in a.assemblies:
        rows += rescan(a, asm)
    if a.out:
        with open(a.out, "w") as o:
            o.write("\t".join(["gene", "assembly", "contig", "strand", "intron", "residue", "phase", "sites", "length",
                               "junction", "chain_score", "rank", "score", "candidate_sites", "candidate_length",
                               "donor_shift", "acceptor_shift", "candidate_intron", "intron_class"]) + "\n")
            for row in rows: o.write("\t".join(map(str, row)) + "\n")


def rescan(a, asm):
    out = []
    r = region(f"{a.regions}/{a.gene}/{asm}.fa")
    gs = r["gs"]
    chain = a.chain or next(l.split("\t")[2] for l in open(a.overrides)
                            if l.startswith(f"{asm}\t{a.gene}\t"))
    ex = [(int(x), int(y)) for x, y in (b.split("-") for b in chain.split(";"))]
    if r["strand"] == "+": blocks = sorted((x - r["start"], y - r["start"] + 1) for x, y in ex)
    else: blocks = sorted((r["end"] - y, r["end"] - x + 1) for x, y in ex)
    if blocks[0][0] < 0 or blocks[-1][1] > len(gs): raise SystemExit("chain outside the region")
    refs = []
    for g in a.references:
        prot = next((s for n, s in read_fasta(f"{a.gold_dir}/{g}.proteins.faa").items() if n.startswith(a.gene + "|")), None)
        if prot: refs.append(prot.rstrip("*").replace("U", "C"))
    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    al.end_deletion_score = al.end_insertion_score = 0

    def score(c):
        p = translate("".join(gs[x:y] for x, y in c))
        body = p[:-1] if p.endswith("*") else p
        if "*" in body: return None, p
        return sum(al.score(body.replace("X", "A"), ref) for ref in refs), p

    s0, p0 = score(blocks)
    print(f"{a.gene} {asm} {r['contig']} {r['strand']}: {len(blocks)} exons, {len(p0.rstrip('*'))} residues, "
          f"score {s0} against {', '.join(a.references)}")
    cds = 0
    for i in range(len(blocks) - 1):
        (a1, b1), (a2, b2) = blocks[i], blocks[i + 1]
        cds += b1 - a1
        k = cds // 3
        print(f"intron {i + 1}: residue {k}, phase {cds % 3}, {gs[b1:b1 + 2]}-{gs[a2 - 2:a2]}, {a2 - b1} nt, "
              f"{p0[max(0, k - 8):k]}|{p0[k:k + 8]}")
        res = []
        for d in range(b1 - a.window, b1 + a.window + 1):
            for acc in range(a2 - a.window, a2 + a.window + 1):
                if d <= a1 or acc >= b2 or acc - d < 20 or (d - b1) % 3 != (acc - a2) % 3: continue
                s, _ = score(blocks[:i] + [(a1, d), (acc, b2)] + blocks[i + 2:])
                if s is None: continue
                res.append((s, d, acc, gs[d:d + 2] + "-" + gs[acc - 2:acc], acc - d))
        res.sort(key=lambda t: (-t[0], abs(t[1] - b1) + abs(t[2] - a2)))
        shown = res[:a.top]
        cls = next((t for t in res if t[3] in CLASSES), None)
        if cls and cls not in shown: shown.append(cls)
        pos = lambda q: (r["start"] + q) if r["strand"] == "+" else (r["end"] - q)
        for rank, t in enumerate(shown, 1):
            lo, hi = min(pos(t[1]), pos(t[2] - 1)), max(pos(t[1]), pos(t[2] - 1))
            out.append([a.gene, asm, r["contig"], r["strand"], i + 1, k, cds % 3, f"{gs[b1:b1 + 2]}-{gs[a2 - 2:a2]}",
                        a2 - b1, f"{p0[max(0, k - 8):k]}|{p0[k:k + 8]}", f"{s0:.0f}", rank if t in res[:a.top] else "class",
                        f"{t[0]:.0f}", t[3], t[4], t[1] - b1, t[2] - a2, f"{lo}-{hi}", "yes" if t[3] in CLASSES else "no"])
            print(f"  score {t[0]:.0f}: {t[3]}, {t[4]} nt (donor {t[1] - b1:+d}, acceptor {t[2] - a2:+d}; "
                  f"intron {min(pos(t[1]), pos(t[2] - 1))}-{max(pos(t[1]), pos(t[2] - 1))})"
                  + (", intron class" if t[3] in CLASSES else ""))
    return out


if __name__ == "__main__":
    main()
