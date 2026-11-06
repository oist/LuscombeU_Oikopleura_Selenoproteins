#!/usr/bin/env python3
"""Coding divergence between individuals and species of O. dioica, gene by gene.

The gold-standard protein of each species (OKI2018, OSKA2016, OdB3; --gold-dir) is aligned with exonerate
(protein2genome, best alignment) to every record of the region of its own assembly and of the other individual of the
species (--pairs, default KUM-M3-7f, AOM-5-5f, Bar2_p4; results/gene_regions/<gene>/<NAME>.fa; the record with the
best score is used, so that a lift that failed on another contig is skipped). The alignment gives, for each residue
of the protein, the codon of the assembly (codons split by an intron are joined; codons at a frameshift or in a gap of
the alignment are left out), so that the two individuals of a species are compared codon by codon through the same
protein, whatever their exon chains and indels; between species, the codons are matched through a global alignment of
the two gold-standard proteins (BLOSUM62). For each pair of assemblies, over the residues with a complete codon (no N)
in both: the fraction of codons that differ at the third position, and the fraction that encode different amino
acids. Output: --out (TSV: gene, a, b, codons, third_position, amino_acid).
"""
import argparse, itertools, os, re, subprocess, sys, tempfile
from Bio import Align
from Bio.Align import substitution_matrices

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from intron_placements import rc, translate   # noqa: E402


def records(path):
    out, head, seq = [], None, []
    for l in open(path):
        if l.startswith(">"):
            if head: out.append((head, "".join(seq).upper()))
            head, seq = l[1:].strip(), []
        else: seq.append(l.strip())
    if head: out.append((head, "".join(seq).upper()))
    res = []
    for h, s in out:
        m = re.match(r"[^|]+\|([^:]+):(\d+)-(\d+)\|([+\-])", h)
        res.append(s if m.group(4) == "+" else rc(s))
    return res


def codons_of(protein, gs):
    """{reference residue (0-based): codon of gs} from the best exonerate alignment, and its score."""
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{protein}\n"); open(f"{d}/t.fa", "w").write(f">t\n{gs}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--bestn", "1", "--showalignment", "no",
                              "--showvulgar", "yes", "--minintron", "20", "--maxintron", "20000",
                              f"{d}/q.fa", f"{d}/t.fa"], capture_output=True, text=True, check=True).stdout
    v = next((l.split() for l in out.splitlines() if l.startswith("vulgar:")), None)
    if v is None or v[8] != "+": return {}, 0
    q, t, part, cod = int(v[2]), int(v[6]), "", {}
    ops = v[10:]
    for i in range(0, len(ops), 3):
        op, ql, tl = ops[i], int(ops[i + 1]), int(ops[i + 2])
        if op in "MC":
            for k in range(ql): cod[q + k] = gs[t + 3 * k:t + 3 * k + 3]
        elif op == "S":
            if ql == 0: part = gs[t:t + tl]                   # bases of a codon before an intron
            else:
                if len(part) + tl == 3: cod[q] = part + gs[t:t + tl]
                part = ""
        if op in "MCGS": q += ql
        t += tl
    return cod, int(v[9])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--pairs", nargs="+", default=["OKI=OKI", "KUM=OKI", "OSKA=OSKA", "AOM=OSKA", "OdB3=OdB3",
                                                     "Bar2=OdB3"], help="assembly=species whose protein is aligned")
    ap.add_argument("--regions", default="results/gene_regions")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--out", default="results/tables/cds_divergence.tsv")
    a = ap.parse_args()
    pal = Align.PairwiseAligner(mode="global")
    pal.substitution_matrix = substitution_matrices.load("BLOSUM62")
    pal.open_gap_score, pal.extend_gap_score = -11, -1
    species = dict(x.split("=") for x in a.pairs)
    rows = []
    for gene in a.genes:
        prot = {}
        for sp in sorted(set(species.values())):
            for h, s in (l.split("\n", 1) for l in open(f"{a.gold_dir}/{sp}.proteins.faa").read().split(">")[1:]):
                if h.startswith(gene + "|"): prot[sp] = s.replace("\n", "").rstrip("*").replace("U", "C")
        cods = {}
        for g, sp in species.items():
            path = f"{a.regions}/{gene}/{g}.fa"
            if sp not in prot or not os.path.exists(path): continue
            best = max((codons_of(prot[sp], gs) for gs in records(path)), key=lambda x: x[1])
            if best[0]: cods[g] = best[0]

        def index_map(sx, sy):   # residue of the protein of sx -> residue of the protein of sy
            if sx == sy: return {k: k for k in range(len(prot[sx]))}
            al = pal.align(prot[sx], prot[sy])[0]
            return {i: j for (xs, xe), (ys, ye) in zip(*al.aligned) for i, j in zip(range(xs, xe), range(ys, ye))}
        for x, y in itertools.combinations([g for g in species if g in cods], 2):
            im = index_map(species[x], species[y])
            common = [(k, im[k]) for k in cods[x] if k in im and im[k] in cods[y]
                      and len(cods[x][k]) == 3 and len(cods[y][im[k]]) == 3 and "N" not in cods[x][k] + cods[y][im[k]]]
            if not common: continue
            d3 = sum(cods[x][k][2] != cods[y][j][2] for k, j in common)
            daa = sum(translate(cods[x][k]) != translate(cods[y][j]) for k, j in common)
            rows.append([gene, x, y, len(common), f"{d3 / len(common):.3f}", f"{daa / len(common):.3f}"])
            print(f"{gene:8s} {x:5s} {y:5s} {len(common):4d} codons: third position {100 * d3 / len(common):5.1f}%, "
                  f"amino acid {100 * daa / len(common):4.1f}%")
    with open(a.out, "w") as o:
        o.write("gene\ta\tb\tcodons\tthird_position\tamino_acid\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
