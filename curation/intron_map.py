#!/usr/bin/env python3
"""Intron positions of one gene in the appendicularian models, on a common protein alignment.

For each gold-standard model of the gene (--gold-dirs, <NAME>.gff3 and <NAME>.proteins.faa; the assemblies of
--assemblies), the coding exons give the position of each intron in the protein (the residue after which it lies, x,
and its phase: 0 between codons, 1 or 2 inside codon x + 1; when the protein was taken from a transcript because the
assembly has an indel in the coding sequence, the positions come from the exonerate alignment of that protein to the
region, as below), its length and its splice sites, read in the region of
make machinery-regions (results/gene_regions/<gene>/<NAME>.fa, first record): the placement of the model and the
equivalent placements (the intron moved across the bases repeated at its two ends, which give the same spliced
sequence; as figures/intron_placements.py). For the assemblies of other individuals (--others NAME=SOURCE), the
gold-standard protein of SOURCE is aligned to the region of NAME with exonerate (protein2genome, best alignment, introns
of 20 nt or more), so that their introns are found independently of the source's exon chain; each of these introns is
moved to its equivalent placement of an intron class nearest to the alignment, when there is one; the position of
an intron in the protein is read from the alignment (frameshifts of the assembly do not shift it). The proteins of all
models and the reference proteins (--references, FASTA: the human protein first) are aligned with MAFFT L-INS-i; each
intron is placed at the alignment column of residue x + 1 (the residue that follows or contains it), and reported with
the residue of the first reference at or before that column. Introns of different models at the same column and phase
are candidates for one intron; a column a few residues away is another intron unless the alignment is uncertain there.
Outputs: --out (TSV: gene, model, source, intron, length, residue, phase, column, reference_residue, sites,
placements; "*" marks the placements of the O. dioica intron classes GT-AG, GC-AG, GA-AG and GG-AG) and, on stdout,
the introns grouped by column and phase with the models that have them.
"""
import argparse, collections, os, re, subprocess, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from intron_placements import region, gold_exons, placements, rc, translate   # noqa: E402

CLASSES = {"GT-AG", "GC-AG", "GA-AG", "GG-AG"}


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].strip(); out[name] = []
        elif name: out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def introns_of(exons, gs, to_class=False):
    """Introns of a chain of gene-strand intervals [a, b) of gs: (residue, phase, length, site, placements). With
    to_class, each intron is first moved to its equivalent placement of an intron class nearest to the alignment
    (for exonerate chains, whose placement among the equivalent ones is arbitrary)."""
    ex = [list(e) for e in exons]
    if to_class:
        for i in range(len(ex) - 1):
            cls = [k for k, site in placements(gs, ex[i][1], ex[i + 1][0]) if site in CLASSES]
            if cls:
                k = min(cls, key=lambda k: (abs(k), k))
                ex[i][1] += k; ex[i + 1][0] += k
    out, cds = [], 0
    for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
        cds += b1 - a1
        pl = placements(gs, b1, a2)
        site = gs[b1:b1 + 2] + "-" + gs[a2 - 2:a2]
        out.append((cds // 3, cds % 3, a2 - b1, site, pl))
    return out


def exonerate_introns(protein, r, to_class=True):
    """Introns of the best exonerate protein2genome alignment of protein on the region (gene strand), with their
    position in the protein read from the alignment (vulgar: residues of the query before the intron, and the bases of
    a codon split by it), so that a frameshift of the assembly inside an exon does not shift the phases of the introns
    after it: (residue, phase, length, site, placements), as introns_of; with to_class, each intron is first moved to
    its equivalent placement of an intron class nearest to the alignment."""
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{protein}\n"); open(f"{d}/t.fa", "w").write(f">t\n{r['gs']}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--bestn", "1", "--showalignment", "no",
                              "--showvulgar", "yes", "--minintron", "20", "--maxintron", "20000",
                              f"{d}/q.fa", f"{d}/t.fa"], capture_output=True, text=True, check=True).stdout
    v = next((l.split() for l in out.splitlines() if l.startswith("vulgar:")), None)
    if v is None: return []
    if v[8] != "+": raise SystemExit("exonerate alignment not on the gene strand of the region")
    q, t, split, res = int(v[2]), int(v[6]), 0, []
    ops = v[10:]
    for i in range(0, len(ops), 3):
        op, ql, tl = ops[i], int(ops[i + 1]), int(ops[i + 2])
        if op == "5": donor, cds_pos = t, 3 * q + split
        if op == "3": res.append((cds_pos, donor, t + tl))
        if op == "S":
            if split == 0 and ql == 0: split = tl          # bases of a codon before the intron
            else: q += ql; split = 0                        # the rest of the codon after it
        elif op in "MCG":
            q += ql
        t += tl
    gs, out = r["gs"], []
    for cds_pos, b1, a2 in res:
        k = 0
        if to_class:
            cls = [kk for kk, site in placements(gs, b1, a2) if site in CLASSES]
            if cls: k = min(cls, key=lambda kk: (abs(kk), kk))
        b, a = b1 + k, a2 + k
        out.append(((cds_pos + k) // 3, (cds_pos + k) % 3, a - b, gs[b:b + 2] + "-" + gs[a - 2:a], placements(gs, b, a)))
    return out


def exonerate_exons(protein, r):
    """Coding exons (gene-strand intervals of r["gs"]) of the best exonerate protein2genome alignment of protein."""
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{protein}\n"); open(f"{d}/t.fa", "w").write(f">t\n{r['gs']}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--bestn", "1", "--showalignment", "no",
                              "--showvulgar", "no", "--showtargetgff", "yes", "--minintron", "20", "--maxintron", "20000",
                              f"{d}/q.fa", f"{d}/t.fa"], capture_output=True, text=True, check=True).stdout
    rows = [l.split("\t") for l in out.splitlines() if l.count("\t") >= 8 and l.split("\t")[2] == "cds"]
    if not rows: return []
    if any(f[6] != "+" for f in rows): raise SystemExit("exonerate alignment not on the gene strand of the region")
    return sorted((int(f[3]) - 1, int(f[4])) for f in rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--gold-dirs", nargs="+", default=["curation/gold_standard_machinery", "curation/gold_standard"])
    ap.add_argument("--assemblies", nargs="+", default=["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"])
    ap.add_argument("--others", nargs="*", default=["KUM=OKI", "AOM=OSKA", "Bar2=OdB3"])
    ap.add_argument("--regions", default="results/gene_regions")
    ap.add_argument("--references", help="FASTA of reference proteins, the first one numbering the columns "
                                         "(default: figures/inputs/reannotation_<gene>.fa, records not from the assemblies)")
    ap.add_argument("--out")
    a = ap.parse_args()
    gene = a.gene
    out_path = a.out or f"results/tables/intron_map_{gene}.tsv"
    models, proteins = {}, {}
    for g in a.assemblies:
        rpath = f"{a.regions}/{gene}/{g}.fa"
        gdir = next((d for d in a.gold_dirs if os.path.exists(f"{d}/{g}.gff3") and
                     any(n.startswith(gene + "|") for n in read_fasta(f"{d}/{g}.proteins.faa"))), None)
        if gdir is None or not os.path.exists(rpath): continue
        r = region(rpath)
        ex = gold_exons(f"{gdir}/{g}.gff3", gene, r)
        prot = next(s for n, s in read_fasta(f"{gdir}/{g}.proteins.faa").items() if n.startswith(gene + "|")).rstrip("*")
        if translate("".join(r["gs"][x:y] for x, y in ex)).rstrip("*") == prot:
            models[g] = ("gold standard", introns_of(ex, r["gs"]))
        else:   # protein from a transcript (the assembly has an indel in the coding sequence): positions from exonerate
            models[g] = ("gold standard, protein from a transcript, positions by exonerate", exonerate_introns(prot, r))
        proteins[g] = prot
    for spec in a.others:
        name, src = spec.split("=")
        rpath = f"{a.regions}/{gene}/{name}.fa"
        if src not in proteins or not os.path.exists(rpath): continue
        r = region(rpath)
        ex = exonerate_exons(proteins[src], r)
        if not ex: continue
        prot = translate("".join(r["gs"][x:y] for x, y in ex))
        models[name] = (f"exonerate of the {src} protein", exonerate_introns(proteins[src], r))
        proteins[name] = prot.rstrip("*")
    refs = a.references or f"figures/inputs/reannotation_{gene}.fa"
    ref = [(n.split()[0].replace(" ", "_"), s.replace("*", "")) for n, s in read_fasta(refs).items()
           if not re.match(r"(OKI|OSKA|OdB3)_", n)]
    if not ref: raise SystemExit(f"no reference protein in {refs}")
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/in.fa", "w") as o:
            for n, s in ref: o.write(f">ref_{n}\n{s}\n")
            for n, s in proteins.items(): o.write(f">{n}\n{s.replace('U', 'C')}\n")
        aln_txt = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", f"{d}/in.fa"],
                                 capture_output=True, text=True, check=True).stdout
    aln = {n: s.replace("\n", "") for n, s in (b.split("\n", 1) for b in aln_txt.split(">")[1:])}
    first = aln[f"ref_{ref[0][0]}"]

    def column(name, res):          # 1-based alignment column of residue res (1-based) of a model
        k = 0
        for j, ch in enumerate(aln[name]):
            if ch != "-":
                k += 1
                if k == res: return j + 1
        return None

    ref_res = lambda col: sum(ch != "-" for ch in first[:col])
    rows, groups = [], collections.defaultdict(list)
    for name, (source, ints) in models.items():
        for i, (res, ph, length, site, pl) in enumerate(ints, 1):
            col = column(name, res + 1)
            plist = ";".join(f"{k}:{s}{'*' if s in CLASSES else ''}" for k, s in pl)
            rows.append([gene, name, source, str(i), str(length), str(res), str(ph), str(col), str(ref_res(col)), site, plist])
            groups[(col, ph)].append(f"{name}:{i}({length} nt, {site})")
    with open(out_path, "w") as o:
        o.write("gene\tmodel\tsource\tintron\tlength\tresidue\tphase\tcolumn\treference_residue\tsites\tplacements\n")
        for r_ in rows: o.write("\t".join(r_) + "\n")
    print(f"{out_path}: {len(rows)} introns of {len(models)} models; reference {ref[0][0]}")
    for (col, ph), members in sorted(groups.items()):
        print(f"column {col} ({ref[0][0]} residue {ref_res(col)}), phase {ph}: " + ", ".join(members))


if __name__ == "__main__":
    main()
