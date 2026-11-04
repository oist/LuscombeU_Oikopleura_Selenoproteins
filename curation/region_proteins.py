#!/usr/bin/env python3
"""Proteins of a gene in the assemblies without a gold-standard model, from the regions of curation/gene_regions.py.

For each region of results/gene_regions/<gene>/regions.tsv whose source is lift:SOURCE or miniprot:SOURCE (the same
locus in the assembly of another individual, or a protein alignment to an assembly without a model), the
gold-standard protein of SOURCE (<gold-dir>/<SOURCE>.proteins.faa) is aligned to the region with exonerate
(protein2genome, introns of --min-intron to --max-intron nt, best alignment only), and the coding sequence of the
alignment is translated. Assembly gaps (runs of N) inside the aligned span are counted. These proteins serve to judge
whether an assembly or a transcript is wrong where they differ (curation/gold_standard.py --other-proteins), since
another assembly of the same species gives the same protein except for polymorphisms.
Outputs in --out (default: the gene directory): proteins.faa (header <gene>|<genome>|<source>|exonerate, stop codons
inside the protein written *), and proteins.tsv: gene, genome, contig, region, source, query, query_span, identity,
similarity, exons, introns (lengths, gene order), sites (dinucleotides of each intron), internal_stops, gaps (N runs
inside the aligned span), protein_length.
"""
import argparse, csv, os, re, subprocess, tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].strip(); out[name] = []
        elif name: out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def gold_protein(path, gene):
    name = None
    for l in open(path):
        if l.startswith(">"): name = l[1:].split("|")[0]
        elif name == gene: return l.strip().rstrip("*").replace("*", "X").replace("U", "X")
    return None


def translate(s):
    s = s.upper()
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def exonerate(query, target, a):
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.faa", "w").write(f">q\n{query}\n")
        open(f"{d}/t.fa", "w").write(f">t\n{target}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--refine", "region", "--bestn", "1",
                              "--minintron", str(a.min_intron), "--maxintron", str(a.max_intron), "--percent", "30",
                              "--showalignment", "no", "--showvulgar", "no", "--showtargetgff", "yes",
                              "--ryo", ">cds %qab %qae %pi %ps\n%tcs\n", f"{d}/q.faa", f"{d}/t.fa"],
                             capture_output=True, text=True).stdout
    m = re.search(r"^>cds (\d+) (\d+) ([\d.]+) ([\d.]+)\n([ACGTNacgtn\n]+)", out, re.M)
    if not m: return None
    gff = [l.split("\t") for l in out.splitlines() if l.count("\t") >= 8 and not l.startswith("#")]
    exons = sorted((int(f[3]), int(f[4])) for f in gff if f[2] == "exon")
    strand = next((f[6] for f in gff if f[2] == "gene"), "+")
    introns = [(b1 + 1, a2 - 1) for (a1, b1), (a2, b2) in zip(exons, exons[1:])]
    if strand == "-": introns = introns[::-1]

    def site(lo, hi):
        s = target[lo - 1:hi]
        if strand == "-": s = s.translate(COMP)[::-1]
        return f"{s[:2]}-{s[-2:]}".upper()
    span = (exons[0][0], exons[-1][1]) if exons else (0, 0)
    return dict(qspan=f"{int(m.group(1)) + 1}-{m.group(2)}", identity=m.group(3), similarity=m.group(4),
                cds=m.group(5).replace("\n", ""), exons=len(exons), introns=",".join(str(hi - lo + 1) for lo, hi in introns),
                sites=",".join(site(lo, hi) for lo, hi in introns),
                gaps=len(re.findall(r"N+", target[span[0] - 1:span[1]].upper())))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", default=["PSTK", "SBP2", "SEPHS_b", "SEPSECS", "EEFSEC"])
    ap.add_argument("--regions", default="results/gene_regions")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--min-intron", type=int, default=20)
    ap.add_argument("--max-intron", type=int, default=25000)
    a = ap.parse_args()
    for gene in a.genes:
        d = f"{a.regions}/{gene}"
        if not os.path.exists(f"{d}/regions.tsv"): print(gene, "no regions"); continue
        rows, fasta = [], []
        for r in csv.DictReader(open(f"{d}/regions.tsv"), delimiter="\t"):
            m = re.match(r"^(lift|miniprot):([^:]+)", r["source"])
            if not m: continue
            src = m.group(2)
            q = gold_protein(f"{a.gold_dir}/{src}.proteins.faa", gene)
            seqs = read_fasta(f"{d}/{r['genome']}.fa")
            key = next((k for k in seqs if k.split("|")[1] == f"{r['contig']}:{r['start']}-{r['end']}" and
                        k.endswith("|" + r["source"])), None)
            if q is None or key is None: continue
            x = exonerate(q, seqs[key], a)
            region = f"{r['contig']}:{r['start']}-{r['end']}"
            if x is None:
                rows.append(dict(gene=gene, genome=r["genome"], contig=r["contig"], region=region, source=r["source"],
                                 query=f"{src}_{gene}", query_span="", identity="", similarity="", exons=0, introns="",
                                 sites="", internal_stops="", gaps="", protein_length=0)); continue
            prot = translate(x["cds"])
            stops = prot.rstrip("*").count("*")
            prot = prot.rstrip("*")
            fasta.append(f">{gene}|{r['genome']}|{r['source']}|exonerate\n{prot}\n")
            rows.append(dict(gene=gene, genome=r["genome"], contig=r["contig"], region=region, source=r["source"],
                             query=f"{src}_{gene}", query_span=x["qspan"], identity=x["identity"],
                             similarity=x["similarity"], exons=x["exons"], introns=x["introns"], sites=x["sites"],
                             internal_stops=stops, gaps=x["gaps"], protein_length=len(prot)))
        cols = ["gene", "genome", "contig", "region", "source", "query", "query_span", "identity", "similarity", "exons",
                "introns", "sites", "internal_stops", "gaps", "protein_length"]
        with open(f"{d}/proteins.tsv", "w") as o:
            o.write("\t".join(cols) + "\n")
            for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
        open(f"{d}/proteins.faa", "w").write("".join(fasta))
        print(gene, len(rows), "regions aligned")


if __name__ == "__main__":
    main()
