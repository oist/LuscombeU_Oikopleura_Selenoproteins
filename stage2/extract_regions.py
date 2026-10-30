#!/usr/bin/env python3
"""Genome regions on the strand of a gene, written to a small FASTA file so that they can be studied without the
genomes (2026-10-31; first use: the SELENOK loci that the UGA/SECIS screens found and Selenoprofiles did not).

Input (--loci, tab-separated, with a header): name, genome, contig, strand, position (1-based forward-strand
coordinate of a reference base, such as the leftmost base of a selenocysteine codon, which on the minus strand is the
last base of the codon on the gene strand), before, after (nucleotides kept
upstream and downstream of it on the gene strand), note. Output: FASTA records named
<name>|<genome>|<contig>:<start>-<end>|<strand>|offset=<0-based offset of the reference base in the record>, the
sequence read on the gene strand (start and end: 1-based, inclusive, forward strand; clipped to the contig).
"""
import argparse, csv, gzip

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path, keep):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in keep: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            elif name in keep: buf.append(l.strip())
    if name in keep: yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    loci = list(csv.DictReader(open(a.loci), delimiter="\t"))
    with open(a.out, "w") as o:
        for g in dict.fromkeys(r["genome"] for r in loci):
            rs = [r for r in loci if r["genome"] == g]
            contigs = dict(read_fasta(sheet[g]["fasta"], {r["contig"] for r in rs}))
            for r in rs:
                seq = contigs.get(r["contig"])
                if seq is None:
                    print(f"{r['name']}: contig {r['contig']} not found in {g}"); continue
                p, st, b, f = int(r["position"]), r["strand"], int(r["before"]), int(r["after"])
                s0 = max(1, p - b if st == "+" else p - f)
                e0 = min(len(seq), p + f if st == "+" else p + b)
                sub = seq[s0 - 1:e0].upper()
                if st == "-": sub = sub.translate(COMP)[::-1]
                off = p - s0 if st == "+" else e0 - p
                o.write(f">{r['name']}|{g}|{r['contig']}:{s0}-{e0}|{st}|offset={off}\n{sub}\n")
                print(f"{r['name']}: {len(sub)} nt, reference base at {off}")


if __name__ == "__main__":
    main()
