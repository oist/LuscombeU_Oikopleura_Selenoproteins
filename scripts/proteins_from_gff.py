#!/usr/bin/env python3
"""Protein sequences from a genome and its GFF3 annotation: one protein per CDS protein_id (RefSeq style), its CDS
segments joined in transcript order, translated from the phase of the first segment with the standard code. The
genome is read one sequence at a time, so memory stays low. Equivalent to `gffread -y` for RefSeq GFF3 files; used when
a genome's protein file was not downloaded.
Usage: proteins_from_gff.py genome.fna.gz annotation.gff3.gz > protein.faa
"""
import collections, gzip, re, sys

CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    genome, gff = sys.argv[1], sys.argv[2]
    cds = collections.defaultdict(list)       # protein_id -> [(seqid, start, end, strand, phase)]
    by_seq = collections.defaultdict(set)
    for l in opener(gff):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        m = re.search(r"protein_id=([^;]+)", f[8])
        if not m: continue
        pid = m.group(1)
        cds[pid].append((f[0], int(f[3]), int(f[4]), f[6], int(f[7]) if f[7] in "012" else 0))
        by_seq[f[0]].add(pid)

    def emit(name, seq):
        for pid in sorted(by_seq.get(name, ())):
            segs = cds[pid]
            if len({s[0] for s in segs}) > 1: continue      # trans-spliced across sequences: skipped
            strand = segs[0][3]
            segs = sorted(segs, key=lambda s: s[1], reverse=(strand == "-"))
            nt = "".join(seq[s[1] - 1:s[2]] if strand == "+" else seq[s[1] - 1:s[2]].translate(COMP)[::-1]
                         for s in segs)
            nt = nt[segs[0][4]:].upper()
            prot = "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3)).rstrip("*")
            sys.stdout.write(f">{pid}\n{prot}\n")

    name, buf = None, []
    for l in opener(genome):
        if l.startswith(">"):
            if name is not None: emit(name, "".join(buf))
            name, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if name is not None: emit(name, "".join(buf))


if __name__ == "__main__":
    main()
