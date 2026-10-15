#!/usr/bin/env python3
"""Spliced alignment of transcripts or ESTs to a genome (minimap2 splice preset through mappy). Writes one row per
primary alignment: query, qlen, qstart, qend, strand, contig, tstart, tend, mapq, blocks (1-based exons
"start-end,..."). Carriage returns and end-of-file characters in the input are ignored."""
import argparse, gzip
import mappy as mp

def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            l = l.replace("\r", "").replace("\x1a", "").rstrip("\n")
            if l.startswith(">"):
                if name: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            elif l: buf.append(l)
    if name: yield name, "".join(buf)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True); ap.add_argument("--transcripts", required=True); ap.add_argument("--output", required=True)
    a = ap.parse_args()
    al = mp.Aligner(a.genome, preset="splice", k=14)
    with open(a.output, "w") as o:
        o.write("query\tqlen\tqstart\tqend\tstrand\tcontig\ttstart\ttend\tmapq\tblocks\n")
        for name, seq in read_fasta(a.transcripts):
            for h in al.map(seq):
                if not h.is_primary: continue
                blocks, t, bs = [], h.r_st, h.r_st
                for length, op in h.cigar:
                    if op in (0, 2, 7, 8): t += length
                    elif op == 3:
                        blocks.append(f"{bs + 1}-{t}"); t += length; bs = t
                blocks.append(f"{bs + 1}-{t}")
                o.write(f"{name}\t{len(seq)}\t{h.q_st}\t{h.q_en}\t{'+' if h.strand == 1 else '-'}\t{h.ctg}\t{h.r_st + 1}\t{h.r_en}\t{h.mapq}\t{','.join(blocks)}\n")

if __name__ == "__main__":
    main()
