#!/usr/bin/env python3
"""Spliced alignment of transcripts or ESTs to a genome (minimap2 splice preset through mappy). Writes one row per
primary alignment: query, qlen, qstart, qend, strand, contig, tstart, tend, mapq, blocks (1-based exons
"start-end,..."). Carriage returns and end-of-file characters in the input are ignored.
Sequences are processed in blocks; finished blocks are kept in <output>.blocks/ so an interrupted run resumes."""
import argparse, gzip, os, shutil
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
    ap.add_argument("--block-size", type=int, default=20000)
    a = ap.parse_args()
    bdir = a.output + ".blocks"; os.makedirs(bdir, exist_ok=True)
    al = None
    def flush(i, block):
        nonlocal al
        path = f"{bdir}/{i:05d}.tsv"
        if os.path.exists(path): return
        if al is None: al = mp.Aligner(a.genome, preset="splice", k=14)
        with open(path + ".part", "w") as o:
            for name, seq in block:
                for h in al.map(seq):
                    if not h.is_primary: continue
                    blocks, t, bs = [], h.r_st, h.r_st
                    for length, op in h.cigar:
                        if op in (0, 2, 7, 8): t += length
                        elif op == 3:
                            blocks.append(f"{bs + 1}-{t}"); t += length; bs = t
                    blocks.append(f"{bs + 1}-{t}")
                    o.write(f"{name}\t{len(seq)}\t{h.q_st}\t{h.q_en}\t{'+' if h.strand == 1 else '-'}\t{h.ctg}\t{h.r_st + 1}\t{h.r_en}\t{h.mapq}\t{','.join(blocks)}\n")
        os.rename(path + ".part", path)
    block, i = [], 0
    for rec in read_fasta(a.transcripts):
        block.append(rec)
        if len(block) == a.block_size:
            flush(i, block); block, i = [], i + 1
    if block: flush(i, block)
    with open(a.output + ".part", "w") as o:
        o.write("query\tqlen\tqstart\tqend\tstrand\tcontig\ttstart\ttend\tmapq\tblocks\n")
        for f in sorted(os.listdir(bdir)):
            if f.endswith(".tsv"):
                with open(f"{bdir}/{f}") as fh: shutil.copyfileobj(fh, o)
    os.rename(a.output + ".part", a.output); shutil.rmtree(bdir)

if __name__ == "__main__":
    main()
