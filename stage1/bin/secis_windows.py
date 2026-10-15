#!/usr/bin/env python3
"""Merge Sec codons from Selenoprofiles (GTF Selenocysteine features of predictions labelled selenocysteine) and
from miniprot (sec_check.py), and write the 3 kb downstream of each Sec codon (strand-aware) for SECIS searches.

Outputs: <prefix>.sec_sites.tsv (contig, position of the first base of the TGA on the + strand, strand, families,
sources) and <prefix>.windows.fa (one sequence per site, named s<row>)."""
import argparse, csv, gzip, re
from Bio.Seq import Seq

def read_sites(gtf, check):
    sites = []
    for l in open(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "Selenocysteine":
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            if gid.endswith("selenocysteine"):
                sites.append((f[0], int(f[3]), f[6], gid.split(":")[-1].split(".")[0], "selenoprofiles"))
    for r in csv.DictReader(open(check), delimiter="\t"):
        sites.append((r["contig"], int(r["start"]), r["strand"], r["family"], "miniprot"))
    merged = []
    for c, p, st, fam, src in sites:
        hit = [m for m in merged if m["contig"] == c and m["strand"] == st and abs(m["pos"] - p) <= 30]
        if hit:
            hit[0]["families"].add(fam); hit[0]["sources"].add(src)
        else:
            merged.append({"contig": c, "pos": p, "strand": st, "families": {fam}, "sources": {src}})
    return merged

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True); ap.add_argument("--gtf", required=True)
    ap.add_argument("--sec-check", required=True); ap.add_argument("--prefix", required=True)
    ap.add_argument("--window", type=int, default=3000)
    a = ap.parse_args()
    sites = read_sites(a.gtf, a.sec_check)
    want = {s["contig"] for s in sites}
    seqs, name, buf = {}, None, []
    with (gzip.open(a.genome, "rt") if a.genome.endswith(".gz") else open(a.genome)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in want: seqs[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            elif name in want: buf.append(l.strip())
    if name in want: seqs[name] = "".join(buf).upper()
    with open(f"{a.prefix}.sec_sites.tsv", "w") as t, open(f"{a.prefix}.windows.fa", "w") as w:
        t.write("site\tcontig\tsec_codon_start\tstrand\tfamilies\tsources\n")
        for i, s in enumerate(sites):
            q, p = seqs[s["contig"]], s["pos"]
            win = q[p - 1:p - 1 + a.window] if s["strand"] == "+" else str(Seq(q[max(0, p + 2 - a.window):p + 2]).reverse_complement())
            t.write(f"s{i}\t{s['contig']}\t{p}\t{s['strand']}\t{','.join(sorted(s['families']))}\t{','.join(sorted(s['sources']))}\n")
            w.write(f">s{i}\n{win}\n")

if __name__ == "__main__":
    main()
