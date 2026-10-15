#!/usr/bin/env python3
"""For each Sec site of a genome, count the transcripts (or ESTs) whose aligned exons contain the whole Sec codon,
and report how far the longest of them extends 3' of the codon."""
import argparse, csv
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", required=True); ap.add_argument("--alignments", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True); ap.add_argument("--output", required=True)
    a = ap.parse_args()
    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    bycontig = defaultdict(list)
    for s in sites: bycontig[s["contig"]].append(s)
    counts = {(s["site"], lab): [0, 0] for s in sites for lab in a.labels}
    for path, lab in zip(a.alignments, a.labels):
        for r in csv.DictReader(open(path), delimiter="\t"):
            for s in bycontig.get(r["contig"], []):
                p = int(s["sec_codon_start"])
                if not (int(r["tstart"]) <= p and p + 2 <= int(r["tend"])): continue
                blocks = [tuple(map(int, b.split("-"))) for b in r["blocks"].split(",")]
                if any(x <= p and p + 2 <= y for x, y in blocks):
                    c = counts[(s["site"], lab)]; c[0] += 1
                    ext = (int(r["tend"]) - p - 2) if s["strand"] == "+" else (p - int(r["tstart"]))
                    c[1] = max(c[1], ext)
    with open(a.output, "w") as o:
        o.write("site\t" + "\t".join(f"{l}_n\t{l}_3prime_extent" for l in a.labels) + "\n")
        for s in sites:
            o.write(s["site"] + "\t" + "\t".join(f"{counts[(s['site'], l)][0]}\t{counts[(s['site'], l)][1]}" for l in a.labels) + "\n")

if __name__ == "__main__":
    main()
