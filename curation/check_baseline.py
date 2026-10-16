#!/usr/bin/env python3
"""Compare the curated baseline (curation/baseline_sec_genes.tsv) with stage 1 Sec sites.

For each curated Sec codon: which stage 1 tools found it (a site on the same contig and strand within 3 nt), the
SECIS score of the Selenoprofiles prediction (at the threshold of the run), and the number of supporting transcripts. Writes results/tables/baseline_vs_stage1.tsv
(one row per curated gene and genome) and results/tables/baseline_vs_stage1_summary.tsv (per genome)."""
import argparse, csv, os
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--baseline", default=f"{here}/baseline_sec_genes.tsv")
    ap.add_argument("--stage1", default=f"{here}/../results/tables/stage1_sec_sites.tsv")
    ap.add_argument("--out-dir", default=f"{here}/../results/tables")
    ap.add_argument("--tolerance", type=int, default=3)
    a = ap.parse_args()
    sites = defaultdict(list)
    for r in csv.DictReader(open(a.stage1), delimiter="\t"):
        sites[(r["genome"], r["contig"], r["strand"])].append(r)
    out = ["genome\tgene\tstatus\tfound_by\tsecis_bits\ttranscripts\toffset"]
    summ = defaultdict(lambda: defaultdict(int))
    genomes_done = {k[0] for k in sites}
    for b in csv.DictReader(open(a.baseline), delimiter="\t"):
        g = b["genome"]
        if g not in genomes_done: continue
        p = int(b["sec_codon"])
        hits = [s for s in sites[(g, b["contig"], b["strand"])] if abs(int(s["sec_codon_start"]) - p) <= a.tolerance]
        bits = lambda s: float(s["secis_bits"]) if s["secis_bits"] not in ("", ".") else 0.0
        h = max(hits, key=bits) if hits else None
        found = h["sources"] if h else "none"
        out.append(f"{g}\t{b['gene']}\t{b['status']}\t{found}\t{h['secis_bits'] if h else ''}\t{h['transcripts'] if h else ''}\t{int(h['sec_codon_start']) - p if h else ''}")
        s = summ[g]; s["curated"] += 1
        if h:
            s["found"] += 1
            s["selenoprofiles"] += h["sources"] in ("selenoprofiles", "both")
            s["miniprot"] += h["sources"] in ("miniprot", "both")
            s["secis"] += bits(h) > 0
            s["transcripts"] += h["transcripts"] not in ("", "0")
    os.makedirs(a.out_dir, exist_ok=True)
    open(f"{a.out_dir}/baseline_vs_stage1.tsv", "w").write("\n".join(out) + "\n")
    cols = ["curated", "found", "selenoprofiles", "miniprot", "secis", "transcripts"]
    lines = ["genome\t" + "\t".join(cols)] + [g + "\t" + "\t".join(str(summ[g][c]) for c in cols) for g in sorted(summ)]
    open(f"{a.out_dir}/baseline_vs_stage1_summary.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
