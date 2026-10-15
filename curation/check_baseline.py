#!/usr/bin/env python3
"""Compare the curated baseline (curation/baseline_sec_genes.tsv) with stage 1 Sec sites.

For each curated Sec codon: which stage 1 tools found it (a site on the same contig and strand within 3 nt), the best
RF00031 score in the 3 kb downstream, and the transcript support. Writes results/tables/baseline_vs_stage1.tsv
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
    out = ["genome\tgene\tstatus\tfound_by\trf00031_bits\ttranscript_support\toffset"]
    summ = defaultdict(lambda: defaultdict(int))
    genomes_done = {k[0] for k in sites}
    for b in csv.DictReader(open(a.baseline), delimiter="\t"):
        g = b["genome"]
        if g not in genomes_done: continue
        p = int(b["sec_codon"])
        hits = [s for s in sites[(g, b["contig"], b["strand"])] if abs(int(s["sec_codon_start"]) - p) <= a.tolerance]
        h = max(hits, key=lambda s: float(s["rf00031_bits"])) if hits else None
        found = h["sources"].replace(",", "+") if h else "none"
        out.append(f"{g}\t{b['gene']}\t{b['status']}\t{found}\t{h['rf00031_bits'] if h else ''}\t{h.get('transcript_support', '') if h else ''}\t{int(h['sec_codon_start']) - p if h else ''}")
        s = summ[g]; s["curated"] += 1
        if h:
            s["found"] += 1
            s["selenoprofiles"] += "selenoprofiles" in h["sources"]
            s["miniprot"] += "miniprot" in h["sources"]
            s["rf00031_ge20"] += float(h["rf00031_bits"]) >= 20
            s["rf00031_ge35"] += float(h["rf00031_bits"]) >= 35
    os.makedirs(a.out_dir, exist_ok=True)
    open(f"{a.out_dir}/baseline_vs_stage1.tsv", "w").write("\n".join(out) + "\n")
    cols = ["curated", "found", "selenoprofiles", "miniprot", "rf00031_ge20", "rf00031_ge35"]
    lines = ["genome\t" + "\t".join(cols)] + [g + "\t" + "\t".join(str(summ[g][c]) for c in cols) for g in sorted(summ)]
    open(f"{a.out_dir}/baseline_vs_stage1_summary.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
