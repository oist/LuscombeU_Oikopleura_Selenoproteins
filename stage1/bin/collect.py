#!/usr/bin/env python3
"""Collect stage 1 outputs of all genomes into small tables under results/tables/:
stage1_sec_sites.tsv     one row per Sec site (Selenoprofiles and/or miniprot) with the best RF00031 SECIS score
                         in the 3 kb downstream
stage1_summary.tsv       per genome: sites by source, sites with an RF00031 SECIS at 20 bits and at the Rfam
                         gathering threshold (35 bits), tRNA-Sec loci from Secmarker
stage1_annotation.tsv    annotation check of the genomes that have an annotation (Tico et al. 2026 classes)"""
import argparse, csv, glob, os, re

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True); ap.add_argument("--genomes", nargs="+", required=True)
    ap.add_argument("--out-dir", required=True); ap.add_argument("--ga", type=float, default=35.0)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sites_out = open(f"{a.out_dir}/stage1_sec_sites.tsv", "w")
    sites_out.write("genome\tsite\tcontig\tsec_codon_start\tstrand\tfamilies\tsources\trf00031_bits\n")
    summ = ["genome\tsites\tselenoprofiles\tminiprot\tboth\trf00031_ge20\trf00031_ge_ga\ttrnasec_loci\ttrnasec_best_bits"]
    annot_rows, annot_header = [], None
    for g in a.genomes:
        d = f"{a.work}/{g}"
        best = {}
        for l in open(f"{d}/{g}.secis_rf00031.tbl"):
            if l.startswith("#"): continue
            f = l.split(); best[f[0]] = max(best.get(f[0], 0.0), float(f[14]))
        rows = list(csv.DictReader(open(f"{d}/{g}.sec_sites.tsv"), delimiter="\t"))
        for r in rows:
            sites_out.write(f"{g}\t{r['site']}\t{r['contig']}\t{r['sec_codon_start']}\t{r['strand']}\t{r['families']}\t{r['sources']}\t{best.get(r['site'], 0.0):.1f}\n")
        src = [set(r["sources"].split(",")) for r in rows]
        trna = []
        for gff in glob.glob(f"{d}/secmarker/**/*.gff", recursive=True) + glob.glob(f"{d}/secmarker/*.gff"):
            for l in open(gff):
                f = l.rstrip("\n").split("\t")
                if len(f) > 8:
                    m = re.search(r"infernal_score:(\S+)", f[8]); trna.append(float(m.group(1)) if m else 0.0)
        summ.append("\t".join(map(str, [g, len(rows), sum("selenoprofiles" in s for s in src), sum("miniprot" in s for s in src),
                                        sum(len(s) == 2 for s in src), sum(best.get(r["site"], 0) >= 20 for r in rows),
                                        sum(best.get(r["site"], 0) >= a.ga for r in rows), len(set(trna)) if trna else 0,
                                        max(trna) if trna else ""])))
        for p in glob.glob(f"{d}/{g}.sec_annotation.tsv"):
            with open(p) as fh:
                r = list(csv.reader(fh, delimiter="\t"))
            if r:
                annot_header = annot_header or ["genome"] + r[0]
                annot_rows += [[g] + x for x in r[1:]]
    sites_out.close()
    open(f"{a.out_dir}/stage1_summary.tsv", "w").write("\n".join(summ) + "\n")
    if annot_header:
        with open(f"{a.out_dir}/stage1_annotation.tsv", "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(annot_header); w.writerows(annot_rows)

if __name__ == "__main__":
    main()
