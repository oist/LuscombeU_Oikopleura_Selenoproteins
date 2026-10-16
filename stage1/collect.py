#!/usr/bin/env python3
"""Collect the stage 1 pipeline outputs into small tables under results/tables/:

stage1_sec_sites.tsv      one row per Sec codon (Selenoprofiles and/or miniprot), with the SECIS score of the
                          Selenoprofiles prediction (RF00031) and, where the samplesheet gave transcripts, their support
stage1_genes_summary.tsv  per genome: Selenoprofiles predictions by label, Sec sites by tool, SECIS support
stage1_trnasec.tsv        per genome: tRNA-Sec loci by tool and support
stage1_annotation.tsv     per Sec codon of the genomes with an annotation: its class in that annotation
"""
import argparse, csv, glob, os


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concat(pattern, suffix, out, first_column="genome"):
    """Concatenate per-genome tables (genome name taken from the file name) into one table."""
    header, lines = None, []
    for path in sorted(glob.glob(pattern)):
        genome = os.path.basename(path)[: -len(suffix)]
        for r in rows(path):
            header = header or [first_column] + [k for k in r if k != "sample"]
            lines.append([genome] + [r[k] for k in header[1:]])
    if header:
        with open(out, "w") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(header); w.writerows(lines)
    return len(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="work/stage1/pipeline", help="pipeline output directory")
    ap.add_argument("--out-dir", default="results/tables")
    ap.add_argument("--prefix", default="stage1")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    o = lambda name: f"{a.out_dir}/{a.prefix}_{name}.tsv"

    header = ["genome", "contig", "sec_codon_start", "sec_codon_end", "strand", "family", "selenoprofiles_label",
              "sources", "secis_bits", "secis_model", "transcripts", "max_3prime_extent"]
    n = 0
    with open(o("sec_sites"), "w") as out:
        out.write("\t".join(header) + "\n")
        for path in sorted(glob.glob(f"{a.outdir}/compare/*.sec_sites.tsv")):
            genome = os.path.basename(path)[: -len(".sec_sites.tsv")]
            tx = {}
            tpath = f"{a.outdir}/transcripts/{genome}.sec_transcripts.tsv"
            if os.path.exists(tpath):
                tx = {(r["contig"], r["start"], r["strand"]): r for r in rows(tpath)}
            for r in rows(path):
                t = tx.get((r["contig"], r["start"], r["strand"]))
                out.write("\t".join([genome, r["contig"], r["start"], r["end"], r["strand"], r["family"],
                                     r["selenoprofiles_label"], r["support"], r["secis_score"], r["secis_model"],
                                     t["transcripts"] if t else "", t["max_3prime_extent"] if t else ""]) + "\n")
                n += 1
    g = concat(f"{a.outdir}/compare/*.genes_summary.tsv", ".genes_summary.tsv", o("genes_summary"))
    t = concat(f"{a.outdir}/compare/*.trnasec_summary.tsv", ".trnasec_summary.tsv", o("trnasec"))
    s = concat(f"{a.outdir}/annotation/*.sec_annotation.tsv", ".sec_annotation.tsv", o("annotation"))
    print(f"{n} Sec sites; {g} genome summaries; {t} tRNA-Sec summaries; {s} annotation rows")


if __name__ == "__main__":
    main()
