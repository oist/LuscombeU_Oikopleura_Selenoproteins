#!/usr/bin/env python3
"""Collect the stage 1 pipeline outputs into small tables under results/tables/:

stage1_sec_sites.tsv      one row per Sec codon (Selenoprofiles and/or miniprot), with the SECIS score of the
                          Selenoprofiles prediction (RF00031) and, where the samplesheet gave transcripts, their support
stage1_genes_summary.tsv  per genome: Selenoprofiles predictions by label, Sec sites by tool, SECIS support
stage1_trnasec.tsv        per genome: tRNA-Sec loci by tool and support
stage1_annotation.tsv     per Sec codon of the genomes with an annotation: its class in that annotation
stage1_predictions.tsv    one row per Selenoprofiles prediction (all labels): family, label, contig, strand, span,
                          number of coding exons and protein length (from selenoprofiles/<genome>.selenoprofiles.gtf
                          and .faa)
With --prefix, the same tables for another run (for example stage2_final for the stage 2 final run).
"""
import argparse, collections, csv, glob, os, re


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


def predictions(outdir, out):
    """One row per Selenoprofiles prediction, from the GTF (coding exons) and the protein FASTA."""
    n = 0
    with open(out, "w") as handle:
        handle.write("genome\tprediction\tfamily\tlabel\tcontig\tstrand\tstart\tend\tcoding_exons\tprotein_length\n")
        for gtf in sorted(glob.glob(f"{outdir}/selenoprofiles/*.selenoprofiles.gtf")):
            genome = os.path.basename(gtf)[: -len(".selenoprofiles.gtf")]
            cds = collections.defaultdict(list)
            where = {}
            for l in open(gtf):
                f = l.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS": continue
                m = re.search(r'gene_id "([^"]+)"', f[8])
                if not m: continue
                cds[m.group(1)].append((int(f[3]), int(f[4])))
                where[m.group(1)] = (f[0], f[6])
            length = {}
            faa = gtf[: -len(".gtf")] + ".faa"
            if os.path.exists(faa):
                name = None
                for l in open(faa):
                    if l.startswith(">"): name = l[1:].split()[0]; length[name] = 0
                    elif name: length[name] += len(l.strip().rstrip("*"))
            for pid in sorted(cds):
                ex = sorted(cds[pid])
                fam, label = pid.split(".")[0], pid.rsplit(".", 1)[-1]
                handle.write("\t".join(map(str, [genome, pid, fam, label, where[pid][0], where[pid][1], ex[0][0], ex[-1][1],
                                                 len(ex), length.get(pid, "")])) + "\n")
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="work/stage1/pipeline", help="pipeline output directory")
    ap.add_argument("--out-dir", default="results/tables")
    ap.add_argument("--prefix", default="stage1")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    o = lambda name: f"{a.out_dir}/{a.prefix}_{name}.tsv"

    header = ["genome", "contig", "sec_codon_start", "sec_codon_end", "strand", "family", "selenoprofiles_label",
              "sources", "secis_bits", "secis_model", "transcripts", "max_3prime_extent", "selenoprofiles_id"]
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
                                     t["transcripts"] if t else "", t["max_3prime_extent"] if t else "",
                                     r.get("selenoprofiles_id", "")]) + "\n")
                n += 1
    g = concat(f"{a.outdir}/compare/*.genes_summary.tsv", ".genes_summary.tsv", o("genes_summary"))
    t = concat(f"{a.outdir}/compare/*.trnasec_summary.tsv", ".trnasec_summary.tsv", o("trnasec"))
    s = concat(f"{a.outdir}/annotation/*.sec_annotation.tsv", ".sec_annotation.tsv", o("annotation"))
    p = predictions(a.outdir, o("predictions"))
    print(f"{n} Sec sites; {g} genome summaries; {t} tRNA-Sec summaries; {s} annotation rows; {p} predictions")


if __name__ == "__main__":
    main()
