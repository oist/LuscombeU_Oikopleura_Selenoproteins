#!/usr/bin/env python3
"""Genomic regions of one gene in several assemblies, written as small FASTA files for inspection outside the
pipeline (for example the splice sites of the PSTK introns in assemblies of other individuals and other species).

Regions, by source (1-based, inclusive coordinates; --flank nt added on each side, clipped to the contig):
  --gold GENOME ...           span of the coding rows of the gene in <gold-dir>/<GENOME>.gff3 (gold-standard model).
  --lift TARGET=SOURCE ...    the gold-standard region of SOURCE (span and --flank) aligned to TARGET with minimap2
                              (-x asm20); every alignment that covers at least --min-cover of the region is written
                              (best first: the locus, then possible copies such as uncollapsed haplotypes), its target
                              span extended by --flank.
  --annotation GENOME=PROTEIN the coding rows of a protein in the annotation of GENOME (samplesheet column
                              "annotation"), and the GFF3 rows of its gene, mRNA, exons and CDS (<out>/<GENOME>.gff3).
  --contig GENOME=CONTIG[:START-END] ...  a region given by hand (the whole contig without a range).
Genome FASTA and annotation paths: --samplesheet files (sample, fasta, annotation columns).
Outputs in --out: <GENOME>.fa (one record per region; header <GENOME>|<contig>:<start>-<end>|<strand>|<source>, the
strand being that of the gene when known: the model strand, or for --lift the source strand carried through the
alignment) and regions.tsv (genome, contig, start, end, strand, source, detail).
"""
import argparse, csv, gzip, os, re, subprocess, tempfile


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_contigs(path, names):
    """{name: sequence} of the named contigs of a FASTA file (read once, keeping only those)."""
    out, name, buf = {}, None, []
    with opener(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in names: out[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            elif name in names: buf.append(l.strip())
    if name in names: out[name] = "".join(buf)
    return out


def gold_span(path, gene):
    """(contig, strand, lo, hi) of the coding rows of a gene in a gold-standard GFF3."""
    rows = []
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "CDS" and re.search(rf"Parent=rna-[^;]*_{re.escape(gene)}(;|$)", f[8]):
            rows.append((f[0], f[6], int(f[3]), int(f[4])))
    if not rows: return None
    return rows[0][0], rows[0][1], min(r[2] for r in rows), max(r[3] for r in rows)


def annotation_rows(path, protein):
    """GFF3 rows of the gene, mRNA, exons and CDS of a protein (RefSeq-style attributes)."""
    lines = [l.rstrip("\n").split("\t") for l in opener(path) if not l.startswith("#")]
    lines = [f for f in lines if len(f) == 9]
    attr = lambda f, k: (re.search(rf"(?:^|;){k}=([^;]+)", f[8]) or [None, None])[1]
    cds = [f for f in lines if f[2] == "CDS" and attr(f, "protein_id") == protein]
    if not cds: return []
    mrna = attr(cds[0], "Parent")
    m = next((f for f in lines if attr(f, "ID") == mrna), None)
    gene = attr(m, "Parent") if m else None
    keep = [f for f in lines if attr(f, "ID") in (gene, mrna) or attr(f, "Parent") == mrna]
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--gold", nargs="*", default=[])
    ap.add_argument("--lift", nargs="*", default=[])
    ap.add_argument("--annotation", nargs="*", default=[])
    ap.add_argument("--contig", nargs="*", default=[])
    ap.add_argument("--samplesheet", nargs="+", default=["stage1/samplesheet.csv", "stage2/run/samplesheet_final.csv"])
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--min-cover", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sheet = {}
    for s in a.samplesheet:
        for r in csv.DictReader(open(s)):
            sheet.setdefault(r["sample"], r)
    os.makedirs(a.out, exist_ok=True)
    regions = []          # genome, contig, start, end, strand, source, detail
    for g in a.gold:
        sp = gold_span(f"{a.gold_dir}/{g}.gff3", a.gene)
        if sp: regions.append([g, sp[0], sp[2] - a.flank, sp[3] + a.flank, sp[1], "gold", "coding span ± flank"])
        else: print(f"{g}: no gold-standard model of {a.gene}")
    for pair in a.lift:
        tgt, src = pair.split("=")
        sp = gold_span(f"{a.gold_dir}/{src}.gff3", a.gene)
        if not sp: print(f"{tgt}: no gold-standard model of {a.gene} in {src}"); continue
        contig, strand, lo, hi = sp
        lo, seq = max(1, lo - a.flank), read_contigs(sheet[src]["fasta"], {contig})[contig]
        hi = min(len(seq), hi + a.flank)
        with tempfile.TemporaryDirectory() as tmp:
            q = f"{tmp}/q.fa"
            with open(q, "w") as o: o.write(f">{src}_{a.gene}\n{seq[lo - 1:hi]}\n")
            paf = subprocess.run(["minimap2", "-x", "asm20", "-t", str(a.threads), sheet[tgt]["fasta"], q],
                                 capture_output=True, text=True, check=True).stdout
        recs = []
        for l in paf.splitlines():
            f = l.split("\t")
            ql, qs, qe, ori, t, ts, te, match, alen = int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
            recs.append((match, ql, qs, qe, ori, t, ts, te, alen))
        # alignments of one locus can be split; merge those on the same target sequence and orientation within 5 kb
        recs.sort(key=lambda r: (r[5], r[4], r[6]))
        merged = []
        for r in recs:
            m = merged[-1] if merged else None
            if m and m[5] == r[5] and m[4] == r[4] and r[6] - m[7] < 5000:
                m[0] += r[0]; m[2] = min(m[2], r[2]); m[3] = max(m[3], r[3]); m[7] = max(m[7], r[7]); m[8] += r[8]
            else: merged.append(list(r))
        merged = [m for m in merged if (m[3] - m[2]) >= a.min_cover * m[1]]
        merged.sort(key=lambda m: -m[0])
        if not merged: print(f"{tgt}: region of {src} {a.gene} not found"); continue
        for k, (match, ql, qs, qe, ori, t, ts, te, alen) in enumerate(merged):
            gstrand = strand if ori == "+" else ("-" if strand == "+" else "+")
            regions.append([tgt, t, ts + 1 - a.flank, te + a.flank, gstrand, f"lift:{src}" + ("" if k == 0 else f":copy{k}"),
                            f"query {qs + 1}-{qe} of {ql} nt aligned, {match} matching bases in {alen} columns"])
    for pair in a.annotation:
        g, prot = pair.split("=")
        rows = annotation_rows(sheet[g]["annotation"], prot)
        if not rows: print(f"{g}: no coding rows for {prot}"); continue
        cds = [f for f in rows if f[2] == "CDS"]
        regions.append([g, cds[0][0], min(int(f[3]) for f in cds) - a.flank, max(int(f[4]) for f in cds) + a.flank,
                        cds[0][6], f"annotation:{prot}", f"{len(cds)} coding rows"])
        with open(f"{a.out}/{g}.gff3", "w") as o:
            o.write("##gff-version 3\n" + "".join("\t".join(f) + "\n" for f in rows))
    for spec in a.contig:
        g, rest = spec.split("=")
        m = re.match(r"^([^:]+)(?::(\d+)-(\d+))?$", rest)
        lo, hi = (int(m.group(2)) - a.flank, int(m.group(3)) + a.flank) if m.group(2) else (1, 10 ** 12)
        regions.append([g, m.group(1), lo, hi, ".", "contig", "region given by hand ± flank" if m.group(2) else "whole contig"])
    by_genome = {}
    for r in regions: by_genome.setdefault(r[0], []).append(r)
    with open(f"{a.out}/regions.tsv", "w") as t:
        t.write("genome\tcontig\tstart\tend\tstrand\tsource\tdetail\n")
        for g, rs in by_genome.items():
            seqs = read_contigs(sheet[g]["fasta"], {r[1] for r in rs})
            with open(f"{a.out}/{g}.fa", "w") as o:
                for r in rs:
                    s = seqs[r[1]]
                    r[2], r[3] = max(1, r[2]), min(len(s), r[3])
                    o.write(f">{g}|{r[1]}:{r[2]}-{r[3]}|{r[4]}|{r[5]}\n")
                    sub = s[r[2] - 1:r[3]]
                    o.write("".join(sub[i:i + 80] + "\n" for i in range(0, len(sub), 80)))
                    t.write("\t".join(map(str, r)) + "\n")
    print(f"{len(regions)} regions written to {a.out}")


if __name__ == "__main__":
    main()
