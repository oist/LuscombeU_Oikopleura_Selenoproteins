#!/usr/bin/env python3
"""Sequences around SELENOI for a wider search of its SECIS element (step 9, follow-up of the absence tests,
2026-10-31): genome regions on the gene strand and transcripts, written to a small FASTA file so that the search
(stage2/secis_stemloop_screen.py) can run without the genomes.

Regions (sense strand of the gene; coordinates 1-based, inclusive, on the forward strand of the contig):
  - O. dioica (OKI2018, OSKA2016, OdB3): from --before nt before to --after nt after the candidate codon placed by
    the transcripts (results/tables/stage2_candidate_transcripts_positions.tsv);
  - every genome of the tblastn table of candidate_genes.py (results/tables/stage2_candidates_genomes.tsv, SELENOI
    rows: appendicularians, other tunicates, amphioxus): the best locus, extended by --before nt upstream and --after
    nt downstream on the gene strand;
  - the stage 1 SELENOI selenocysteine sites outside the appendicularians (results/tables/stage1_sec_sites.tsv:
    vertebrates, C. intestinalis), from --before nt before to --after nt after the codon;
  - the SELENOI transcripts of the three O. dioica transcript sets (results/tables/stage2_candidate_transcripts.tsv,
    all isoforms), in the orientation of their open reading frame.
Names: <genome>|<kind>|<contig>:<start>-<end>|<strand>|<note>; the table (--table) gives, for each region, the
offset of the codon in the sequence (0-based), when there is one.
"""
import argparse, csv, gzip, os

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s): return s.translate(COMP)[::-1]


def read_fasta(path, keep):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in keep: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            elif name in keep: buf.append(l.strip())
    if name in keep: yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--positions", default="results/tables/stage2_candidate_transcripts_positions.tsv")
    ap.add_argument("--transcripts", default="results/tables/stage2_candidate_transcripts.tsv")
    ap.add_argument("--genomes-table", default="results/tables/stage2_candidates_genomes.tsv")
    ap.add_argument("--stage1", default="results/tables/stage1_sec_sites.tsv")
    ap.add_argument("--gene", default="SELENOI")
    ap.add_argument("--before", type=int, default=5000)
    ap.add_argument("--after", type=int, default=8000)
    ap.add_argument("--out", default="results/sequences/selenoi_regions.fa")
    ap.add_argument("--table", default="results/sequences/selenoi_regions.tsv")
    a = ap.parse_args()
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    # (genome, kind, contig, strand, lo, hi, codon or None, note); lo, hi 1-based inclusive, codon 1-based start
    want = []
    for r in csv.DictReader(open(a.positions), delimiter="\t"):
        if r["gene"] != a.gene: continue
        c = int(r["codon"])
        want.append((r["genome"], "codon", r["contig"], r["strand"], c, c + 2, c,
                     f"candidate codon placed by transcripts ({r['transcripts']})"))
    for r in csv.DictReader(open(a.genomes_table), delimiter="\t"):
        if r["candidate"] != a.gene: continue
        c = int(r["position"]) if r.get("position") else None
        want.append((r["genome"], "tblastn_locus", r["contig"], r["strand"], int(r["from"]), int(r["to"]), c,
                     f"best tblastn locus, codon {r['codon'] or 'not aligned'}"))
    for r in csv.DictReader(open(a.stage1), delimiter="\t"):
        if r["family"] != a.gene or r["genome"] in ("OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery",
                                                     "Fbor"): continue
        c = int(r["sec_codon_start"])
        want.append((r["genome"], "stage1_site", r["contig"], r["strand"], c, c + 2, c,
                     f"stage 1 selenocysteine site ({r['sources']}, label {r['selenoprofiles_label']})"))
    by_genome = {}
    for w in want: by_genome.setdefault(w[0], []).append(w)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    rows = []
    with open(a.out, "w") as o:
        for g, ws in by_genome.items():
            if g not in sheet:
                print(f"{g}: not in the samplesheet"); continue
            contigs = dict(read_fasta(sheet[g]["fasta"], {w[2] for w in ws}))
            seen = set()
            for _, kind, contig, st, lo, hi, codon, note in ws:
                seq = contigs.get(contig)
                if seq is None:
                    print(f"{g} {contig}: contig not found"); continue
                s0 = max(1, (lo - a.before) if st == "+" else (lo - a.after))
                e0 = min(len(seq), (hi + a.after) if st == "+" else (hi + a.before))
                if (contig, s0, e0, st) in seen: continue
                seen.add((contig, s0, e0, st))
                sub = seq[s0 - 1:e0].upper()
                sub = sub if st == "+" else rc(sub)
                off = "" if codon is None else ((codon - s0) if st == "+" else (e0 - codon - 2))
                name = f"{g}|{kind}|{contig}:{s0}-{e0}|{st}"
                o.write(f">{name}|{note.replace(' ', '_')}\n{sub}\n")
                rows.append((name, g, kind, contig, s0, e0, st, len(sub), off,
                             sub[off:off + 3] if off != "" else "", note))
        txs = [r for r in csv.DictReader(open(a.transcripts), delimiter="\t") if r["gene"] == a.gene]
        for g in sorted({r["genome"] for r in txs}):
            path = sheet[g].get("transcripts", "")
            if not path or not os.path.exists(path): continue
            strands = {r["transcript"]: r["transcript_strand"] for r in txs if r["genome"] == g}
            for k, s in read_fasta(path, set(strands)):
                s = s.upper() if strands[k] == "+" else rc(s.upper())
                name = f"{g}|transcript|{k}|{strands[k]}"
                o.write(f">{name}|transcript_set\n{s}\n")
                rows.append((name, g, "transcript", k, 1, len(s), strands[k], len(s), "", "", "transcript"))
    with open(a.table, "w") as o:
        o.write("name\tgenome\tkind\tcontig\tstart\tend\tstrand\tlength\tcodon_offset\tcodon\tnote\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} sequences -> {a.out}")


if __name__ == "__main__":
    main()
