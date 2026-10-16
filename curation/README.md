# Curation

Baseline of the three *O. dioica* assemblies, written as data:

| File | Contents |
|---|---|
| `baseline_sec_genes.tsv` | 18 Sec genes per assembly: Sec codon (lowest coordinate of the TGA), strand, Sec context of the protein alignment, status (`sec`, `uncertain`), notes on assembly errors |
| `machinery_pstk.tsv` | PSTK protein alignment spans (RefSeq model truncated in OKI2018) |
| `machinery_proteins.faa` | curated *O. dioica* machinery proteins from transcripts, added to the appendicularian Selenoprofiles profiles (header: family, genome, source): PSTK from the OKI2018 Trinity transcriptome and from the OSKA2016 transcript assembly; SEPHS2, with its Sec, from the OSKA2016 transcript assembly |
| `candidates.tsv` | candidate Sec genes not counted in the baseline (thioredoxin-like protein 1) |

The Sec codons were placed by aligning curated reference proteins and checking the Sec context; transcripts and
cross-assembly comparison support them (see the stage 1 tables). Assembly errors are recorded, not corrected.

`check_baseline.py` (run by `make curation`) reports, for each curated Sec codon, which stage 1 tools found it, the
RF00031 score and the transcript support: `results/tables/baseline_vs_stage1.tsv` and `_summary.tsv`.
