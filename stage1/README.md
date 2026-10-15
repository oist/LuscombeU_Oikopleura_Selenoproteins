# Stage 1 — published tools and parameters

Run with `make stage1` from the top directory (or `make -C stage1 GENOMES="OKI OSKA"`). Per genome
(`run_genome.sh`, outputs in `work/stage1/<id>/`, each step skipped when its output exists):

| Step | Tool | Output |
|---|---|---|
| 1 | Selenoprofiles 4.5.8, profiles `eukarya,machinery` | `<id>.selenoprofiles.gtf` |
| 2 | miniprot 0.11 (`-I`), queries = Sec-containing profile sequences (`bin/profile_queries.py`), Sec codons with `bin/sec_check.py` | `<id>.sec_check.tsv` |
| 3 | Sec sites merged (`bin/secis_windows.py`), Rfam RF00031 with Infernal 1.1.5 on the 3 kb downstream of each Sec codon | `<id>.sec_sites.tsv`, `<id>.secis_rf00031.tbl` |
| 4 | Secmarker 1.0.0 | `secmarker/trnasec.gff` |
| 5 | Annotation check, classes of Ticó et al. 2026 (`bin/sec_annotation_check.py`), when an annotation exists | `<id>.sec_annotation.tsv` |
| 6 | Transcript evidence for genomes listed in `config/transcripts.tsv`: spliced alignment (minimap2 through mappy, `bin/map_transcripts.py`), Sec codon inside an exon (`bin/transcript_support.py`) | `<tid>.alignments.tsv`, `<id>.transcript_support.tsv` |

`bin/collect.py` writes `results/tables/stage1_sec_sites.tsv`, `stage1_summary.tsv` and `stage1_annotation.tsv`.

`SELENOPROFILES_CACHE=<dir>` reuses `<id>.selenoprofiles.gtf` files from an earlier run of the same version and
profiles instead of re-running Selenoprofiles.

Scripts `profile_queries.py`, `sec_check.py` and `sec_annotation_check.py` come from the selenoprotein annotation
pipeline (`selenoproteins`, branch `public-data`).

