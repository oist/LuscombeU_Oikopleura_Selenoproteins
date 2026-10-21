# Result tables

Small tab-separated tables written by the stages and read by the figure scripts and the manuscript. Large outputs stay
in `work/`.

| Table | Written by | Contents |
|---|---|---|
| `stage1_*.tsv` | `make stage1` (`stage1/collect.py`) | stage 1 outputs of the pipeline: genes, selenocysteine sites, annotation classes, tRNA-Sec, machinery, RF00031 scores |
| `stage1_secis_scores*.tsv`, `stage1_rf00031_scores.tsv` | `make secis-stage1`, `make stage1-figures` | RF00031 scores of 3′ sequences and random windows |
| `baseline_vs_stage1*.tsv` | `make curation` | curated *O. dioica* baseline against stage 1 |
| `stage2_secis_model_rounds.tsv`, `stage2_secis_scores.tsv` | `make secis-model`, `make secis-stage2` | appendicularian SECIS model: bootstrap rounds and scores |
| `stage2_machinery_profiles_test.tsv` | `make machinery-test` | appendicularian machinery profiles against stage 1 |
| `stage2_machinery_fragments.tsv` | `make machinery-fragments` | PSTK fragments discarded by Selenoprofiles |
| `stage2_machinery_census.{tsv,md}` | `make machinery-census` | machinery per gene and genome |
| `stage2_absence_{pstk,secp43}.*.tsv` | `make machinery-absence` | absence tests: candidates, controls, summary |
| `stage2_machinery_rbh.tsv` (`stage1_machinery_rbh.tsv` when regenerated) | `make machinery-rbh` | reciprocal best hits against the human proteome |
| `stage2_machinery_lba.tsv` | `make machinery-lba` | supports and AU tests of the machinery trees, all families |
| `stage2_machinery_table.{md,tsv}` | `make machinery-table` | Table MACHINERY of the manuscript: census with orthology evidence |
