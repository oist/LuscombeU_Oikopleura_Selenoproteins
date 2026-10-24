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
| `stage2_machinery_domains.tsv`, `stage2_machinery_domain_conflicts.tsv` | `make machinery-domains` | domain analysis of the machinery trees: supports and AU tests per domain; strongly supported conflicts between the two domain trees |
| `stage2_machinery_domains_pmsf.tsv`, `stage2_machinery_domain_conflicts_pmsf.tsv` | `make machinery-domains-pmsf` | the same under LG+C60+F+G (PMSF) |
| `stage2_machinery_table.{md,tsv}` | `make machinery-table` | Table MACHINERY of the manuscript: census with orthology evidence |
| `stage2_selenoprotein_rbh.tsv` | `make selenoprotein-rbh` | reciprocal best hits of the stage 1 selenoprotein predictions against the human proteome |
| `stage2_selenoprotein_lba.tsv` | `make selenoprotein-pmsf` | supports of the largest appendicularian clade in the selenoprotein trees, maximum likelihood and PMSF |
| `stage2_selenoprotein_orthology.{tsv,md}` | `make selenoprotein-orthology` | first-pass orthology of each appendicularian selenoprotein clade |
| `stage2_secis_hits.tsv`, `stage2_secis_hits_summary.tsv`, `stage2_secis_hits_families.tsv` | `make secis-hits` | step 1: validated selenocysteine hits with their RF00031 and current-model SECIS scores; per genome and per family |
| `stage2_secis_v2_windows.tsv`, `stage2_secis_v2_rounds.tsv`, `stage2_secis_v2_summary.tsv`, `stage2_secis_v2_heldout.tsv` | `make secis-v2` | step 3: positive and cysteine windows; bootstrap rounds per fold; per fold and model, held-out windows detected and negatives above threshold, with pooled rows; held-out scores of each positive window |
| `gold_standard.tsv` | `make gold-standard` | step 2: gold-standard gene model of each curated *O. dioica* selenoprotein gene: evidence (transcript, annotated mRNA or prediction), protein source (assembly or transcript, with identities to the other assemblies), completeness, introns and non-canonical splice sites, UTRs, differences between transcripts and assembly, public annotation status and the change it needs |
