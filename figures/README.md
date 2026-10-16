# Figures

Figure scripts read `results/tables/` and, where noted, the pipeline outputs in `work/`. Run them with `make` (see
the `secis-*` targets of the Makefile); they need Infernal and matplotlib, available in the image.

| Script | Make target | Figure |
|---|---|---|
| `stage1_secis_scores.py` | `secis-stage1` | RF00031 scores of the Selenoprofiles 3′ sequences and of random genomic windows, per genome (`stage1_secis_scores.*`; table `results/tables/stage1_secis_scores.tsv`) |
| `stage1_secis_scores.py` with the stage 2 model | `secis-stage2` | the same scores with the appendicularian SECIS model, with the curated *O. dioica* Sec genes and new random windows (`stage2_secis_scores.*`) |

`--max` turns off the `cmsearch` filters, so that every sequence gets its best score instead of "no hit" (slower);
`--replot` redraws a figure from its table.
