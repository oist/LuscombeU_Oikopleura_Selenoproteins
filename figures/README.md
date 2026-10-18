# Figures

Figure scripts read `results/tables/` and, where noted, the pipeline outputs in `work/`. Run them with `make` (see
the `secis-*` targets of the Makefile); they need Infernal and matplotlib, available in the image.

| Script | Make target | Figure |
|---|---|---|
| `stage1_secis_scores.py` | `secis-stage1` | RF00031 scores of the Selenoprofiles 3′ sequences and of random genomic windows, per genome (`stage1_secis_scores.*`; table `results/tables/stage1_secis_scores.tsv`) |
| `stage1_secis_scores.py` with the stage 2 model | `secis-stage2` | the same scores with the appendicularian SECIS model, with the curated *O. dioica* Sec genes and new random windows (`stage2_secis_scores.*`) |

`--max` turns off the `cmsearch` filters, so that every sequence gets its best score instead of "no hit" (slower);
`--replot` redraws a figure from its table.
| `stage1_overview.py` | `stage1-figures` | stage 1 per genome (18 genomes): selenocysteine predictions and their RF00031 support, fraction in correctly annotated genes, tRNA-Sec loci and machinery genes (`stage1_overview.*`) |
| `stage1_secis_scores.py --all` | `stage1-figures` | RF00031 scores of the 3′ sequences of all eighteen genomes, with random windows where the genome is available (`stage1_secis_scores_all.*`) |
| `phylogeny_collapsed.py` | after `machinery-phylogeny` | display versions of the machinery trees: nodes with ultrafast bootstrap below 50 collapsed, clades of one taxonomic group (`config/taxon_groups.tsv`) drawn as one tip, tunicate sequences kept (`phylogeny_<family>_collapsed.*`) |
