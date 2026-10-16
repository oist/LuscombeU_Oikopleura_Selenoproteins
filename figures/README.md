# Figures

Figure scripts read `results/tables/` and, where noted, the pipeline outputs in `work/`.

| Script | Figure | Inputs |
|---|---|---|
| `stage1_secis_scores.py` | RF00031 scores of the Selenoprofiles 3′ sequences and of random genomic windows, per genome | `work/stage1/pipeline/selenoprofiles/*.three_prime.fa`, `data/genomes/`, `data/references/RF00031.cm`; writes `results/tables/stage1_secis_scores.tsv` |
| `stage1_secis_scores.py --cm work/stage2/secis_model/appendicularia_secis.cm --baseline curation/baseline_sec_genes.tsv --seed 1000 --table results/tables/stage2_secis_scores.tsv --figure figures/stage2_secis_scores --model-name 'appendicularian SECIS model'` | the same scores with the stage 2 model, with the curated *O. dioica* Sec genes and new random windows | as above, plus the model built by `stage2/build_secis_model.py` |
