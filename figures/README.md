# Figures

Figure scripts read `results/tables/` and, where noted, the pipeline outputs in `work/`.

| Script | Figure | Inputs |
|---|---|---|
| `stage1_secis_scores.py` | RF00031 scores of the Selenoprofiles 3′ sequences and of random genomic windows, per genome | `work/stage1/pipeline/selenoprofiles/*.three_prime.fa`, `data/genomes/`, `data/references/RF00031.cm`; writes `results/tables/stage1_secis_scores.tsv` |
