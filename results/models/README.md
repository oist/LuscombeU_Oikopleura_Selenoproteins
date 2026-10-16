# Models

| File | Content |
|---|---|
| `appendicularia_secis.cm` | appendicularian SECIS covariance model (Infernal 1.1.5, calibrated), built by `make secis-model` |
| `appendicularia_secis.sto` | its seed alignment (73 sequences) |
| `appendicularia_secis.seeds.fa` | the seed sequences; their names give genome, prediction and position in its 3′ sequence |

The model file records its build date and command line, so it changes at each build; to check that a rebuild
reproduces this one, compare the seed sequences (`appendicularia_secis.seeds.fa`) and
`results/tables/stage2_secis_model_rounds.tsv`.
