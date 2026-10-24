# Models

| File | Content |
|---|---|
| `appendicularia_secis.cm` | appendicularian SECIS covariance model (Infernal 1.1.5, calibrated), built by `make secis-model` |
| `appendicularia_secis.sto` | its seed alignment (73 sequences) |
| `machinery_profiles/` | appendicularian Selenoprofiles profiles of the Sec machinery (SEPHS2, pstk, SecS, eEFsec, SBP2, secp43), with the standard names: aligned sequences (`.fa`), configuration (`.fa.config`) and the list of added sequences (`members.tsv`); built by `make machinery-profiles` |
| `appendicularia_machinery_profiles.tar.gz` | the same profiles with their precomputed `.profile_data`, as passed to the pipeline (`--selenoprofiles_overrides`) and deposited on Zenodo |
| `appendicularia_secis.seeds.fa` | the seed sequences; their names give genome, prediction and position in its 3′ sequence |
| `oikopleura_secis.cm`, `.sto`, `.seeds.fa` | SECIS model trained on the six *Oikopleura* assemblies (step 3; `make secis-v2`), its seed alignment and seeds (window names give genome and gene or prediction, and the position in the window after the stop codon) |

The model file records its build date and command line, so it changes at each build; to check that a rebuild
reproduces this one, compare the seed sequences (`appendicularia_secis.seeds.fa`) and
`results/tables/stage2_secis_model_rounds.tsv`.
