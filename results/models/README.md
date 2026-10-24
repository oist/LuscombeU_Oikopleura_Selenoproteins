# Models

| File | Content |
|---|---|
| `appendicularia_secis.cm` | first appendicularian SECIS covariance model (Infernal 1.1.5, calibrated), built by `make secis-model` from the stage 1 predictions; replaced by `oikopleura_secis.cm` for the SECIS of genes, still used by the stage 1 SECIS figures and `make secis-hits` |
| `appendicularia_secis.sto` | its seed alignment (73 sequences) |
| `machinery_profiles/` | appendicularian Selenoprofiles profiles of the Sec machinery (SEPHS2, pstk, SecS, eEFsec, SBP2, secp43), with the standard names: aligned sequences (`.fa`), configuration (`.fa.config`) and the list of added sequences (`members.tsv`); built by `make machinery-profiles` |
| `appendicularia_machinery_profiles.tar.gz` | the same profiles with their precomputed `.profile_data`, as passed to the pipeline (`--selenoprofiles_overrides`) and deposited on Zenodo |
| `appendicularia_secis.seeds.fa` | the seed sequences; their names give genome, prediction and position in its 3′ sequence |
| `oikopleura_secis.cm`, `.sto`, `.seeds.fa` | SECIS model in use (step 3; `make secis-v2`), trained on the windows after the stop codons of the gold-standard genes of the three *O. dioica* assemblies and of the stage 1 predictions of three other *Oikopleura* species (76 seeds, 66 consensus positions); its seed alignment and seeds (window names give genome and gene or prediction, and the position in the window after the stop codon). Threshold for the 3′ windows of genes: 20 bits |

The model files record their build date and command line, so they change at each build; to check that a rebuild
reproduces them, compare the seed sequences (`*.seeds.fa`) and the rounds tables
(`results/tables/stage2_secis_model_rounds.tsv`, `results/tables/stage2_secis_v2_rounds.tsv`).
