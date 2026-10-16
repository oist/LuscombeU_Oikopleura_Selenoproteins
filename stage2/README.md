# Stage 2

## Appendicularian SECIS model

`build_secis_model.py` bootstraps a SECIS covariance model from the stage 1 outputs: the 3′ sequences that
Selenoprofiles wrote after its predictions labelled selenocysteine, in seven oikopleurid genomes (OKI2018, OdB3,
*O. longicauda*, *O. albicans*, *O. vanhoeffeni*, *B. stygius*, *M. erythrocephalus*). OSKA2016 and *F. borealis* are
held out. Seeds of round 0 are the RF00031 hits of at least 15 bits; each round aligns the seeds to the current model,
rebuilds and calibrates it, and keeps as next seeds the training hits above the best hit in 200 random genomic
windows per training genome. The rounds stop when the seed set no longer changes.

Inputs: `work/stage1/pipeline/selenoprofiles/*.three_prime.fa` (stage 1), `data/genomes/` and
`data/references/RF00031.cm` (`make fetch`).

### Reproduce (on a Slurm cluster, in the image built from `environment.yml`)

```bash
sbatch [--partition=<partition>] --cpus-per-task=8 scripts/run_in_container.sh make secis
```

`make secis` runs three targets, which can also be run separately:

| Target | Output | Time on 1 CPU (measured) |
|---|---|---|
| `secis-stage1` | RF00031 scores, `results/tables/stage1_secis_scores.tsv`, `figures/stage1_secis_scores.*` | 4 min |
| `secis-model` | `work/stage2/secis_model/`; the final model, seed alignment, seeds and rounds copied to `results/models/` and `results/tables/stage2_secis_model_rounds.tsv` | about 50 min (5 rounds, then a full calibration of about 25 min) |
| `secis-stage2` | scores of the new model, `results/tables/stage2_secis_scores.tsv`, `figures/stage2_secis_scores.*` | 4 min |

The scripts use as many CPUs as Slurm allocates (`SLURM_CPUS_PER_TASK`). A rebuild with the same inputs should give
the same rounds as `results/tables/stage2_secis_model_rounds.tsv` and the same seed sequences as
`results/models/appendicularia_secis.seeds.fa`: the background windows come from fixed random seeds, and Infernal's
calibration uses its fixed default seed.

### Replicate with other background windows

```bash
sbatch --cpus-per-task=8 scripts/run_in_container.sh \
    make secis-model BACKGROUND_SEED=200 SECIS_OUT=work/stage2/secis_model_seed200
```

Replicates write only to their `SECIS_OUT` directory; the files in `results/` are those of the default run
(`BACKGROUND_SEED=100`). Compare their `rounds.tsv` and `final.seeds.fa` with the default ones.

Other options of `build_secis_model.py` (`--help`): training genomes, seed threshold of round 0, number of background
windows, maximum number of rounds, calibration length of the intermediate rounds.

## Appendicularian profiles of the Sec machinery

`build_machinery_profiles.py` extends each standard Selenoprofiles profile of the machinery set (SEPHS2, pstk, SecS,
eEFsec, SBP2, secp43) with appendicularian proteins: the stage 1 predictions labelled "homologue" or "selenocysteine"
that are at least half the median profile length, and the curated transcript proteins of
`curation/machinery_proteins.faa`. They are added with `mafft --add --keeplength`, so the columns and the standard
configuration remain valid, and the profile is rebuilt with `selenoprofiles build`. OdB3 contributes no sequence and
serves as the test genome. The profiles keep the standard names, so that they can replace the standard ones.

`test_machinery_profiles.py` runs Selenoprofiles with these profiles on the nine genomes and compares the
predictions, by label, with stage 1 (`results/tables/stage2_machinery_profiles_test.tsv`).

```bash
sbatch [--partition=<partition>] --cpus-per-task=4 scripts/run_in_container.sh make machinery
```

| Target | Output | Time on 1 CPU (measured) |
|---|---|---|
| `machinery-profiles` | `work/stage2/machinery_profiles/`; copied to `results/models/machinery_profiles/` and `results/models/appendicularia_machinery_profiles.tar.gz` | 2 min |
| `machinery-test` | `work/stage2/machinery_test/` (resumable per genome); `results/tables/stage2_machinery_profiles_test.tsv` | about 75 min for the nine genomes |

| `machinery-fragments` | PSTK predictions discarded by Selenoprofiles in the genomes without a PSTK prediction, classified by blastp against the PSTK found elsewhere and KTI12 (`curation/machinery_outgroups.faa`); `results/tables/stage2_machinery_fragments.tsv` | about 5 min per genome |
| `machinery-census` | census per gene and genome, `results/tables/stage2_machinery_census.tsv` and `.md` (the table of the manuscript) | seconds |

The committed test table was produced with the same profile alignments under the names `<family>_app`, before the
profiles took the standard names; a rerun of `make machinery-test` regenerates it with the current files.
