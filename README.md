# Selenoproteins of *Oikopleura dioica* and other appendicularians

Code to reproduce the census of selenoprotein genes in appendicularian genomes, from public data, in two stages:

1. **Stage 1** runs published tools with their default or published parameters.
2. **Curation** records, as reviewed tables, the decisions that turn stage 1 output into the baseline
   selenoproteome of *O. dioica* (three assemblies: OKI2018_I68, OSKA2016, OdB3).
3. **Stage 2** uses parameters justified by stage 1 and the baseline: an *Oikopleura* SECIS covariance model,
   lineage reference proteins and profile HMMs, and re-calling on all genomes.
4. **Figures** and supplementary tables are built from the stage 1 and stage 2 result tables.

Data are not stored here: `scripts/fetch.sh` downloads them (see `config/`). Large outputs go to `work/`, which is
ignored by git; small result tables used for the figures are kept in `results/tables/`.

## Run

```bash
conda env create -f environment.yml && conda activate oikopleura-selenoproteins
git config core.hooksPath .githooks   # refuses commits of files larger than 1 MB
make fetch                            # needs network; some sources are still placeholders (see below)
make stage1
make curation
make stage2
make figures
```

## Data

| Set | Source |
|---|---|
| *O. dioica* OKI2018_I68 genome and RefSeq annotation | GCF_907165135.1 |
| *O. dioica* OSKA2016 genome and annotation | to be filled in `config/genomes.tsv` |
| *O. dioica* OdB3 genome and annotation | GCA_000209535.1 |
| Other appendicularians (Naville et al. 2019) | GCA_004367895.1, GCA_004367875.1, GCA_004367855.1, GCA_004367955.1, GCA_004367975.1, GCA_004368075.1 |
| Other tunicates (SECIS non-target controls) | GCF_018327825.1, GCA_054131345.1, GCF_964204865.1, GCF_965202585.1 |
| OKI2018 transcriptome (Bliznina et al. 2021), OSKA2016 transcripts (GEO GSE64421), OdB3 ESTs (OikoBase) | to be filled in `config/transcripts.tsv` |
| Rfam RF00031, Selenoprofiles profiles v1.3.0 | `config/references.tsv` |

## Layout

| Path | Contents |
|---|---|
| `config/` | data manifests (accessions, URLs, checksums) |
| `scripts/` | download and repository utilities |
| `stage1/`, `stage2/` | analysis steps, each with its own Makefile |
| `curation/` | baseline decisions as tables, with evidence |
| `figures/` | figure scripts |
| `results/tables/` | small result tables used by the figures |
| `docs/` | plan and methods notes |

Related tool: the selenoprotein annotation pipeline (`selenoproteins`, Nextflow) — cited, not nested here.
