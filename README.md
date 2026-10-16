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

## Before running

All the software is listed, with pinned versions, in `environment.yml`. There are three ways to provide it.

1. **Container image (recommended; no installation).** Build it once with Wave from `environment.yml` (commands at the
   top of the file), then fill in the image addresses in `scripts/run_in_container.sh` and
   `conf/pipeline_single_container.config`. Run each step through the wrapper, which uses Singularity or Apptainer:
   `scripts/run_in_container.sh make curation`. Directories outside the repository that `data/` links point to must
   be bound (`export SINGULARITY_BIND=<directory>`).
2. **Environment modules.** On a cluster that provides the tools as modules, load them before running the scripts.
   Module names vary between sites; for example, `scripts/fetch.sh` needs the NCBI datasets command-line tool
   (at OIST: `ml bioinfo-ugrp-modules Other/ncbi-datasets-cli`), curl, unzip, tar and gzip.
3. **conda-forge installer** (Miniforge), where allowed: `conda env create -f environment.yml`, then
   `conda activate oikopleura-selenoproteins`.

Also: `git config core.hooksPath .githooks` refuses commits of files larger than 1 MB.

## Run

From the repository root (prefix each `make` with `scripts/run_in_container.sh` when using the image):

```bash
make fetch          # needs network; some sources are still placeholders (see Data)
PROFILE=<profile> stage1/run.sh   # the pipeline (Nextflow)
make stage1
make curation
make stage2
make figures
```

### On a Slurm cluster

Scripts that run as batch jobs carry generic `#SBATCH` headers (job name prefixed with `selenoproteins-`, CPUs,
memory, time limit, log file named `<job name>-<job id>.log` in the submission directory). Submit them from the repository root and give site-specific
options (partition, account) on the command line, which also overrides the header values:

```bash
sbatch [--partition=<data transfer partition>] scripts/fetch.sh
PROFILE=<profile> sbatch [options] stage1/run.sh                # Nextflow driver for the stage 1 pipeline run
sbatch [options] scripts/run_in_container.sh make stage1       # collects its outputs
```

### The pipeline in the same image

`conf/pipeline_single_container.config` makes every process of `oist/LuscombeU_Selenoproteins_Pipeline` use the image
built from `environment.yml` instead of the per-module images:

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r <release> -profile <profile> \
    -c conf/pipeline_single_container.config --input <samplesheet> --outdir <directory>
```

## Data

| Set | Source |
|---|---|
| *O. dioica* OKI2018_I68 genome and RefSeq annotation | GCF_907165135.1 |
| *O. dioica* OSKA2016 genome and annotation | ANISEED archives `Oidioi_genome_fasta_gff3.zip` and `Oidioi_transcripts_fasta_gff3.zip` (`config/archives.tsv`); the ANISEED download links are currently broken: place copies in `data/downloads/` (checked by md5); a Zenodo copy is planned |
| *O. dioica* OdB3 genome and annotation | GCA_000209535.1 |
| Other appendicularians (Naville et al. 2019) | *Oikopleura longicauda* GCA_004367895.1, *Oikopleura albicans* GCA_004367875.1, *Oikopleura vanhoeffeni* GCA_004367855.1, *Bathochordaeus stygius* GCA_004367955.1, *Mesochordaeus erythrocephalus* GCA_004367975.1, *Fritillaria borealis* GCA_004368075.1 |
| Other tunicates (SECIS non-target controls) | *Ciona intestinalis* GCF_018327825.1, *Ciona savignyi* GCA_054131345.1, *Styela clava* GCF_964204865.1, *Thalia democratica* GCF_965202585.1 |
| OSKA2016 transcripts (Wang et al. 2015) | GEO GSE64421 supplementary file `GSE64421_transcripts.trinity.fasta.gz` |
| OKI2018 transcriptome (Bliznina et al. 2021) | not yet public; to be added to Zenodo record 4604144 |
| OdB3 cDNA (Denoeud et al. 2010) | Genoscope `cDNA.fa.bz2` (OikoBase, which also distributed them, is offline) |
| OdB3 operon annotation (Denoeud et al. 2010) | Genoscope `operons_GGB.gff.bz2` (`config/extra_files.tsv`) |
| Rfam RF00031, Selenoprofiles profiles v1.3.0 | `config/references.tsv` |

**Assembly name.** The OKI2018 assembly is named `OKI2018_I68` at NCBI (GCF_907165135.1) because of a typo: it is
the `OKI2018_I69` assembly of Bliznina et al. 2021, where I69 identifies the sequenced *O. dioica* individual. This
repository uses `OKI2018_I68`, the NCBI name, throughout.

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
