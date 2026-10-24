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
   (at OIST: `ml bioinfo-ugrp-modules Other/ncbi-datasets-cli`), curl, unzip, tar, gzip, and `formatdb` from
   blast-legacy to format the Selenoprofiles UniRef50 tag database; without `formatdb`, run it inside the image:
   `scripts/run_in_container.sh scripts/fetch.sh`.
3. **conda-forge installer** (Miniforge), where allowed: `conda env create -f environment.yml`, then
   `conda activate oikopleura-selenoproteins`.

Also: `git config core.hooksPath .githooks` refuses commits of files larger than 1 MB.

## Run

### Where Nextflow writes its intermediate files: `NXF_WORK`

The pipeline runs (`stage1/run.sh`, and later stage 2) write every intermediate file of every task to a Nextflow
work directory: copies of the genomes, alignments, search results. By default it is `work/<stage>/nextflow` inside
the repository. **On a cluster, set `NXF_WORK` to the file system meant for such files** (often a scratch or
high-performance file system, as required by the site), and keep it for the whole run, including `-resume`:

```bash
export NXF_WORK=/path/to/scratch/selenoproteins-work      # before stage1/run.sh, or in front of the command
```

Only the work directory goes there; the pipeline outputs used by the analysis stay in `work/<stage>/pipeline`.
Delete the work directory once the run is complete and its outputs collected.

### Steps

From the repository root (prefix each `make` with `scripts/run_in_container.sh` when using the image):

```bash
make fetch                    # needs network; some sources are still placeholders (see Data)
NXF_WORK=<work directory> PROFILE=<profile> stage1/run.sh   # the pipeline (Nextflow), all genomes in one run
make stage1                   # collects its outputs into results/tables and results/pipeline_info
make stage1-figures
make curation                 # checks the baseline (curation/) against stage 1
make secis                    # appendicularian SECIS model and its scores
make machinery                # appendicularian machinery profiles, their test, fragments, census, absence tests
make machinery-rbh            # reciprocal best hits against the human proteome
make machinery-phylogeny      # trees of PSTK, SEPHS and the SECIS-binding proteins
make machinery-lba            # long-branch analyses, AU tests, C60 tree figures and Table MACHINERY
make machinery-domains        # trees and AU tests of each domain of PSTK, SEPHS and the SECIS-binding proteins
make selenoprotein-rbh        # best hits of the selenoprotein predictions against the human proteome
make selenoprotein-phylogeny  # trees of the selenoprotein families with appendicularian sequences
make selenoprotein-pmsf       # the same under LG+C60+F+G (PMSF)
make selenoprotein-orthology  # first-pass orthology of each appendicularian selenoprotein clade
```

The Makefile documents each target (inputs, outputs and the cluster submissions that worked); `stage2/README.md`
describes the stage 2 steps.

### On a Slurm cluster

Scripts that run as batch jobs carry generic `#SBATCH` headers (job name prefixed with `selenoproteins-`, CPUs,
memory, time limit, log file named `<job name>-<job id>.log` in the submission directory). Submit them from the repository root and give site-specific
options (partition, account) on the command line, which also overrides the header values:

```bash
sbatch [--partition=<data transfer partition>] scripts/fetch.sh
NXF_WORK=<work directory> PROFILE=<profile> sbatch [options] stage1/run.sh   # Nextflow driver (stage 1)
sbatch [options] scripts/run_in_container.sh make stage1                       # collects its outputs
```

Batch jobs write their log in the submission directory and their results into the repository. On clusters where the
storage that holds the repository is read-only for compute nodes (at OIST, `/bucket`), the job can write neither: its
log is lost and nothing changes. Run such jobs from a clone on a file system that compute nodes can write (at OIST,
`/flash`), with the large inputs linked from the main copy, then commit and push from the clone:

```bash
cd /flash/<unit>/<user> && git clone git@github.com:oist/LuscombeU_Oikopleura_Selenoproteins.git && cd LuscombeU_Oikopleura_Selenoproteins
ln -s <main copy>/data data                                    # read-only inputs
mkdir -p work/stage1 && ln -s <main copy>/work/stage1/pipeline work/stage1/pipeline
ml singularity                                                 # at OIST; sbatch passes the environment to the job
sbatch [options] scripts/run_in_container.sh make stage1-figures
git add results figures && git commit -m "..." && git push    # then git pull in the main copy
```

The Nextflow driver (`stage1/run.sh`) writes its history and log in the repository root and publishes the pipeline
outputs under `work/`: as a batch job it needs the repository on a file system that compute nodes can write; if the
repository is on storage that is read-only for compute nodes, run the driver on a login or interactive node.

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
| Other tunicates (stage 1 comparison, SECIS non-target controls, machinery trees) | *Ciona intestinalis* GCF_018327825.1, *Ciona robusta* (KH) GCF_000224145.3, *Ciona savignyi* GCA_054131345.1, *Halocynthia roretzi* GCA_013436055.1, *Styela clava* GCF_964204865.1, *Thalia democratica* GCF_965202585.1; for the machinery trees also *Botryllus schlosseri* GCF_051294905.1, *Phallusia mammillata* GCF_965637545.1, *Pegea* sp. GCA_034620365.1, *Salpa thompsoni* GCA_047495815.1 |
| Cephalochordates (machinery trees) | *Branchiostoma floridae* GCF_000003815.2, *Branchiostoma lanceolatum* GCF_035083965.1 |
| Vertebrate controls (stage 1; human proteome for the reciprocal best hits) | human GRCh38.p14 GCF_000001405.40, mouse GRCm39 GCF_000001635.27, zebrafish GRCz11 GCF_000002035.6, with their RefSeq annotations and proteins |
| OSKA2016 transcripts (Wang et al. 2015) | GEO GSE64421 supplementary file `GSE64421_transcripts.trinity.fasta.gz` |
| OKI2018 transcriptome (Bliznina et al. 2021) | Zenodo, <https://doi.org/10.5281/zenodo.23472675> (`config/transcripts.tsv`) |
| OdB3 cDNA (Denoeud et al. 2010) | Genoscope `cDNA.fa.bz2` (OikoBase, which also distributed them, is offline) |
| OdB3 operon annotation (Denoeud et al. 2010) | Genoscope `operons_GGB.gff.bz2` (`config/extra_files.tsv`) |
| Rfam RF00031 and RF01852; Selenoprofiles data (profiles, Gene Ontology, UniRef50 tag database) from the Selenoprofiles server | `config/references.tsv` |
| SPS1/SPS2 proteins of Mariotti et al. 2015 (SEPHS tree) | translated from their Supplemental Material S8 into `curation/mariotti2015_sps.faa` |

**Assembly name.** The OKI2018 assembly is named `OKI2018_I68` at NCBI (GCF_907165135.1) because of a typo: it is
the `OKI2018_I69` assembly of Bliznina et al. 2021, where I69 identifies the sequenced *O. dioica* individual. This
repository uses `OKI2018_I68`, the NCBI name, throughout.

## Layout

| Path | Contents |
|---|---|
| `config/` | data manifests (accessions, URLs, checksums) |
| `scripts/` | download and repository utilities |
| `Makefile` | every step, from download to figures and tables |
| `stage1/`, `stage2/` | analysis scripts (stage 1 driver and collection; stage 2 models, census, orthology) |
| `curation/` | baseline decisions as tables, with evidence; reference proteins |
| `figures/` | figure scripts and figures |
| `results/tables/` | small result tables used by the figures and the manuscript |
| `results/trees/` | machinery trees: sequence sets, removals, alignments, IQ-TREE reports; `C20/`, `lba/` and `domains/` for the other models, the long-branch analyses and the domain analysis; `selenoproteins/` for the selenoprotein families |
| `results/models/` | SECIS models (`oikopleura_secis.cm` in use) and machinery profiles |
| `manuscript/` | manuscript (Markdown with Pandoc citations) and bibliography |
| `docs/` | plan and decisions |

Related tool: the selenoprotein annotation pipeline (`selenoproteins`, Nextflow) — cited, not nested here.
