# Stage 1 — the pipeline with its default parameters

The samplesheet lists eighteen genomes, analysed in one run: the nine appendicularian assemblies (three of
*O. dioica*, with their annotations and transcripts, and six other species), six other tunicates (*Ciona
intestinalis*, *C. robusta* and *C. savignyi*, *Halocynthia roretzi*, *Styela clava*, *Thalia democratica*) and three vertebrates as
positive controls (human, mouse, zebrafish), with their RefSeq annotations where available.

Stage 1 runs the selenoprotein pipeline (`oist/LuscombeU_Selenoproteins_Pipeline`, release 0.2.0, set in `run.sh`) on
all these genomes, with the parameters a new user would get:

- Selenoprofiles 4.5.8, profiles `eukarya,machinery`;
- miniprot alignment of the Sec-containing profile sequences, Sec codons by Sec-context matching;
- Rfam RF00031 (SECIS) on the 3′ sequence of each Selenoprofiles prediction, at the Rfam gathering threshold;
- Secmarker and cmsearch (Rfam RF01852) for tRNA-Sec;
- for the genomes with an annotation, the annotation check (classes of Ticó et al. 2026) and the transcript
  support of each Sec codon for the *O. dioica* assemblies (OKI2018 transcriptome, OSKA2016 transcripts, OdB3 cDNA).

| File | Contents |
|---|---|
| `samplesheet.csv` | genomes, annotations and transcripts, as paths under `data/` (from `scripts/fetch.sh`) |
| `params.yaml` | pipeline parameters: local copies of the Selenoprofiles profiles and of the two Rfam models |
| `run.sh` | Nextflow driver (SBATCH headers); `PROFILE=<profile> stage1/run.sh` or `PROFILE=<profile> sbatch stage1/run.sh` |
| `collect.py` | pipeline outputs → `results/tables/stage1_*.tsv` (run by `make stage1`) |

Before the real run, run the pipeline's own test with the same release and profile; it exercises every process on
small public genomes in a few minutes and catches problems with container images or the site configuration:

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r 0.2.0 -profile test,<profile> --outdir <test directory>
```

Run from the repository root:

```bash
export NXF_WORK=/path/to/scratch/selenoproteins-stage1    # Nextflow work directory: see below
PROFILE=singularity stage1/run.sh            # or an institutional profile; resumes the previous run by default
make stage1                                  # collect the tables
make curation                                # compare with the baseline
```

**Work directory.** Nextflow keeps the intermediate files of every task in `NXF_WORK`, by default
`work/stage1/nextflow` inside the repository. On a cluster, set it to the file system meant for such files and keep
the same value between runs, since `run.sh` resumes from the cached results by default (`RESUME=false` to
recompute everything); `run.sh` prints the directory it uses. The pipeline outputs, which `make stage1` reads,
go to `work/stage1/pipeline`. Neither directory is tracked by git.
