# Stage 1 — the pipeline with its default parameters

Stage 1 runs the selenoprotein pipeline (`oist/LuscombeU_Selenoproteins_Pipeline`, release 0.1.0, set in `run.sh`) on
the nine appendicularian genomes, with the parameters a new user would get:

- Selenoprofiles 4.5.8, profiles `eukarya,machinery`;
- miniprot alignment of the Sec-containing profile sequences, Sec codons by Sec-context matching;
- Rfam RF00031 (SECIS) on the 3′ sequence of each Selenoprofiles prediction, at the Rfam gathering threshold;
- Secmarker and cmsearch (Rfam RF01852) for tRNA-Sec;
- for the three *O. dioica* assemblies, the annotation check (classes of Ticó et al. 2026) and the transcript
  support of each Sec codon (OKI2018 transcriptome, OSKA2016 transcripts, OdB3 cDNA).

| File | Contents |
|---|---|
| `samplesheet.csv` | genomes, annotations and transcripts, as paths under `data/` (from `scripts/fetch.sh`) |
| `params.yaml` | pipeline parameters: local copies of the Selenoprofiles profiles and of the two Rfam models |
| `run.sh` | Nextflow driver (SBATCH headers); `PROFILE=<profile> stage1/run.sh` or `PROFILE=<profile> sbatch stage1/run.sh` |
| `collect.py` | pipeline outputs → `results/tables/stage1_*.tsv` (run by `make stage1`) |

Run from the repository root:

```bash
PROFILE=singularity stage1/run.sh            # or an institutional profile; add -resume to continue a run
make stage1                                  # collect the tables
make curation                                # compare with the baseline
```

The pipeline outputs go to `work/stage1/pipeline` and its work directory to `work/stage1/nextflow` (override with
`NXF_WORK`); neither is tracked by git.
