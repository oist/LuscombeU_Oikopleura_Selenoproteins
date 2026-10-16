#!/usr/bin/env bash
#SBATCH --job-name=selenoproteins-stage1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x-%j.log
# Stage 1: run the selenoprotein pipeline on the nine appendicularian genomes with its default parameters.
# This job only drives Nextflow; the tasks run where the Nextflow profile sends them.
# usage, from the repository root:
#   PROFILE=<profile> stage1/run.sh [extra Nextflow options]
#   PROFILE=<profile> sbatch [options] stage1/run.sh [extra Nextflow options]       e.g. -resume
# environment: PROFILE (required: container engine or institutional profile); NXF_WORK (work directory, default
#   work/stage1/nextflow); PIPELINE_REVISION (default: the revision analysed below). To run every task in the image
#   built from environment.yml, add: -c conf/pipeline_single_container.config
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
cd "$R"
: "${PROFILE:?set PROFILE to a Nextflow profile (for example singularity, or an institutional profile)}"
PIPELINE_REVISION=${PIPELINE_REVISION:-0.1.0}
export NXF_WORK=${NXF_WORK:-$R/work/stage1/nextflow}
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r "$PIPELINE_REVISION" -profile "$PROFILE" \
    -params-file stage1/params.yaml "$@"
