#!/usr/bin/env bash
#SBATCH --job-name=selenoproteins-stage1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x-%j.log
# Stage 1: run the selenoprotein pipeline with its default parameters on all the genomes of stage1/samplesheet.csv:
# nine appendicularian assemblies, six other tunicates and three vertebrates (positive controls), in one run.
# This job only drives Nextflow; the tasks run where the Nextflow profile sends them.
# usage, from the repository root:
#   NXF_WORK=<work directory> PROFILE=<profile> stage1/run.sh [extra Nextflow options]
#   NXF_WORK=<work directory> PROFILE=<profile> sbatch [options] stage1/run.sh [extra Nextflow options]
# environment:
#   NXF_WORK  Nextflow work directory, holding the intermediate files of every task. Default: work/stage1/nextflow
#             in the repository. On a cluster, point it to the file system meant for such files (scratch), and
#             keep the same value when resuming.
#   PROFILE   required: container engine or institutional Nextflow profile.
#   PIPELINE_REVISION  default: the release analysed (below).
#   RESUME    default true: the run resumes from the cached results of the previous run (Nextflow -resume), which
#             changes nothing on a first run. Set RESUME=false to recompute every task.
# To run every task in the image built from environment.yml, add: -c conf/pipeline_single_container.config
# Nextflow writes its history (.nextflow/) and log in the repository root and publishes the outputs to
# work/stage1/pipeline: when this script runs as a batch job, the repository must be on a file system that compute
# nodes can write. Otherwise, run it on a login or interactive node.
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
cd "$R"
: "${PROFILE:?set PROFILE to a Nextflow profile (for example singularity, or an institutional profile)}"
PIPELINE_REVISION=${PIPELINE_REVISION:-0.2.0}
export NXF_WORK=${NXF_WORK:-$R/work/stage1/nextflow}
echo "Nextflow work directory (NXF_WORK): $NXF_WORK"
echo "Pipeline outputs: $R/work/stage1/pipeline"
resume=()
if [[ ${RESUME:-true} == true && " $* " != *" -resume "* ]]; then resume=(-resume); fi
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r "$PIPELINE_REVISION" -profile "$PROFILE" \
    -params-file stage1/params.yaml ${resume[@]+"${resume[@]}"} "$@"
