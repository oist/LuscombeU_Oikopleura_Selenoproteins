#!/usr/bin/env bash
#SBATCH --job-name=analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.log
# Run a command of this repository inside the image built from environment.yml (Singularity or Apptainer).
# usage, from the repository root:   scripts/run_in_container.sh make stage1
#                                    sbatch [options] scripts/run_in_container.sh make curation
# environment: IMAGE = a local .sif file, or a docker:// or oras:// address (default below). An https:// address of
# Seqera Containers can be downloaded first: curl -L -o analysis.sif <address>; IMAGE=analysis.sif
# Directories outside the repository that data/ links point to must be bound: export SINGULARITY_BIND=<dir>
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
IMAGE=${IMAGE:-PLACEHOLDER_DOCKER_IMAGE_ADDRESS}
[[ $IMAGE == PLACEHOLDER* ]] && { echo "Set IMAGE: the image address is not filled in yet" >&2; exit 1; }
RUNNER=$(command -v singularity || command -v apptainer) || { echo "singularity or apptainer not found" >&2; exit 1; }
cd "$R"
exec "$RUNNER" exec --bind "$R" "$IMAGE" "$@"
