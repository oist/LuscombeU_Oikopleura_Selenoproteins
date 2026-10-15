#!/usr/bin/env bash
#SBATCH --job-name=analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.log
# Run a command of this repository inside the image built from environment.yml (Singularity or Apptainer).
# usage, from the repository root:   scripts/run_in_container.sh make stage1
#                                    sbatch [options] scripts/run_in_container.sh make curation
# environment: IMAGE = a local .sif file, or a docker:// or oras:// address (default below: the image built from
# environment.yml, converted by Singularity on first use). The same image as a ready-made Singularity file:
#   curl -L -o oikopleura-selenoproteins.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7a5982e08142390038579f334f8683c3a4cd9c922d0c4e84227b5eb2c65d421/data
#   export IMAGE=$PWD/oikopleura-selenoproteins.sif
# Directories outside the repository that data/ links point to must be bound: export SINGULARITY_BIND=<dir>
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
IMAGE=${IMAGE:-docker://community.wave.seqera.io/library/oikopleura-selenoproteins:60a1e963fa7d21f6}
RUNNER=$(command -v singularity || command -v apptainer) || { echo "singularity or apptainer not found" >&2; exit 1; }
cd "$R"
exec "$RUNNER" exec --bind "$R" "$IMAGE" "$@"
