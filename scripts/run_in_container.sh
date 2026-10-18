#!/usr/bin/env bash
#SBATCH --job-name=selenoproteins-analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.log
# Run a command of this repository inside the image built from environment.yml (Singularity or Apptainer).
# usage, from the repository root:   scripts/run_in_container.sh make stage1
#                                    sbatch [options] scripts/run_in_container.sh make curation
# The job writes its log (%x-%j.log) in the submission directory and its results in the repository: both must be
# writable from compute nodes (at OIST, submit from a clone on /flash, not from /bucket; see README).
# environment: IMAGE = a local .sif file, or a docker:// or oras:// address (default below: the image built from
# environment.yml, converted by Singularity on first use). The same image as a ready-made Singularity file:
#   curl -L -o oikopleura-selenoproteins.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7a5982e08142390038579f334f8683c3a4cd9c922d0c4e84227b5eb2c65d421/data
#   export IMAGE=$PWD/oikopleura-selenoproteins.sif
# Links from data/ or work/<stage>/<dir> to other directories are bound automatically; other directories outside the
# repository must be bound with SINGULARITY_BIND=<dir>.
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
# The image 60a1e963fa7d21f6 predates IQ-TREE and ClipKIT (machinery phylogenies): rebuild it from environment.yml
# with Wave before running make machinery-phylogeny in the container.
IMAGE=${IMAGE:-docker://community.wave.seqera.io/library/oikopleura-selenoproteins:60a1e963fa7d21f6}
RUNNER=$(command -v singularity || command -v apptainer) || {
    echo "singularity or apptainer not found: load it before running or submitting this script (at OIST: ml singularity;" >&2
    echo "sbatch passes the loaded environment to the job)" >&2; exit 1; }
cd "$R"
# bind the repository and the targets of links to inputs kept elsewhere (for example a clone on /flash whose data/
# and work/stage1/pipeline link to the main copy)
binds=$R
for l in data work/*/*; do
    if [[ -L $l ]]; then binds+=",$(readlink -f "$l")"; fi
done
exec "$RUNNER" exec --bind "$binds" "$IMAGE" "$@"
