#!/usr/bin/env bash
#SBATCH --job-name=tiara
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.log
# Tiara classification (eukarya, bacteria, archaea, prokarya, organelle, unknown) of the contigs of each genome, for
# the contamination screen of the census (step 9; read by stage2/census_contigs.py). Runs outside the image of
# environment.yml, whose Python is too recent for Tiara: TIARA_IMAGE is a Singularity image or address with Tiara,
# for example one built by Wave from conda-forge, as for environment.yml:
#   wave --conda-package conda-forge::tiara=1.0.3 --freeze --await --singularity
# usage, from the repository root (at OIST, ml singularity first):
#   TIARA_IMAGE=<image> sbatch --partition=compute scripts/run_tiara.sh OKI OSKA OdB3 KUM AOM Bar2 Olon Oalb Ovan Bsty Mery Fbor
# Outputs work/stage2/census/tiara/<genome>.tiara.tsv (genomes already done are skipped). Tiara classifies sequences
# of at least 3,000 nt (its default).
set -euo pipefail
: "${TIARA_IMAGE:?set TIARA_IMAGE to a Singularity image or address with Tiara}"
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
RUNNER=$(command -v singularity || command -v apptainer) || { echo "singularity or apptainer not found (at OIST: ml singularity)" >&2; exit 1; }
cd "$R"
out=work/stage2/census/tiara
mkdir -p "$out"
tmp=$(mktemp -d "$R/work/stage2/census/tmp.XXXXXX")
trap 'rm -f "${tmp:?}"/*.fa; rmdir "${tmp:?}"' EXIT
binds=$R
for l in data data/genomes/*; do
    if [[ -L $l ]]; then binds+=",$(readlink -f "$l")"; fi
done
for g in "$@"; do
    if [[ -s $out/$g.tiara.tsv ]]; then echo "$g: done before"; continue; fi
    zcat -f "data/genomes/$g/genome.fna.gz" > "$tmp/$g.fa"
    "$RUNNER" exec --bind "$binds" "$TIARA_IMAGE" tiara -i "$tmp/$g.fa" -o "$out/$g.tiara.tsv.part" \
        -t "${SLURM_CPUS_PER_TASK:-1}"
    mv "$out/$g.tiara.tsv.part" "$out/$g.tiara.tsv"
    rm -f "$tmp/$g.fa"
    echo "$g: $(($(wc -l < "$out/$g.tiara.tsv") - 1)) sequences classified"
done
