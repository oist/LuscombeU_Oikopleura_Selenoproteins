#!/bin/bash -l
#SBATCH --job-name=selenoproteins-miniprot-v011
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.log
# make miniprot-junctions with miniprot 0.11, the version of the pipeline (module miniprot/align: image
# miniprot:0.11--he4a0461_2), instead of the container's 0.18: same queries (gold-standard proteins), same options
# (-I --gff). miniprot runs from its own image, outside the analysis container; the table is made inside it.
# usage, from the repository root (a clone writable from compute nodes):   sbatch -p compute scripts/miniprot_v011.sh
# environment: MP_SIF = a local miniprot 0.11 image (default: downloaded once into the work directory from MP_URL).
set -euo pipefail
command -v singularity >/dev/null || ml singularity
W=work/stage2/miniprot_junctions_v0.11
MP_URL=${MP_URL:-https://depot.galaxyproject.org/singularity/miniprot:0.11--he4a0461_2}
MP_SIF=${MP_SIF:-$W/miniprot_0.11.sif}
N=${SLURM_CPUS_PER_TASK:-1}
mkdir -p "$W"
[ -s "$MP_SIF" ] || singularity pull "$MP_SIF" "$MP_URL"
echo "miniprot version: $(singularity exec "$MP_SIF" miniprot --version)"
scripts/run_in_container.sh python3 stage2/miniprot_junctions.py --queries-only --work "$W"
for g in OKI OSKA OdB3; do
  fa=$(awk -F, -v g="$g" '$1 == g {print $2}' stage1/samplesheet.csv)
  [ -s "$W/$g.mpi" ] || singularity exec "$MP_SIF" miniprot -t "$N" -d "$W/$g.mpi" "$fa"
  singularity exec "$MP_SIF" miniprot -t "$N" -I --gff "$W/$g.mpi" "$W/queries.faa" > "$W/$g.gff.part"
  mv "$W/$g.gff.part" "$W/$g.gff"
done
scripts/run_in_container.sh python3 stage2/miniprot_junctions.py --from-gff "$W" --work "$W" \
  --table results/tables/stage2_miniprot_junctions_v0.11.tsv
