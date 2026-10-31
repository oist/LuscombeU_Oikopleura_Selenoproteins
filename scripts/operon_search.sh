#!/bin/bash -l
#SBATCH --job-name=selenoproteins-operon-search
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.log
# DIAMOND search of the operon-partner proteins (make operon-proteins) against the Swiss-Prot database of the OIST
# Bioinformatics User Group (https://github.com/oist/BioinfoUgrp/blob/master/DataBases.md), outside the container.
# usage, from the repository root:   sbatch -p compute scripts/operon_search.sh
set -euo pipefail
# SPROT_MODULE: the Swiss-Prot DIAMOND database module (DB/diamondDB/sprot/2020_06, 2021-11-28 or 2022-07 on Deigo in
# 2026-10; each directory also holds the NCBI taxonomy files it was built with, for --taxonlist), DIAMOND_MODULE the
# DIAMOND version that reads it
SPROT_MODULE=${SPROT_MODULE:-DB/diamondDB/sprot/2022-07}
DIAMOND_MODULE=${DIAMOND_MODULE:-Other/DIAMOND/2.0.15.153}
ml bioinfo-ugrp-modules
ml "$DIAMOND_MODULE" "$SPROT_MODULE"
echo "modules: $DIAMOND_MODULE $SPROT_MODULE; DIAMONDDB=$DIAMONDDB"
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
make operon-search
