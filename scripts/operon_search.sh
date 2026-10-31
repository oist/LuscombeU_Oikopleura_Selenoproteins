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
ml bioinfo-ugrp-modules 2>/dev/null || true
ml Other/DIAMOND/2.0.4.142 DB/diamondDB/sprot/2020_06
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
make operon-search
