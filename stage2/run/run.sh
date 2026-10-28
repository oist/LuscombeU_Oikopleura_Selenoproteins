#!/usr/bin/env bash
#SBATCH --job-name=selenoproteins-stage2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x-%j.log
# Stage 2: run the selenoprotein pipeline with the settings of one stage 2 run (stage2/run/<run>.yaml):
#   default        stage 1 settings on the three new assemblies (KUM, AOM, Bar2)
#   eval           stage 2 settings, evaluation profiles E0, twelve appendicularian assemblies
#   eval_miniprot  as eval, miniprot with -C 0.5 (run it after eval: Nextflow reuses the Selenoprofiles tasks)
#   loo_OKI, loo_OSKA, loo_OdB3   profiles without one gold-standard assembly, on that assembly and its species' new one
#   final          final profiles, for the census
# Inputs made by "make stage2-inputs" (profiles in work/stage2/profiles/, SECIS model in work/stage2/models/) and
# "make fetch". This job only drives Nextflow; the tasks run where the Nextflow profile sends them.
# usage, from the repository root:
#   NXF_WORK=<work directory> PROFILE=<profile> stage2/run/run.sh <run> [extra Nextflow options]
#   NXF_WORK=<work directory> PROFILE=<profile> sbatch [options] stage2/run/run.sh <run> [extra Nextflow options]
# environment: as stage1/run.sh (NXF_WORK, default work/stage2/nextflow; PROFILE, required; PIPELINE_REVISION;
# RESUME, default true).
set -euo pipefail
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/../.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
cd "$R"
RUN=${1:?usage: stage2/run/run.sh <run> [extra Nextflow options]; runs: the names of stage2/run/*.yaml}
shift
[[ -f stage2/run/$RUN.yaml ]] || { echo "no stage2/run/$RUN.yaml" >&2; exit 1; }
: "${PROFILE:?set PROFILE to a Nextflow profile (for example singularity, or an institutional profile)}"
PIPELINE_REVISION=${PIPELINE_REVISION:-0.2.0}
export NXF_WORK=${NXF_WORK:-$R/work/stage2/nextflow}
configs=()
case $RUN in
    default) ;;
    eval_miniprot) configs=(-c stage2/run/stage2.config -c stage2/run/miniprot_C05.config) ;;
    *) configs=(-c stage2/run/stage2.config) ;;
esac
for f in $(awk '$1 == "secis_cm:" || $1 == "selenoprofiles_overrides:" {print $2}' "stage2/run/$RUN.yaml"); do
    [[ -e $f ]] || { echo "missing input $f: run make stage2-inputs first" >&2; exit 1; }
done
echo "Run: $RUN; Nextflow work directory (NXF_WORK): $NXF_WORK"
echo "Pipeline outputs: $(awk '$1 == "outdir:" {print $2}' "stage2/run/$RUN.yaml")"
resume=()
if [[ ${RESUME:-true} == true && " $* " != *" -resume "* ]]; then resume=(-resume); fi
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r "$PIPELINE_REVISION" -profile "$PROFILE" \
    -params-file "stage2/run/$RUN.yaml" ${configs[@]+"${configs[@]}"} ${resume[@]+"${resume[@]}"} "$@"
