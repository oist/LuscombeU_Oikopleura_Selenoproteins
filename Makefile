# Reproduce the analysis: make fetch ; (stage1/run.sh) ; make stage1 ; make curation ; make secis
# Each target writes its outputs under work/ (large, not in git) and small tables under results/tables/ (in git).
SHELL := /bin/bash
GENOMES_BASELINE := OKI OSKA OdB3
GENOMES_COMPARISON := Olon Oalb Ovan Bsty Mery Fbor
GENOMES_NONTARGET := Cint Csav Scla Tdem

.PHONY: all fetch stage1 curation secis secis-stage1 secis-model secis-stage2 check clean
all: stage1 curation secis

fetch:
	scripts/fetch.sh

# Stage 1: the pipeline with its default parameters. Run stage1/run.sh first (Nextflow; see stage1/README.md);
# this target collects its outputs into results/tables/.
stage1:
	@test -d work/stage1/pipeline/compare || { echo "No pipeline output in work/stage1/pipeline: run stage1/run.sh first"; exit 1; }
	python3 stage1/collect.py --outdir work/stage1/pipeline --out-dir results/tables

# Baseline: reviewed decisions in curation/*.tsv, checked against stage 1 outputs
curation:
	python3 curation/check_baseline.py

# SECIS (needs Infernal and matplotlib: run in the image, for example
#   sbatch --cpus-per-task=8 scripts/run_in_container.sh make secis)
# Stage 1: RF00031 scores of the Selenoprofiles 3' sequences against random genomic windows
secis-stage1:
	python3 figures/stage1_secis_scores.py --table results/tables/stage1_secis_scores.tsv --figure figures/stage1_secis_scores

# Stage 2: appendicularian SECIS model, bootstrapped from the stage 1 outputs. The default run reproduces the model in
# results/models/; to replicate with other background windows, change BACKGROUND_SEED and SECIS_OUT, for example
#   make secis-model BACKGROUND_SEED=200 SECIS_OUT=work/stage2/secis_model_seed200
SECIS_OUT ?= work/stage2/secis_model
BACKGROUND_SEED ?= 100
secis-model:
	python3 stage2/build_secis_model.py --out $(SECIS_OUT) --background-seed $(BACKGROUND_SEED)
	@if [ "$(SECIS_OUT)" = work/stage2/secis_model ] && [ "$(BACKGROUND_SEED)" = 100 ]; then \
	    mkdir -p results/models && cp $(SECIS_OUT)/appendicularia_secis.cm $(SECIS_OUT)/appendicularia_secis.sto results/models/ && \
	    cp $(SECIS_OUT)/final.seeds.fa results/models/appendicularia_secis.seeds.fa && \
	    cp $(SECIS_OUT)/rounds.tsv results/tables/stage2_secis_model_rounds.tsv; fi

# Scores of the stage 2 model, with new random windows and the curated O. dioica Sec genes
secis-stage2:
	python3 figures/stage1_secis_scores.py --cm $(SECIS_OUT)/appendicularia_secis.cm --model-name 'appendicularian SECIS model' \
	    --baseline curation/baseline_sec_genes.tsv --seed 1000 \
	    --table results/tables/stage2_secis_scores.tsv --figure figures/stage2_secis_scores

secis: secis-stage1 secis-model secis-stage2

check:
	scripts/check_sizes.sh

clean:
	rm -rf work
