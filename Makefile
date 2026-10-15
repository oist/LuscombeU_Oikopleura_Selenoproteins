# Reproduce the analysis: make fetch ; make stage1 ; make curation ; make stage2 ; make figures
# Each target writes its outputs under work/ (large, not in git) and small tables under results/tables/ (in git).
SHELL := /bin/bash
GENOMES_BASELINE := OKI OSKA OdB3
GENOMES_COMPARISON := Olon Oalb Ovan Bsty Mery Fbor
GENOMES_NONTARGET := Cint Csav Scla Tdem

.PHONY: all fetch stage1 curation stage2 figures check clean
all: stage1 curation stage2 figures

fetch:
	scripts/fetch.sh

# Stage 1: published tools with their default or published parameters
stage1:
	$(MAKE) -C stage1 GENOMES="$(GENOMES_BASELINE) $(GENOMES_COMPARISON)"

# Baseline: reviewed decisions in curation/*.tsv, checked against stage 1 outputs
curation:
	python3 curation/check_baseline.py

# Stage 2: parameters justified by stage 1 and the baseline
stage2:
	$(MAKE) -C stage2 GENOMES="$(GENOMES_BASELINE) $(GENOMES_COMPARISON)" NONTARGET="$(GENOMES_NONTARGET)"

figures:
	$(MAKE) -C figures

check:
	scripts/check_sizes.sh

clean:
	rm -rf work
