# Stage 1: published tools and parameters. Called from the top-level Makefile with GENOMES="...".
SHELL := /bin/bash
R := $(abspath ..)
GENOMES ?= OKI OSKA OdB3 Olon Oalb Ovan Bsty Mery Fbor
DONE := $(addprefix $(R)/work/stage1/,$(addsuffix /done,$(GENOMES)))

.PHONY: all
all: $(R)/results/tables/stage1_sec_sites.tsv

$(R)/work/stage1/queries.faa:
	mkdir -p $(dir $@)
	python3 bin/profile_queries.py $(R)/data/references/selenoprofiles_data $@

$(R)/work/stage1/%/done: $(R)/work/stage1/queries.faa
	./run_genome.sh $* && touch $@

$(R)/results/tables/stage1_sec_sites.tsv: $(DONE)
	python3 bin/collect.py --work $(R)/work/stage1 --genomes $(GENOMES) --out-dir $(R)/results/tables
