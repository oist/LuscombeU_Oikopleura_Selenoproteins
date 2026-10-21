# Reproduce the analysis: make fetch ; (stage1/run.sh) ; make stage1 ; make curation ; make secis
# Each target writes its outputs under work/ (large, not in git) and small tables under results/tables/ (in git).
SHELL := /bin/bash
GENOMES_BASELINE := OKI OSKA OdB3
GENOMES_COMPARISON := Olon Oalb Ovan Bsty Mery Fbor
GENOMES_NONTARGET := Cint Csav Scla Tdem

.PHONY: all fetch stage1 curation secis secis-stage1 secis-model secis-stage2 stage1-figures machinery machinery-profiles machinery-test machinery-fragments machinery-census machinery-absence machinery-rbh machinery-phylogeny machinery-lba machinery-lba-figures machinery-table check clean
all: stage1 curation secis

fetch:
	scripts/fetch.sh

# Stage 1: the pipeline with its default parameters. Run stage1/run.sh first (Nextflow; see stage1/README.md);
# this target collects its outputs into results/tables/.
stage1:
	@test -d work/stage1/pipeline/compare || { echo "No pipeline output in work/stage1/pipeline: run stage1/run.sh first"; exit 1; }
	python3 stage1/collect.py --outdir work/stage1/pipeline --out-dir results/tables
	mkdir -p results/pipeline_info && cp data/MD5SUMS results/pipeline_info/stage1_input_MD5SUMS.txt
	cp $$(ls -t work/stage1/pipeline/pipeline_info/params_*.json | head -1) results/pipeline_info/stage1_params.json
	cp $$(ls -t work/stage1/pipeline/pipeline_info/execution_trace_*.txt | head -1) results/pipeline_info/stage1_execution_trace.txt
	cp work/stage1/pipeline/pipeline_info/*software_mqc_versions.yml results/pipeline_info/stage1_software_versions.yml


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

# Stage 1 figures on all eighteen genomes (random windows only for the genomes in data/genomes)
stage1-figures:
	python3 figures/stage1_overview.py
	python3 figures/stage1_secis_scores.py --all --table results/tables/stage1_secis_scores_all.tsv --figure figures/stage1_secis_scores_all

# Machinery: appendicularian Selenoprofiles profiles, and their test on the genomes (Selenoprofiles in the image)
#   sbatch --cpus-per-task=4 scripts/run_in_container.sh make machinery
machinery-profiles:
	python3 stage2/build_machinery_profiles.py --out work/stage2/machinery_profiles
	@mkdir -p results/models/machinery_profiles && cp work/stage2/machinery_profiles/*.fa* work/stage2/machinery_profiles/members.tsv results/models/machinery_profiles/ && \
	    tar -czf results/models/appendicularia_machinery_profiles.tar.gz -C work/stage2/machinery_profiles $$(cd work/stage2/machinery_profiles && ls *.fa *.fa.config *.fa.profile_data)

machinery-test:
	python3 stage2/test_machinery_profiles.py --profiles-dir work/stage2/machinery_profiles --out work/stage2/machinery_test
	cp work/stage2/machinery_test/comparison.tsv results/tables/stage2_machinery_profiles_test.tsv

machinery-fragments:
	python3 stage2/machinery_fragments.py --profiles-dir work/stage2/machinery_profiles --test-dir work/stage2/machinery_test --out work/stage2/machinery_fragments
	cp work/stage2/machinery_fragments/fragments.tsv results/tables/stage2_machinery_fragments.tsv

machinery-census:
	python3 stage2/machinery_census.py --test-dir work/stage2/machinery_test --fragments work/stage2/machinery_fragments/fragments.tsv

# Absence tests for machinery genes missing from appendicularian genomes (tblastn, then comparison with related families)
machinery-absence:
	python3 stage2/machinery_absence.py --gene secp43 --outgroup-pattern 'TIA1|TIAR|[Nn]ucleolysin' --outgroup-exclude 'TRNAU1AP|SECp43|tRNA selenocysteine'
	python3 stage2/machinery_absence.py --gene pstk --outgroup-pattern 'KTI12|[Kk]iller toxin.insensitive' --outgroup-faa curation/machinery_outgroups.faa
	for g in secp43 pstk; do for f in summary candidates controls; do cp work/stage2/machinery_absence/$$g.$$f.tsv results/tables/stage2_absence_$$g.$$f.tsv; done; done

# Reciprocal best hits of the machinery predictions against the human RefSeq proteome: stage 1 (standard profiles,
# all genomes) and stage 2 (appendicularian profiles)
machinery-rbh:
	python3 stage2/machinery_rbh.py --table results/tables/stage1_machinery_rbh.tsv
	python3 stage2/machinery_rbh.py --predictions work/stage2/machinery_test/*.faa --table results/tables/stage2_machinery_rbh.tsv

# Phylogenies of PSTK (with KTI12), SEPHS (SPS1/SPS2) and the SECISBP2 family (needs IQ-TREE and ClipKIT)
# Runs in the current shell (no job submission); on a cluster, wrap it in one job, for example
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=24:00:00 scripts/run_in_container.sh make machinery-phylogeny
# The scripts use SLURM_CPUS_PER_TASK threads.
# Starts afresh: the trees in results/trees (committed provisional versions) are removed first; within one run the
# script resumes from the IQ-TREE checkpoints.
machinery-phylogeny:
	rm -rf results/trees/C20 && rm -f results/trees/*.aln results/trees/*.trimmed results/trees/*.treefile results/trees/*.iqtree results/trees/*.ckp.gz results/trees/*.log
	python3 stage2/machinery_phylogeny.py --families pstk SEPHS2 --max-per-group 3
	python3 stage2/machinery_phylogeny.py --families SBP2 --max-per-group 3 --min-fraction 0.2 --clipkit-mode gappy --clipkit-gaps 0.5
	mkdir -p results/trees/C20 && cp results/trees/pstk.trimmed results/trees/SEPHS2.trimmed results/trees/C20/ 2>/dev/null || true
	python3 stage2/machinery_phylogeny.py --families pstk SEPHS2 --max-per-group 3 --model LG+C20+F+G --out results/trees/C20 --figures results/trees/C20
	python3 figures/phylogeny_collapsed.py
	python3 figures/phylogeny_collapsed.py --families pstk SEPHS2 --trees results/trees/C20 --figures results/trees/C20

# Long-branch analyses of the machinery trees (after machinery-phylogeny): PMSF LG+C60, slow-fast, taxon removal,
# AU test; several hours to days on 8 CPUs
# Starts afresh. Run it only after machinery-phylogeny has finished, for example with
#   sbatch --dependency=afterok:<job id of machinery-phylogeny> ... make machinery-lba
# To resume, or to refresh the summaries after a change of the script, run python3 stage2/machinery_lba.py directly:
# finished analyses (IQ-TREE report present) are kept, interrupted ones are redone, and a constrained tree that no
# longer satisfies the constraint is recomputed.
machinery-lba:
	rm -rf results/trees/lba
	python3 stage2/machinery_lba.py
	$(MAKE) machinery-lba-figures machinery-table

# Collapsed display trees from the LG+C60+F+G (PMSF) trees of the long-branch analyses
machinery-lba-figures:
	python3 figures/phylogeny_collapsed.py --trees results/trees/lba --tree-suffix .pmsf --suffix _pmsf --model-label "LG+C60+F+G (PMSF)"

# Table MACHINERY of the manuscript: census, reciprocal best hits, trees and AU tests (seconds; reads results/ only)
machinery-table:
	python3 stage2/machinery_table.py

machinery: machinery-profiles machinery-test machinery-fragments machinery-census machinery-absence

check:
	scripts/check_sizes.sh

clean:
	rm -rf work
