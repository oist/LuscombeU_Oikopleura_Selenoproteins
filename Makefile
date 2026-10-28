# Reproduce the analysis: make fetch ; (stage1/run.sh) ; make stage1 ; make curation ; make secis
# Each target writes its outputs under work/ (large, not in git) and small tables under results/tables/ (in git).
SHELL := /bin/bash
GENOMES_BASELINE := OKI OSKA OdB3
GENOMES_COMPARISON := Olon Oalb Ovan Bsty Mery Fbor
GENOMES_NONTARGET := Cint Csav Scla Tdem

.PHONY: all fetch stage1 curation secis secis-stage1 secis-model secis-stage2 stage1-figures machinery machinery-profiles machinery-test machinery-fragments machinery-census machinery-absence machinery-rbh machinery-phylogeny machinery-lba machinery-lba-figures machinery-table machinery-scheme reannotation-pages gold-standard secis-hits secis-v2 secis-scan uga-screen secis-description secis-orphans candidate-genes txnl1-homologs txnl1-paralogs txnl1-introns sec-positions tandem-stops selenoprotein-rbh selenoprotein-phylogeny selenoprotein-pmsf selenoprotein-orthology stage2-inputs stage2-truth stage2-benchmark machinery-domains machinery-domains-pmsf machinery-domain-figures check clean
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
	python3 stage2/build_profiles.py --families SEPHS2 pstk SecS eEFsec SBP2 secp43 --gold-genomes \
	    --stage1-genomes OKI OSKA OdB3 Olon Oalb Ovan Bsty Mery Fbor --held-out OdB3 --out work/stage2/machinery_profiles
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

# Domain analysis of the machinery trees (after machinery-phylogeny): one tree per domain (domains defined on human
# proteins, curation/machinery_domains.tsv), supports and AU tests of each domain's tree against the other domain's data
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=1-00:00:00 scripts/run_in_container.sh make machinery-domains
machinery-domains:
	rm -rf results/trees/domains
	python3 stage2/machinery_domains.py
	$(MAKE) machinery-domain-figures

# The same under LG+C60+F+G with PMSF (after machinery-domains, whose maximum-likelihood trees are the guide trees;
# keeps finished analyses, so it can be resubmitted after an interruption)
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=2-00:00:00 scripts/run_in_container.sh make machinery-domains-pmsf
machinery-domains-pmsf:
	python3 stage2/machinery_domains.py --pmsf --table results/tables/stage2_machinery_domains_pmsf.tsv \
	  --conflicts results/tables/stage2_machinery_domain_conflicts_pmsf.tsv
	$(MAKE) machinery-domain-figures DOMAIN_SUFFIX=.pmsf

# Collapsed display trees of the domains (PDF only; seconds); DOMAIN_SUFFIX=.pmsf for the PMSF trees
DOMAIN_TREES = pstk:kinase pstk:C_terminal SEPHS2:N_terminal SEPHS2:C_terminal SBP2:upstream_of_L7Ae SBP2:L7Ae
DOMAIN_SUFFIX =
machinery-domain-figures:
	for x in $(DOMAIN_TREES); do f=$${x%%:*}; d=$${x#*:}; \
	  if [ -n "$(DOMAIN_SUFFIX)" ]; then m="LG+C60+F+G (PMSF)"; \
	  else m=$$(sed -n 's/^Best-fit model according to BIC: //p' results/trees/domains/$$f.$$d.iqtree); fi; \
	  l=$$(echo $$d | tr _ ' ' | sed 's/^\([CN]\) /\1-/'); \
	  python3 figures/phylogeny_collapsed.py --families $$f --trees results/trees/domains --tree-suffix .$$d$(DOMAIN_SUFFIX) \
	    --suffix _domain_$${d}$(subst .,_,$(DOMAIN_SUFFIX)) --model-label "$$l domain, $$m" --formats pdf || exit 1; \
	done

# Step 2 of the plan of 2026-10-24: gold-standard gene models of the curated O. dioica selenoprotein genes in OKI2018,
# OSKA2016 and OdB3 (transcripts aligned with minimap2, cached in work/stage2/gold_standard; minutes)
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=4:00:00 scripts/run_in_container.sh make gold-standard
gold-standard:
	python3 curation/gold_standard.py

# Step 1 of the plan of 2026-10-24: SECIS scores of the validated selenocysteine hits (seconds; reads results/ only)
secis-hits:
	python3 stage2/secis_validated_hits.py

# Step 3 of the plan of 2026-10-24: SECIS model trained on the six Oikopleura assemblies (gold-standard stop codons of
# step 2, stage 1 predictions of O. longicauda, O. albicans and O. vanhoeffeni), evaluated leaving out one family or
# one species at a time against RF00031, RF00031 merged and the current model (a few hours), for example
#   sbatch -p compute --cpus-per-task=16 --mem=32G --time=8:00:00 scripts/run_in_container.sh make secis-v2
secis-v2:
	python3 stage2/secis_model_training.py --jobs 4

# Step 4 of the plan of 2026-10-24: UGA/SECIS screens. secis-scan: genome-wide SECIS scans (Oikopleura model) of the
# six Oikopleura assemblies, the three held-out appendicularians (B. stygius, M. erythrocephalus, F. borealis; for the
# final census) and four non-target tunicates, and of window-shuffled copies, and SECIS searches of the
# transcript sets (O. dioica annotations and transcript assemblies; human and C. intestinalis RefSeq transcripts as
# controls, also with RF00031); one to a few hours, for example
#   sbatch -p compute --cpus-per-task=16 --mem=32G --time=12:00:00 scripts/run_in_container.sh make secis-scan
# uga-screen: genome-scan threshold, transcript and comparative screens, blastp, tables (minutes); it stops when an
# output of secis-scan is missing. Both at once, the second waiting for the first:
#   jid=$(sbatch --parsable -p compute --cpus-per-task=16 --mem=32G --time=12:00:00 scripts/run_in_container.sh make secis-scan)
#   sbatch --dependency=afterok:$jid -p compute --cpus-per-task=8 --mem=16G --time=4:00:00 scripts/run_in_container.sh make uga-screen
secis-scan:
	python3 stage2/uga_secis_screen.py scan --jobs 4

uga-screen:
	python3 stage2/uga_secis_screen.py screen

# Step 5: candidate genes outside the baseline (curation/candidates.tsv and .faa: TXNL1.b, the KKSU family): the rule
# fixed in advance in the three O. dioica assemblies, tblastn against the appendicularian, tunicate and amphioxus
# genomes (codon at the selenocysteine position, SECIS after the stop codon) and blastp against RefSeq proteomes
# (about an hour), for example
#   sbatch -p compute --cpus-per-task=16 --mem=32G --time=6:00:00 scripts/run_in_container.sh make candidate-genes
# Description of the Oikopleura SECIS model (seeds, base pairs against decoys, changes between O. dioica orthologs)
# and its scores in genomes and controls (reads repository files only; about 7 minutes with 8 CPUs)
secis-description:
	python3 stage2/secis_model_description.py
	python3 figures/secis_model_figures.py

# Genome-scan SECIS hits at the genome threshold without a selenoprotein gene upstream (classes of
# results/tables/stage2_secis_scan_classes.tsv): annotation, blastx of the upstream region, copies in the genome,
# conservation in the other Oikopleura assemblies (genome databases built in work/stage2/secis_orphans and removed)
secis-orphans:
	python3 stage2/secis_orphan_hits.py

candidate-genes:
	python3 stage2/candidate_genes.py --jobs 4

# TXNL1 homologs in the three O. dioica assemblies: reciprocal best hits with human TXNL1, next hits of O. dioica
# TXNL1.b in its own assembly, hits of the thioredoxin and PITH parts of human TXNL1 (a few minutes), for example
#   sbatch -p compute --cpus-per-task=4 --mem=16G --time=1:00:00 scripts/run_in_container.sh make txnl1-homologs
txnl1-homologs:
	python3 stage2/txnl1_homologs.py

# The two TXNL1 genes of O. dioica (TXNL1.a, full-length, and TXNL1.b, selenocysteine, PITH region only) in the
# other genomes and proteomes (structure, codon at the selenocysteine position, SECIS), a tree of the PITH regions
# rooted on PITHD1, and the neighbourhood of the full-length copy (GPR107-like gene; about an hour), for example
#   sbatch -p compute --cpus-per-task=16 --mem=32G --time=6:00:00 scripts/run_in_container.sh make txnl1-paralogs
txnl1-paralogs:
	python3 stage2/txnl1_paralogs.py

# Position of the selenocysteine codon in the gold-standard models (exon junctions, last exon, distance to the stop
# codon), against the other codons of each gene (repository files only; seconds)
sec-positions:
	python3 stage2/sec_codon_positions.py

# Genes ending with a TGA followed closely by another stop codon (a path to a C-terminal selenoprotein), and
# SECIS-like hits after them in the genome scans of step 4 (needs work/stage2/uga_screen/scan; under an hour), e.g.
#   sbatch -p compute --cpus-per-task=4 --mem=32G --time=4:00:00 scripts/run_in_container.sh make tandem-stops
tandem-stops:
	python3 stage2/tandem_stops.py

# Intron positions and phases of TXNL1.b and TXNL1.a on the TXNL1.a protein (repository files only; seconds)
txnl1-introns:
	python3 stage2/txnl1_introns.py


# Orthology of the appendicularian selenoproteins, with the approach used for the machinery: reciprocal best hits,
# trees, LG+C60+F+G (PMSF) trees, then a first-pass summary per appendicularian clade (stage2/README.md)
SELENO_TREES = DIO GPX MsrA SELENOF SELENOK SELENOL SELENOO SELENOS SELENOT SELENOU SELENOW TXNRD

# Reciprocal best hits of the stage 1 selenoprotein predictions (all genomes) against the human RefSeq proteome
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=12:00:00 scripts/run_in_container.sh make selenoprotein-rbh
selenoprotein-rbh:
	python3 stage2/machinery_rbh.py --families selenoproteins --table results/tables/stage2_selenoprotein_rbh.tsv

# Trees of the families of SELENO_TREES that have appendicularian sequences (stage 1 predictions; human RefSeq
# proteins of each family as references); starts afresh, removing results/trees/selenoproteins
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=2-00:00:00 scripts/run_in_container.sh make selenoprotein-phylogeny
selenoprotein-phylogeny:
	rm -rf results/trees/selenoproteins
	python3 stage2/machinery_phylogeny.py --families $(SELENO_TREES) --max-per-group 3 --min-fraction 0.3 \
	  --stage2 work/stage1/pipeline/selenoprofiles --out results/trees/selenoproteins --figures figures/selenoproteins --formats pdf
	python3 figures/phylogeny_collapsed.py --families $(SELENO_TREES) --trees results/trees/selenoproteins \
	  --groups-dir results/trees/selenoproteins --figures figures/selenoproteins --formats pdf

# LG+C60+F+G (PMSF) trees of the same families, after selenoprotein-phylogeny (keeps finished trees, so it can be
# resubmitted), for example
#   sbatch --dependency=afterok:<job id> -p compute --cpus-per-task=8 --mem=32G --time=2-00:00:00 scripts/run_in_container.sh make selenoprotein-pmsf
selenoprotein-pmsf:
	python3 stage2/machinery_lba.py --families $(SELENO_TREES) --trees results/trees/selenoproteins --analyses pmsf \
	  --out results/trees/selenoproteins/lba --table results/tables/stage2_selenoprotein_lba.tsv
	python3 figures/phylogeny_collapsed.py --families $(SELENO_TREES) --trees results/trees/selenoproteins/lba \
	  --groups-dir results/trees/selenoproteins --tree-suffix .pmsf --suffix _pmsf --model-label "LG+C60+F+G (PMSF)" \
	  --figures figures/selenoproteins --formats pdf

# First-pass orthology of each appendicularian clade, from the trees and best hits (seconds; reads results/ only)
selenoprotein-orthology:
	python3 stage2/selenoprotein_orthology.py

# Schematic of the selenocysteine machinery (figures/machinery_scheme.{pdf,svg}; seconds, no input)
machinery-scheme:
	python3 figures/machinery_scheme.py

# Figure REANNOTATION (two A4 pages), from the gold standard and figures/inputs (written by
# figures/reannotation_inputs.sh, which needs www.ebi.ac.uk and the genomes). Genes grouped by family, families kept
# together in reading order: page 1, glutathione peroxidases, deiodinases, then thioredoxin reductase, SELENOO and
# SEPHS2; page 2, the Rdx family (SELENOW paralogs and SELENOT), then SELENOF, SELENOK and SELENOS.
REANNOTATION_PAGE1 = GPX_PAR GPX_chr1a GPX_chr1b GPX_chr2 DIO_PAR DIO_chr1 TXNRD SELENOO SEPHS2
REANNOTATION_PAGE2 = SELENOW_PAR SELENOW_chr1a SELENOW_chr2 SELENOW_pair1 SELENOW_pair2 SELENOT SELENOF SELENOK SELENOS
# genes added to the baseline in step 5 (2026-10-26), drawn in their own figure in the same style; TXNL1.b with its
# paralog TXNL1.a (figures/inputs/reannotation_TXNL1.b.paralogs.*, written by make txnl1-paralogs)
REANNOTATION_CANDIDATES = TXNL1.b KKSU
reannotation-pages:
	python3 figures/reannotation_page.py --genes $(REANNOTATION_PAGE1) --name reannotation_page
	python3 figures/reannotation_page.py --genes $(REANNOTATION_PAGE2) --name reannotation_page2
	python3 figures/reannotation_page.py --genes $(REANNOTATION_CANDIDATES) --name reannotation_candidates --fit \
	  --sources "human: UniProt; C. intestinalis: RefSeq; O. albicans: open stretch of the comparative screen on its genome (accessions in the paragraphs)"

# Table MACHINERY of the manuscript: census, reciprocal best hits, trees and AU tests (seconds; reads results/ only)
machinery-table:
	python3 stage2/machinery_table.py

# Stage 2 (steps 6 and 7 of the plan of 2026-10-24). Inputs, in the image (Selenoprofiles; about 1 h on 4 CPUs):
#   sbatch [--partition=<partition>] --cpus-per-task=4 scripts/run_in_container.sh make stage2-inputs stage2-truth
# - profiles: E0 (gold standard of the three O. dioica assemblies; machinery from their stage 1 predictions), LOO_<g>
#   (the same without one assembly), final (also the validated stage 1 predictions of the three other Oikopleura
#   species, and the curated proteins of other species); B. stygius, M. erythrocephalus and F. borealis held out;
# - the Oikopleura SECIS model with a gathering threshold of 20 bits (the pipeline then searches it with --cut_ga).
# Then the pipeline runs (stage2/run/run.sh <run>; see stage2/README.md), then make stage2-benchmark.
STAGE2_PROFILES = work/stage2/profiles
STAGE2_BUILDS = E0 LOO_OKI LOO_OSKA LOO_OdB3 final
stage2-inputs:
	@mkdir -p work/stage2/models
	python3 stage2/cm_thresholds.py --ga 20 results/models/oikopleura_secis.cm work/stage2/models/oikopleura_secis.cm
	python3 stage2/build_profiles.py --stage1-genomes OKI OSKA OdB3 --out $(STAGE2_PROFILES)/E0
	for g in OKI OSKA OdB3; do \
	    python3 stage2/build_profiles.py --stage1-genomes OKI OSKA OdB3 --held-out Bsty Mery Fbor $$g \
	        --out $(STAGE2_PROFILES)/LOO_$$g || exit 1; done
	python3 stage2/build_profiles.py --stage1-genomes OKI OSKA OdB3 Olon Oalb Ovan --validated-genomes Olon Oalb Ovan \
	    --candidates curation/candidates.faa --out $(STAGE2_PROFILES)/final \
	    --archive work/stage2/appendicularia_selenoprofiles.tar.gz
	@mkdir -p results/models/stage2_profiles && rm -f results/models/stage2_profiles/*.fa results/models/stage2_profiles/*.fa.config && \
	    cp $(STAGE2_PROFILES)/final/*.fa $(STAGE2_PROFILES)/final/*.fa.config results/models/stage2_profiles/
	awk 'NR == 1 {print "build\t" $$0} FNR > 1 {n = split(FILENAME, p, "/"); print p[n - 1] "\t" $$0}' \
	    $(foreach b,$(STAGE2_BUILDS),$(STAGE2_PROFILES)/$(b)/members.tsv) > results/tables/stage2_profile_members.tsv

# Reference sets: gold standard, gold models lifted to the new assemblies of the same species (minimap2), validated
# stage 1 predictions of the other Oikopleura species (minutes)
stage2-truth:
	python3 stage2/truth_sets.py

# Stage 1 and stage 2 settings against the reference sets (seconds; reads the pipeline outputs of each run present)
stage2-benchmark:
	python3 stage2/benchmark.py
	@mkdir -p results/pipeline_info && for d in work/stage2/pipeline/*/pipeline_info; do \
	    r=$$(basename $$(dirname $$d)); \
	    cp $$(ls -t $$d/params_*.json | head -1) results/pipeline_info/stage2_$${r}_params.json; \
	    cp $$(ls -t $$d/execution_trace_*.txt | head -1) results/pipeline_info/stage2_$${r}_execution_trace.txt; \
	    cp $$d/*software_mqc_versions.yml results/pipeline_info/stage2_$${r}_software_versions.yml; done

machinery: machinery-profiles machinery-test machinery-fragments machinery-census machinery-absence

check:
	scripts/check_sizes.sh

clean:
	rm -rf work
