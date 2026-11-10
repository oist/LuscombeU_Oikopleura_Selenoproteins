# Reproduce the analysis: make fetch ; (stage1/run.sh) ; make stage1 ; make curation ; make secis
# Each target writes its outputs under work/ (large, not in git) and small tables under results/tables/ (in git).
SHELL := /bin/bash
GENOMES_BASELINE := OKI OSKA OdB3
GENOMES_COMPARISON := Olon Oalb Ovan Bsty Mery Fbor
GENOMES_NONTARGET := Cint Csav Scla Tdem

.PHONY: all fetch stage1 curation secis secis-stage1 secis-model secis-stage2 stage1-figures machinery machinery-profiles machinery-test machinery-fragments machinery-census machinery-absence machinery-rbh gold-standard-machinery machinery-phylogeny machinery-lba machinery-lba-figures machinery-table machinery-clade-supports machinery-reannotation machinery-reannotation-inputs sl-transcripts machinery-regions region-proteins pstk-regions pstk-one-leaf pstk-no-outgroup pstk-enriched pstk-trees pstk-lba-domains pstk-side-trees pstk-summaries sbp2-trees sbp2-lba-domains sbp2-summaries sbp2-figures pstk-introns pstk-rates intron-map junction-rescan cds-divergence intron-presence intron-pairing pstk-structures pstk-intron-scheme pstk-intron-columns intron-content eefsec-insertion machinery-scheme reannotation-pages gold-standard secis-hits secis-v2 secis-scan uga-screen secis-description secis-orphans candidate-genes txnl1-homologs txnl1-paralogs txnl1-introns sec-positions tandem-stops selenoprotein-rbh selenoprotein-phylogeny selenoprotein-pmsf selenoprotein-au selenoprotein-orthology selenoprotein-table selenoprotein-absence stage2-inputs stage2-truth stage2-benchmark stage2-diagnose selenok-diagnose selenoi-regions selenok-regions selenok-models stage2-release-profiles selenoi-secis operons operons-figure human-uniprot operon-proteins operon-search operon-partners miniprot-junctions miniprot-junctions-pipeline miniprot-junctions-v011-table ms-peptides selenoi-junction stage2-bundle stage2-collect census-contigs census machinery-domains machinery-domains-pmsf machinery-domain-figures check clean
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
	python3 stage2/machinery_rbh.py --predictions $(SELENO_PREDICTIONS)/*.selenoprofiles.faa --table results/tables/stage2_machinery_rbh.tsv

# Phylogenies of PSTK (with KTI12), SEPHS (SPS1/SPS2) and the SECISBP2 family (needs IQ-TREE and ClipKIT)
# Inputs (decision of 2026-11-01): the gold-standard proteins of the three O. dioica assemblies (selenoprotein gold
# standard for SEPHS2, machinery gold standard for the other genes; make gold-standard-machinery first) and the
# predictions of the final stage 2 run for the other appendicularians; B. stygius and M. erythrocephalus (held out,
# fragmented assemblies) left out, F. borealis kept (the only fritillariid).
# Runs in the current shell (no job submission); on a cluster, wrap it in one job, for example
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=24:00:00 scripts/run_in_container.sh make machinery-phylogeny
# The scripts use SLURM_CPUS_PER_TASK threads.
# Starts afresh: the trees in results/trees (committed provisional versions) are removed first; within one run the
# script resumes from the IQ-TREE checkpoints.
# SECIS-binding proteins (2026-11-09): deuterostome proteins only, rooted on mollusc and annelid proteins, with the
# proteins of curation/extra_tree_proteins.faa (more ascidians, a third amphioxus, shark, frog, bird, lampreys and
# hagfishes: BLASTP searches by Charles; 2026-11-10: both paralogs of chicken, Xenopus tropicalis and coelacanth, two
# hemichordates, two echinoderms); the profile sequences of other groups (non-animals, ecdysozoans, flatworms,
# unclassified accessions) and the truncated proteins of curation/tree_exclusions.tsv are left out (recorded in
# results/trees/SBP2.removed.tsv)
SBP2_TREE_OPTS = --max-per-group 3 --min-fraction 0.2 --clipkit-mode gappy --clipkit-gaps 0.5 \
  --keep-groups Mammals Birds Reptiles Amphibians "Ray-finned fishes" "Lobe-finned fishes" "Cartilaginous fishes" \
  "Jawless fishes" Cephalochordates "Other tunicates" Appendicularians Echinoderms Hemichordates Molluscs Annelids \
  Brachiopods --outgroup-groups Molluscs Annelids Brachiopods
machinery-phylogeny:
	rm -rf results/trees/C20 && rm -f results/trees/*.aln results/trees/*.trimmed results/trees/*.treefile results/trees/*.iqtree results/trees/*.ckp.gz results/trees/*.log
	python3 stage2/machinery_phylogeny.py --families pstk SEPHS2 --max-per-group 3 $(MACHINERY_TREE_INPUTS)
	python3 stage2/machinery_phylogeny.py --families SBP2 $(SBP2_TREE_OPTS) $(MACHINERY_TREE_INPUTS)
	mkdir -p results/trees/C20 && cp results/trees/pstk.trimmed results/trees/SEPHS2.trimmed results/trees/C20/ 2>/dev/null || true
	python3 stage2/machinery_phylogeny.py --families pstk SEPHS2 --max-per-group 3 --model LG+C20+F+G --out results/trees/C20 --figures results/trees/C20 $(MACHINERY_TREE_INPUTS)
	python3 figures/phylogeny_collapsed.py
	python3 figures/phylogeny_collapsed.py --families pstk SEPHS2 --trees results/trees/C20 --figures results/trees/C20

# PSTK tree with one leaf per appendicularian species (robustness of the appendicularian placement to related leaves):
# OdB3 as the only O. dioica assembly, without the second O. vanhoeffeni locus (pseudogene); maximum likelihood and
# LG+C20+F+G (hours; after or alongside machinery-phylogeny)
pstk-one-leaf:
	rm -rf results/trees/pstk_one_leaf && mkdir -p results/trees/pstk_one_leaf
	python3 stage2/machinery_phylogeny.py --families pstk --max-per-group 3 $(MACHINERY_TREE_INPUTS) \
	  --exclude-sequences '^OKI_' '^OSKA_' 'pseudo' --out results/trees/pstk_one_leaf --figures results/trees/pstk_one_leaf
	mkdir -p results/trees/pstk_one_leaf/C20 && cp results/trees/pstk_one_leaf/pstk.trimmed results/trees/pstk_one_leaf/C20/
	python3 stage2/machinery_phylogeny.py --families pstk --max-per-group 3 --model LG+C20+F+G $(MACHINERY_TREE_INPUTS) \
	  --exclude-sequences '^OKI_' '^OSKA_' 'pseudo' --out results/trees/pstk_one_leaf/C20 --figures results/trees/pstk_one_leaf/C20

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
# PSTK tree without its outgroup KTI12 (2026-11-04): where the appendicularian branch goes when the outgroup cannot
# attract it; local, about 15 min on 2 CPUs (needs results/trees/pstk.trimmed and .treefile)
# PSTK trees with more eukaryotic PSTK and with archaeal PSTK (2026-11-04): sequences of the BLASTP searches of
# curation/pstk_blast selected by stage2/pstk_enrichment.py (selection.tsv), trees aside from the main one; local,
# a few hours: run on Deigo (sbatch, 8 CPUs), not in a session container; then make machinery-clade-supports
pstk-enriched:
	rm -rf results/trees/pstk_enriched results/trees/pstk_enriched_archaea
	python3 stage2/pstk_enrichment.py --out results/trees/pstk_enriched

# PSTK analyses again after a change of the PSTK sequence set only (2026-11-07: PSTK of B. floridae and B. belcheri
# added from curation/extra_tree_proteins.faa, see stage2/machinery_phylogeny.py --extra-proteins). The PSTK
# outputs of each analysis are removed; the scripts keep the finished analyses of SEPHS2 and SBP2, so the tables stay
# complete. On Deigo: pstk-trees (hours), then pstk-one-leaf, pstk-lba-domains (a day or more) and pstk-side-trees
# in parallel, then pstk-summaries (minutes).
pstk-trees:
	rm -f results/trees/pstk.aln results/trees/pstk.trimmed results/trees/pstk.treefile results/trees/pstk.iqtree \
	  results/trees/pstk.ckp.gz results/trees/pstk.log results/trees/C20/pstk.*
	python3 stage2/machinery_phylogeny.py --families pstk --max-per-group 3 $(MACHINERY_TREE_INPUTS)
	mkdir -p results/trees/C20 && cp results/trees/pstk.trimmed results/trees/C20/
	python3 stage2/machinery_phylogeny.py --families pstk --max-per-group 3 --model LG+C20+F+G --out results/trees/C20 \
	  --figures results/trees/C20 $(MACHINERY_TREE_INPUTS)
	python3 figures/phylogeny_collapsed.py
	python3 figures/phylogeny_collapsed.py --families pstk SEPHS2 --trees results/trees/C20 --figures results/trees/C20

pstk-lba-domains:
	rm -f results/trees/lba/pstk.* results/trees/domains/pstk.*
	python3 stage2/machinery_lba.py
	python3 stage2/machinery_domains.py
	python3 stage2/machinery_domains.py --pmsf --table results/tables/stage2_machinery_domains_pmsf.tsv \
	  --conflicts results/tables/stage2_machinery_domain_conflicts_pmsf.tsv

pstk-side-trees: pstk-no-outgroup pstk-enriched

pstk-summaries:
	$(MAKE) machinery-lba-figures machinery-domain-figures
	$(MAKE) machinery-domain-figures DOMAIN_SUFFIX=.pmsf
	$(MAKE) machinery-clade-supports pstk-rates
	$(MAKE) machinery-table   # after the clade supports, which give its SECIS-binding protein row

# SECIS-binding protein analyses again after a change of its sequence set only (2026-11-09: deuterostome set rooted on
# protostomes, SBP2_TREE_OPTS; 2026-11-10: truncated proteins removed, full-length ones added). The SBP2 outputs of each analysis are removed; the scripts keep the finished analyses
# of PSTK and SEPHS2. On Deigo: sbp2-trees (an hour or two), then sbp2-lba-domains (a day or more), then
# sbp2-summaries (minutes).
sbp2-trees:
	rm -f results/trees/SBP2.aln results/trees/SBP2.trimmed results/trees/SBP2.treefile results/trees/SBP2.iqtree \
	  results/trees/SBP2.ckp.gz results/trees/SBP2.log
	python3 stage2/machinery_phylogeny.py --families SBP2 $(SBP2_TREE_OPTS) $(MACHINERY_TREE_INPUTS)
	python3 figures/phylogeny_collapsed.py --families SBP2

sbp2-lba-domains:
	rm -f results/trees/lba/SBP2.* results/trees/domains/SBP2.*
	python3 stage2/machinery_lba.py
	python3 stage2/machinery_domains.py
	python3 stage2/machinery_domains.py --pmsf --table results/tables/stage2_machinery_domains_pmsf.tsv \
	  --conflicts results/tables/stage2_machinery_domain_conflicts_pmsf.tsv

sbp2-summaries:
	$(MAKE) machinery-lba-figures machinery-domain-figures
	$(MAKE) machinery-domain-figures DOMAIN_SUFFIX=.pmsf
	$(MAKE) machinery-clade-supports
	$(MAKE) machinery-table   # after the clade supports, which give its SECIS-binding protein row

# SECIS-binding protein figures again (seconds), after a change of the drawing code only (2026-11-10: rooting on the
# outgroup as a clade; best human hit dropped from the labels of invertebrate proteins): full and collapsed
# maximum-likelihood tree, PMSF collapsed tree and the domain trees, SECIS-binding proteins only
sbp2-figures:
	python3 stage2/machinery_phylogeny.py --families SBP2 --draw-only
	python3 figures/phylogeny_collapsed.py --families SBP2
	python3 figures/phylogeny_collapsed.py --families SBP2 --trees results/trees/lba --tree-suffix .pmsf --suffix _pmsf \
	  --model-label "LG+C60+F+G (PMSF)"
	$(MAKE) machinery-domain-figures DOMAIN_TREES="SBP2:upstream_of_L7Ae SBP2:L7Ae"
	$(MAKE) machinery-domain-figures DOMAIN_TREES="SBP2:upstream_of_L7Ae SBP2:L7Ae" DOMAIN_SUFFIX=.pmsf

pstk-no-outgroup:
	rm -rf results/trees/lba_no_outgroup
	python3 stage2/machinery_lba.py --families pstk --analyses taxa --exclude KTI12 --out results/trees/lba_no_outgroup \
	  --table results/tables/stage2_machinery_lba_no_outgroup.tsv

machinery-lba-figures:
	python3 figures/phylogeny_collapsed.py --trees results/trees/lba --tree-suffix .pmsf --suffix _pmsf --model-label "LG+C60+F+G (PMSF)"

# Domain analysis of the machinery trees (after machinery-phylogeny): one tree per domain (domains defined on human
# proteins, curation/machinery_domains.tsv), supports and AU tests of each domain's tree against the other domain's data
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=1-00:00:00 scripts/run_in_container.sh make machinery-domains
machinery-domains:
	rm -rf results/trees/domains
	python3 stage2/machinery_domains.py
	$(MAKE) machinery-domain-figures

# The same under LG+C60+F+G with PMSF (after machinery-domains, whose maximum-likelihood trees are the guide trees, and
# after machinery-lba, whose full-length PMSF trees it compares; keeps finished analyses, so it can be resubmitted after
# an interruption)
#   sbatch --dependency=afterok:<machinery-domains job>:<machinery-lba job> ...
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

# Machinery genes of the three O. dioica assemblies and of O. longicauda, O. albicans and O. vanhoeffeni (decisions of
# 2026-11-01): anchors from the machinery predictions of the final stage 2 run and human reference proteins
# (curation/machinery_anchors.py); a first pass on the O. dioica assemblies (curation/gold_standard.py --anchor sense,
# candidates ranked by similarity to the human proteins; the transcript alignments of gold-standard are reused); miniprot
# alignments of the human proteins and of the first-pass proteins on the six genomes (curation/machinery_projections.py);
# the final pass on the six genomes, with these alignments as candidate models and the first-pass proteins of the other
# genomes added to the references; exon chains of curation/machinery_overrides.tsv set by the curator, models of
# curation/machinery_incomplete.tsv marked incomplete (minutes)
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=2:00:00 scripts/run_in_container.sh make gold-standard-machinery
GOLD_MACHINERY = curation/gold_standard_machinery
MACHINERY_GOLD_GENOMES = OKI OSKA OdB3 Olon Oalb Ovan
MACHINERY_WORK = work/stage2/gold_standard_machinery
gold-standard-machinery:
	python3 curation/machinery_anchors.py --gtf-dir $(SELENO_PREDICTIONS) --genomes $(MACHINERY_GOLD_GENOMES) \
	  --human data/genomes/Hsap --references curation/machinery_references.faa --extra curation/machinery_extra_anchors.tsv \
	  --out curation/machinery_genes.tsv
	python3 curation/gold_standard.py --anchor sense --genomes OKI OSKA OdB3 --baseline curation/machinery_genes.tsv \
	  --stage1 $(SELENO_PREDICTIONS) --references curation/machinery_references.faa --overrides curation/machinery_overrides.tsv \
	  --out $(MACHINERY_WORK)/first_pass --table $(MACHINERY_WORK)/first_pass.tsv > /dev/null
	python3 curation/machinery_projections.py --queries curation/machinery_references.faa $(MACHINERY_WORK)/first_pass/*.proteins.faa \
	  --genomes $(MACHINERY_GOLD_GENOMES) --out $(MACHINERY_WORK)/projections
	python3 curation/gold_standard.py --anchor sense --genomes $(MACHINERY_GOLD_GENOMES) --baseline curation/machinery_genes.tsv \
	  --stage1 $(SELENO_PREDICTIONS) --references curation/machinery_references.faa \
	  --reference-proteins $(MACHINERY_WORK)/first_pass --projections $(MACHINERY_WORK)/projections \
	  --overrides curation/machinery_overrides.tsv --incomplete curation/machinery_incomplete.tsv \
	  --other-proteins $(wildcard results/gene_regions/*/proteins.faa) \
	  --out $(GOLD_MACHINERY) --table results/tables/gold_standard_machinery.tsv

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

# Step 5: candidate genes outside the baseline (curation/candidates.tsv and .faa: TXNL1.b, the KKSU family; SELENOI
# from the absence tests of step 9, its codons located through its transcripts aligned to the O. dioica assemblies,
# stage2/candidate_transcripts.py, which reads the transcript alignments of make gold-standard): the rule
# fixed in advance in the three O. dioica assemblies, tblastn against the appendicularian, tunicate and amphioxus
# genomes (codon at the selenocysteine position, SECIS after the stop codon) and blastp against RefSeq proteomes
# (about an hour; minutes when the genome searches are cached), for example
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
	python3 stage2/candidate_transcripts.py
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

# Appendicularian predictions of the selenoprotein trees. Round 2 (step 8, from 2026-10-29): stage 2 final run, with
# the O. dioica genes as gold-standard proteins; first pass (2026-10-24): make selenoprotein-rbh selenoprotein-phylogeny
# selenoprotein-orthology SELENO_PREDICTIONS=work/stage1/pipeline/selenoprofiles GOLD=
SELENO_PREDICTIONS = work/stage2/pipeline/final/selenoprofiles
GOLD = curation/gold_standard
# Inputs of the machinery trees (machinery-phylogeny; decision of 2026-11-01)
MACHINERY_TREE_INPUTS = --stage2 $(SELENO_PREDICTIONS) --gold $(GOLD) $(GOLD_MACHINERY) --gold-genomes $(MACHINERY_GOLD_GENOMES) \
  --exclude-genomes Bsty Mery
# Best hits: SELENO_PREDICTIONS for the genomes they hold, stage 1 predictions for the other genomes
SELENO_RBH = $(sort $(wildcard $(SELENO_PREDICTIONS)/*.selenoprofiles.faa) \
  $(filter-out $(addprefix work/stage1/pipeline/selenoprofiles/,$(notdir $(wildcard $(SELENO_PREDICTIONS)/*.selenoprofiles.faa))), \
    $(wildcard work/stage1/pipeline/selenoprofiles/*.selenoprofiles.faa)))

# Reciprocal best hits of the selenoprotein predictions (all genomes) against the human RefSeq proteome
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=12:00:00 scripts/run_in_container.sh make selenoprotein-rbh
selenoprotein-rbh:
	python3 stage2/machinery_rbh.py --families selenoproteins --predictions $(SELENO_RBH) \
	  --table results/tables/stage2_selenoprotein_rbh.tsv

# Trees of the families of SELENO_TREES that have appendicularian sequences (SELENO_PREDICTIONS, GOLD); human RefSeq
# proteins of each family (and MIEN1 for SELENOW) as references; starts afresh, removing results/trees/selenoproteins
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=2-00:00:00 scripts/run_in_container.sh make selenoprotein-phylogeny
selenoprotein-phylogeny:
	rm -rf results/trees/selenoproteins
	python3 stage2/machinery_phylogeny.py --families $(SELENO_TREES) --max-per-group 3 --min-fraction 0.3 \
	  --stage2 $(SELENO_PREDICTIONS) $(if $(GOLD),--gold $(GOLD)) \
	  --out results/trees/selenoproteins --figures figures/selenoproteins --formats pdf
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

# AU tests of the positions of stage2/au_hypotheses.tsv (step 8), under the ML model and LG+C60+F+G (PMSF), after
# selenoprotein-pmsf (keeps finished analyses); clades resolved in results/trees/selenoproteins/lba/<family>.hypotheses.tsv
#   sbatch -p compute --cpus-per-task=8 --mem=32G --time=1-00:00:00 scripts/run_in_container.sh make selenoprotein-au
SELENO_AU = DIO SELENOF SELENOK SELENOO SELENOS SELENOT SELENOW
selenoprotein-au:
	python3 stage2/machinery_lba.py --families $(SELENO_AU) --trees results/trees/selenoproteins --analyses pmsf hypotheses \
	  --out results/trees/selenoproteins/lba --table results/tables/stage2_selenoprotein_lba.tsv

# First-pass orthology of each appendicularian clade, from the trees and best hits (seconds; reads results/ only)
selenoprotein-orthology:
	python3 stage2/selenoprotein_orthology.py $(if $(GOLD),,--annotation results/tables/stage1_annotation.tsv)

# Absence tests of the selenoprotein families found in other tunicates but in no appendicularian (step 9; census):
# the selenocysteine forms predicted in stage 1 (other tunicates, cephalochordates, vertebrates; and F. borealis for
# SELENOL) searched in the appendicularian genomes and transcriptomes (tblastn), each candidate compared with them and
# with human proteins of related families (outgroups). MSRB: MSRB1-type references only (selenocysteine forms, no
# profile), against MSRB2 and MSRB3, since every appendicularian has cysteine MSRB of the MSRB2/MSRB3 type.
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=1-00:00:00 scripts/run_in_container.sh make selenoprotein-absence
SELENO_ABSENT = SELENON SELENOH MSRB FrnE SELENOL SELENOI
selenoprotein-absence:
	python3 stage2/machinery_absence.py --gene SELENON --labels selenocysteine --outgroup-human CALM1 --out work/stage2/selenoprotein_absence
	python3 stage2/machinery_absence.py --gene SELENOH --labels selenocysteine --outgroup-human SELENOW,SELENOV,SELENOT --out work/stage2/selenoprotein_absence
	python3 stage2/machinery_absence.py --gene MSRB --labels selenocysteine --no-profile --outgroup-human MSRB2,MSRB3 --out work/stage2/selenoprotein_absence
	python3 stage2/machinery_absence.py --gene FrnE --labels selenocysteine --outgroup-human GSTK1 --out work/stage2/selenoprotein_absence
	python3 stage2/machinery_absence.py --gene SELENOL --labels selenocysteine --outgroup-human TXN --out work/stage2/selenoprotein_absence
	python3 stage2/machinery_absence.py --gene SELENOI --labels selenocysteine --outgroup-human CEPT1,CHPT1 --out work/stage2/selenoprotein_absence
	for g in $(SELENO_ABSENT); do for f in summary candidates controls; do \
	  cp work/stage2/selenoprotein_absence/$$g.$$f.tsv results/tables/stage2_absence_$$g.$$f.tsv; done; done

# Table SELENO-ORTHOLOGY of the manuscript (round 2: gold-standard genes), after selenoprotein-orthology (seconds)
selenoprotein-table:
	python3 stage2/selenoprotein_table.py

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

# Genomic regions of PSTK outside the gold standard, to read the splice sites of its introns by hand (seconds):
# gold-standard assemblies and the O. albicans and O. vanhoeffeni models, the assemblies of other individuals of
# the three cryptic species (region aligned from the gold standard), the C. intestinalis RefSeq gene, and the
# fragments of O. longicauda, M. erythrocephalus and B. stygius (loci of results/tables/stage2_absence_pstk.candidates.tsv)
pstk-regions:
	python3 curation/gene_regions.py --gene PSTK --gold OKI OSKA OdB3 Oalb Ovan --lift KUM=OKI AOM=OSKA Bar2=OdB3 \
	  --annotation Cint=XP_002120985.2 Crob=XP_002120985.1 --contig Olon=SCLD01107789.1:4496-5011 Olon=SCLD01051426.1:5-604 \
	  Mery=SCLF01145532.1:2326-2851 Bsty=SCLE01192721.1:1101-1292 --flank 3000 --out results/gene_regions/PSTK

# Regions of the other machinery genes under review (2026-11-03), as for PSTK: gold-standard models with 3 kb on each
# side, the same regions in the assemblies of other individuals (minimap2), and miniprot alignments of an O. dioica
# protein to the other genomes (second SEPHS copy in the assemblies without a model; independent exon chains in the
# other species); results/gene_regions/<gene>/ (FASTA, regions.tsv, <genome>.miniprot.gff3); miniprot indexes in
# work/gene_regions. About 30 min on 4 CPUs:
#   sbatch -p compute --cpus-per-task=4 --mem=16G --time=2:00:00 scripts/run_in_container.sh make machinery-regions
machinery-regions:
	python3 curation/gene_regions.py --gene SBP2 --gold OKI OSKA OdB3 Olon Oalb Ovan --lift KUM=OKI AOM=OSKA Bar2=OdB3 \
	  --protein Olon=OdB3 Oalb=OdB3 Ovan=OdB3 --flank 3000 --out results/gene_regions/SBP2
	python3 curation/gene_regions.py --gene SBP2_2 --gold Olon --flank 3000 --out results/gene_regions/SBP2_2
	python3 curation/gene_regions.py --gene SEPHS_b --gold OKI Olon Oalb Ovan --lift KUM=OKI \
	  --protein KUM=OKI OSKA=OKI AOM=OKI OdB3=OKI Bar2=OKI Olon=OKI Oalb=OKI Ovan=OKI --flank 3000 --out results/gene_regions/SEPHS_b
	python3 curation/gene_regions.py --gene SEPSECS --gold OKI OSKA OdB3 Olon Oalb Ovan --lift KUM=OKI AOM=OSKA Bar2=OdB3 \
	  --protein Olon=OKI Oalb=OKI Ovan=OKI --flank 3000 --out results/gene_regions/SEPSECS
	python3 curation/gene_regions.py --gene EEFSEC --gold OKI OSKA OdB3 Olon Oalb Ovan --lift KUM=OKI AOM=OSKA Bar2=OdB3 \
	  --protein Olon=OKI Oalb=OKI Ovan=OKI --flank 3000 --out results/gene_regions/EEFSEC

# Proteins of the machinery genes in the assemblies of other individuals and in the alignments to assemblies without
# a model (exonerate on the regions of make pstk-regions and machinery-regions); local, minutes. They serve
# gold-standard-machinery (--other-proteins) to tell assembly errors from transcript errors.
region-proteins:
	python3 curation/region_proteins.py --genes PSTK SBP2 SBP2_2 SEPHS_b SEPSECS EEFSEC

# Intron positions of PSTK across chordates, on the coordinates of human PSTK (seconds): RefSeq models of
# vertebrates, amphioxus and tunicates (proteins of results/trees/pstk.faa; the Thalia democratica transcript found by
# tblastn of the RefSeq transcripts with the C. intestinalis protein) and the gold-standard models of the Oikopleura
# species; then curation/intron_columns.py (local, MAFFT) places the introns on an alignment
pstk-introns:
	python3 curation/intron_positions.py --gene PSTK --reference Hsap=NP_001350460.1 \
	  --models Mmus=NP_001034623.1 Drer=XP_009296681.1 Drer=XP_009305284.1 Blan=XP_066278120.1 Cint=XP_002120985.2 \
	  Crob=XP_002120985.1 Pmam=XP_079949190.1 Scla=XP_039273489.2 Bsch=XP_080181757.1 Tdem=XM_082946145.1 \
	  --gold OKI OSKA OdB3 Oalb Ovan --out results/tables/pstk_intron_positions.tsv

# How fast PSTK evolves (local, seconds): distances of PSTK and of KTI12 between the same chordates in the PSTK trees
# (whole alignment and its two parts, selected model and PMSF), and identities of the human, Ciona and O. dioica
# proteins of the genes with references in figures/inputs
pstk-rates:
	python3 stage2/pstk_rates.py

# Intron positions of SEPSECS and PSTK in the appendicularian models on one protein alignment (local, seconds; needs
# MAFFT, exonerate and the regions of make machinery-regions): gold-standard models, and the assemblies of other
# individuals by exonerate of the gold-standard protein of the same species
INTRON_MAP_GENES = SEPSECS PSTK EEFSEC
intron-map:
	for g in $(INTRON_MAP_GENES); do python3 curation/intron_map.py --gene $$g; done

# Junctions of the curator's EEFSEC chains of the species without transcripts, rescanned against the OSKA2016 and OdB3
# proteins (seconds; evidence of curation/machinery_overrides.tsv; needs make machinery-regions)
# Presence of each O. dioica gold-standard intron in the assembly of another individual of the same species
# (KUM-M3-7f, AOM-5-5f, Bar2_p4): introns that are insertions absent from the other assembly (minutes; genomes)
intron-presence:
	python3 curation/intron_presence.py --introns results/tables/intron_content.tsv --out results/tables/intron_presence.tsv

# Paired intron ends and splice sites in the O. dioica gold-standard introns, with and without the insertions of
# make intron-presence (seconds)
intron-pairing:
	python3 curation/intron_pairing.py --out results/tables/intron_pairing.tsv

# Coding divergence between the two individuals of each O. dioica species, codon by codon (seconds; needs make
# machinery-regions)
cds-divergence:
	python3 curation/cds_divergence.py --genes EEFSEC PSTK SEPSECS SEPHS_b SBP2 --out results/tables/cds_divergence.tsv

junction-rescan:
	python3 curation/junction_rescan.py --gene EEFSEC --assemblies Olon Oalb Ovan --out results/tables/junction_rescan_EEFSEC.tsv

# Intron-exon structures of PSTK in the species of the manuscript on a protein alignment (seconds; needs make
# pstk-introns and pstk-regions): figures/pstk_structures.{pdf,png,svg}, results/tables/pstk_structures.tsv
pstk-structures:
	python3 figures/pstk_structures.py --out figures/pstk_structures --table results/tables/pstk_structures.tsv

# PSTK introns as gene structures on the residues of human PSTK (local, seconds; MAFFT, exonerate, legacy BLAST; needs
# pstk-introns and pstk-regions): one row per group of identical structures, and for O. longicauda and the giant
# larvaceans the tblastn HSPs of the O. vanhoeffeni protein with the best class-site intron at each junction:
# figures/pstk_intron_scheme.{pdf,png,svg}, results/tables/pstk_intron_scheme.tsv
pstk-intron-scheme:
	python3 figures/pstk_intron_scheme.py --out figures/pstk_intron_scheme --table results/tables/pstk_intron_scheme.tsv

pstk-intron-columns:
	python3 curation/intron_columns.py --positions results/tables/pstk_intron_positions.tsv --out results/tables/pstk_intron_columns.tsv \
	  --alignment results/trees/pstk_introns.aln

# Content of the introns of the gold-standard models (OKI2018, OSKA2016, OdB3; machinery models of O. albicans):
# secondary structure and pairing of the two ends (ViennaRNA, against dinucleotide shuffles) and copies elsewhere in
# the genome (legacy BLASTN), with the PSTK introns marked (minutes to an hour; formats BLAST databases in work/)
intron-content:
	python3 curation/intron_content.py --genomes OKI OSKA OdB3 Oalb --focus PSTK --out results/tables/intron_content.tsv

# The 9-kb intron of OdB3 EEFSEC as an insertion: empty site in Bar2_p4 (results/gene_regions/EEFSEC, make
# machinery-regions), target-site duplication, terminal repeats, reading frames, copies in the six O. dioica
# assemblies (BLAST databases shared with make intron-content)
eefsec-insertion:
	python3 curation/insertion.py --genome OdB3 --intron FN653026.1:876889-885870:+ --empty results/gene_regions/EEFSEC/Bar2.fa \
	  --copies-in OdB3 Bar2 OKI KUM OSKA AOM \
	  --out results/tables/eefsec_insertion.tsv --copies results/tables/eefsec_insertion_copies.tsv

# Reannotation panels of the machinery genes (one gene per figure, aligned on the start codon of the gold-standard
# model, no SECIS), from the machinery gold standard and figures/inputs (figures/machinery_reannotation_inputs.sh,
# which needs www.ebi.ac.uk and the genomes); panel C, equivalent placements of introns, when
# figures/inputs/intron_placements_<gene>.tsv exists (figures/intron_placements.py, from results/gene_regions)
MACHINERY_FIGURE_GENES = PSTK SBP2 SEPHS_b SEPSECS EEFSEC
# inputs of the machinery panels (figures/machinery_reannotation_inputs.sh): needs the EBI APIs (ENA, InterPro,
# dbfetch) and the OKI2018, OSKA2016 and OdB3 genomes; minutes. Then make machinery-reannotation (local, seconds)
#   scripts/run_in_container.sh make machinery-reannotation-inputs
machinery-reannotation-inputs:
	OKI_GENOME=data/genomes/OKI/genome.fna.gz OSKA_GENOME=data/genomes/OSKA/genome.fna.gz ODB3_GENOME=data/genomes/OdB3/genome.fna.gz \
	  sh figures/machinery_reannotation_inputs.sh

# Spliced leader at the 5' end of the transcripts drawn in the machinery panels (seconds; reads the transcript files of
# stage1/samplesheet.csv): results/tables/sl_transcripts.tsv, read by make machinery-reannotation ("SL" marks)
#   scripts/run_in_container.sh make sl-transcripts
sl-transcripts:
	python3 curation/sl_transcripts.py --genes $(MACHINERY_FIGURE_GENES) --out results/tables/sl_transcripts.tsv

machinery-reannotation:
	for g in $(MACHINERY_FIGURE_GENES); do python3 figures/reannotation_page.py --genes $$g --columns 1 --origin start \
	  --no-secis --gold $(GOLD_MACHINERY) --fit --name machinery_$$g --sources "human: UniProt; Ciona: RefSeq" \
	  $$([ -f figures/inputs/intron_placements_$$g.tsv ] && echo --side-panel figures/inputs/intron_placements_$$g.tsv \
	  --side-names KUM=KUM-M3-7f AOM=AOM-5-5f Bar2=Bar2_p4) \
	  $$([ $$g = PSTK ] && echo --tree-panel results/trees/pstk_enriched_archaea/pstk.pmsf.treefile --tree-groups results/trees/pstk_enriched_archaea/pstk.groups.tsv); done

# Table MACHINERY of the manuscript: census, reciprocal best hits, trees and AU tests (seconds; reads results/ only)
machinery-table:
	python3 stage2/machinery_table.py

# supports of the other clades cited in the text, from all machinery trees (full length, C20, long-branch analyses,
# one leaf per species, domains); local, seconds; after machinery-lba, machinery-domains(-pmsf) and pstk-one-leaf
machinery-clade-supports:
	python3 stage2/machinery_clade_supports.py

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

# Which filter of a profile removes a prediction (Selenoprofiles with one profile per job; about 1 h on 4 CPUs):
#   sbatch [--partition=<partition>] --cpus-per-task=4 scripts/run_in_container.sh make stage2-diagnose
stage2-diagnose:
	python3 stage2/diagnose_refiltering.py

# SELENOK loci of the UGA/SECIS screens in the other Oikopleura species (and M. erythrocephalus), which the final
# SELENOK profile does not find (make selenok-diagnose: no candidate prediction in O. albicans, O. vanhoeffeni and
# M. erythrocephalus; an "unaligned" fragment in O. longicauda): regions around them, to build their gene models
# (seconds to minutes), for example
#   sbatch -p compute --cpus-per-task=1 --mem=8G --time=1:00:00 scripts/run_in_container.sh make selenok-regions
selenok-regions:
	@mkdir -p results/sequences
	python3 stage2/extract_regions.py --loci curation/selenok_loci.tsv --out results/sequences/selenok_regions.fa

# Gene models of these SELENOK genes: genewise of the O. dioica gold-standard SELENOK proteins, joined to the exon of
# the screen codon or extended through TGA (reads results/sequences/selenok_regions.fa; seconds)
selenok-models:
	python3 stage2/region_gene_models.py --regions results/sequences/selenok_regions.fa --family SELENOK \
	  --queries curation/gold_standard/OKI.proteins.faa curation/gold_standard/OSKA.proteins.faa curation/gold_standard/OdB3.proteins.faa \
	  --faa results/sequences/selenok_models.faa --table results/tables/stage2_selenok_models.tsv

# Release profiles for the deposit: the final profiles plus the SELENOK models of O. longicauda, O. albicans and
# O. vanhoeffeni (make selenok-models; the M. erythrocephalus model is left out, held-out genome), built as the final
# ones (Selenoprofiles in the image; minutes); test with make selenok-diagnose SELENOK_DIAG_BUILDS="final release"
stage2-release-profiles:
	python3 stage2/build_profiles.py --stage1-genomes OKI OSKA OdB3 Olon Oalb Ovan --validated-genomes Olon Oalb Ovan \
	    --candidates curation/candidates.faa --curated curation/machinery_proteins.faa results/sequences/selenok_models.faa \
	    --out $(STAGE2_PROFILES)/release --archive work/stage2/appendicularia_selenoprofiles_release.tar.gz
	awk 'NR == 1 {print "build\t" $$0} NR > 1 {print "release\t" $$0}' $(STAGE2_PROFILES)/release/members.tsv \
	    > results/tables/stage2_profile_members_release.tsv

# SELENOI regions for a wider SECIS search (genome regions around the gene in the O. dioica assemblies, the other
# appendicularians, tunicates, amphioxus and vertebrates, and the O. dioica transcripts; a few minutes, reads the
# genomes), for example
#   sbatch -p compute --cpus-per-task=1 --mem=8G --time=1:00:00 scripts/run_in_container.sh make selenoi-regions
selenoi-regions:
	python3 stage2/selenoi_regions.py

# Wider search for the SECIS element of SELENOI, from stem-loops rather than from a SECIS model, with covariance
# models without filters, core patterns, a comparison of five Oikopleura species and calibrations on the gold-standard
# 3' sequences (reads results/sequences/selenoi_regions.fa and repository files; about 15 minutes on 4 CPUs), for
#   sbatch -p compute --cpus-per-task=4 --mem=8G --time=2:00:00 scripts/run_in_container.sh make selenoi-secis
selenoi-secis:
	python3 stage2/selenoi_secis_search.py

# Step 9 follow-up 2: why the final SELENOK profile misses the SELENOK genes that the UGA/SECIS screens found in
# O. longicauda, O. albicans, O. vanhoeffeni and M. erythrocephalus (and how it does in the other appendicularians and
# OKI2018): candidate predictions before filtering, with their proteins (minutes), for example
#   sbatch -p compute --cpus-per-task=4 --mem=16G --time=4:00:00 scripts/run_in_container.sh make selenok-diagnose
SELENOK_DIAG_GENOMES = OKI Olon Oalb Ovan Bsty Mery Fbor
SELENOK_DIAG_BUILDS = final
selenok-diagnose:
	python3 stage2/diagnose_refiltering.py --jobs $(foreach b,$(SELENOK_DIAG_BUILDS),$(foreach g,$(SELENOK_DIAG_GENOMES),$(g):SELENOK:$(b))) \
	  --out work/stage2/diagnose_selenok --table results/tables/stage2_selenok_diagnosis.tsv \
	  --faa results/sequences/stage2_selenok_diagnosis.faa

# Step 10: selenoprotein genes, operons and SECIS elements in the three O. dioica assemblies (public annotations,
# gold-standard models, transcript alignments of make gold-standard; minutes), for example
#   sbatch -p compute --cpus-per-task=4 --mem=16G --time=2:00:00 scripts/run_in_container.sh make operons
operons:
	python3 stage2/operons.py

# Figure OPERONS from the tables of make operons (seconds; matplotlib only)
operons-figure:
	python3 figures/operons_figure.py

# Genes sharing an operon with a selenoprotein gene (stage2/operon_partners.py): best human hits and redox/selenium
# flags, against the other genes in operons. Three steps:
#   1. make operon-proteins (container): proteins of all annotated loci of the three assemblies
#   2. make operon-search, outside the container, with the Swiss-Prot DIAMOND database of the OIST Bioinformatics
#      User Group (https://github.com/oist/BioinfoUgrp/blob/master/DataBases.md; scripts/operon_search.sh loads
#      Other/DIAMOND/2.0.15.153 and DB/diamondDB/sprot/2022-07), human entries only
#   3. make operon-partners (container): tables
# Optional, before step 3: make human-uniprot (a node with internet access, at OIST the login node; about 2 MB) for the
# UniProt keywords and GO terms of the human proteins; without it, the flags come from the protein names only.
OPERON_WORK = work/stage2/operon_partners
HUMAN_UNIPROT_QUERY = %28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29
human-uniprot:
	mkdir -p data/references
	curl -sSfL --retry 4 --retry-delay 20 -D data/references/uniprot_human_sp.headers.txt -o data/references/uniprot_human_sp.tsv \
	  "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession%2Cid%2Cgene_primary%2Cprotein_name%2Ckeyword%2Cgo_id&query=$(HUMAN_UNIPROT_QUERY)"
	gzip -f data/references/uniprot_human_sp.tsv
	echo "$$(zcat data/references/uniprot_human_sp.tsv.gz | tail -n +2 | wc -l) entries"
	md5sum data/references/uniprot_human_sp.tsv.gz | tee data/references/uniprot_human_sp.md5
	date -u +%Y-%m-%dT%H:%MZ > data/references/uniprot_human_sp.date
	grep -i "x-uniprot-release" data/references/uniprot_human_sp.headers.txt || true

operon-proteins:
	python3 stage2/operon_partners.py --proteins-only

operon-search:
	@test -n "$${DIAMONDDB:-}" || { echo "DIAMONDDB is not set: load a DB/diamondDB/sprot module (scripts/operon_search.sh)" >&2; exit 1; }
	ls -l $${DIAMONDDB}
	db=$$(ls $${DIAMONDDB}/*.dmnd | head -n 1); echo "database: $$db"; diamond --version; \
	for g in OKI OSKA OdB3; do \
	  q=$(OPERON_WORK)/$$g.proteins.faa; o=$(OPERON_WORK)/$$g.diamond.tsv; \
	  opts="--more-sensitive --evalue 1e-5 --threads $${SLURM_CPUS_PER_TASK:-1} --outfmt 6 qseqid sseqid pident evalue bitscore stitle"; \
	  if diamond blastp --db $$db --query $$q --taxonlist 9606 --max-target-seqs 5 $$opts --out $$o.part; then \
	    mv $$o.part $$o; \
	  else \
	    echo "--taxonlist failed: all of Swiss-Prot, keeping the _HUMAN entries" >&2; \
	    diamond blastp --db $$db --query $$q --max-target-seqs 200 $$opts --out $$o.all \
	    && awk -F'\t' 'index($$2, "_HUMAN") || index($$6, "_HUMAN")' $$o.all > $$o || exit 1; \
	  fi; \
	done

operon-partners:
	python3 stage2/operon_partners.py

# miniprot on the gold-standard proteins (selenocysteine as U, and as C), compared with the gold-standard junctions near
# the selenocysteine codon; why miniprot misses SELENOT and SEPHS2 (minutes), for example
#   sbatch -p compute --cpus-per-task=8 --mem=16G --time=2:00:00 scripts/run_in_container.sh make miniprot-junctions
miniprot-junctions:
	python3 stage2/miniprot_junctions.py

# The same comparison on the alignments of the stage 2 final pipeline run itself (miniprot 0.11, profile queries;
# seconds); with miniprot 0.11 on the gold-standard proteins: sbatch -p compute scripts/miniprot_v011.sh
PIPELINE_FINAL ?= work/stage2/pipeline/final
# table from the miniprot 0.11 alignments of scripts/miniprot_v011.sh (rebuild without realigning)
miniprot-junctions-v011-table:
	python3 stage2/miniprot_junctions.py --from-gff work/stage2/miniprot_junctions_v0.11 \
	  --work work/stage2/miniprot_junctions_v0.11 --table results/tables/stage2_miniprot_junctions_v0.11.tsv
miniprot-junctions-pipeline:
	@test -d $(PIPELINE_FINAL)/miniprot || { echo "no $(PIPELINE_FINAL)/miniprot: set PIPELINE_FINAL to the outdir of the stage 2 final run" >&2; exit 1; }
	@q=$$(ls $(PIPELINE_FINAL)/summary/selenoprotein_queries.faa $(PIPELINE_FINAL)/profile/*.faa 2>/dev/null | head -n 1); \
	test -n "$$q" || { echo "no query FASTA (selenoprotein_queries.faa) in $(PIPELINE_FINAL)/summary" >&2; exit 1; }; \
	echo "queries: $$q"; ls $(PIPELINE_FINAL)/miniprot; \
	python3 stage2/miniprot_junctions.py --from-gff $(PIPELINE_FINAL)/miniprot --queries $$q \
	  --work work/stage2/miniprot_junctions_pipeline --table results/tables/stage2_miniprot_junctions_pipeline.tsv

# Peptides for a mass-spectrometry validation (planning, in silico digestion of the gold-standard proteins; seconds)
ms-peptides:
	python3 stage2/ms_peptides.py

# SELENOI: placement and splice sites of the junction that removes the TGA segment (minimap2 on the sequences of
# make selenoi-regions, which are in the repository; seconds)
selenoi-junction:
	python3 stage2/selenoi_junction.py

# Step 9: census of the appendicularian selenoproteins from the stage 2 final run (stage2/README.md, "Census").
# 1. stage2-collect (minutes): tables of the final run (results/tables/stage2_final_*.tsv) and its machinery census
# 2. census-contigs (about an hour): composition of the contigs with predictions (Tiara classes if a Tiara output
#    is placed in work/stage2/census/tiara/; not run, decision of 2026-10-30 in docs/plan.md)
# 3. census (seconds; reads results/ only, after selenoprotein-rbh, -phylogeny and -pmsf of round 2)
CENSUS_GENOMES = OKI OSKA OdB3 KUM AOM Bar2 Olon Oalb Ovan Bsty Mery Fbor
stage2-collect:
	python3 stage1/collect.py --outdir work/stage2/pipeline/final --out-dir results/tables --prefix stage2_final
	mkdir -p work/stage2/census/machinery
	for g in $(CENSUS_GENOMES); do ln -sf ../../pipeline/final/selenoprofiles/$$g.selenoprofiles.gtf work/stage2/census/machinery/$$g.gtf; done
	python3 stage2/machinery_census.py --trnasec results/tables/stage2_final_trnasec.tsv --test-dir work/stage2/census/machinery \
	  --fragments work/stage2/machinery_fragments/fragments.tsv \
	  --table results/tables/stage2_final_machinery_census.tsv --markdown results/tables/stage2_final_machinery_census.md

#   sbatch -p compute --cpus-per-task=2 --mem=32G --time=12:00:00 scripts/run_in_container.sh make census-contigs
census-contigs:
	python3 stage2/census_contigs.py --genomes $(CENSUS_GENOMES)

census:
	python3 stage2/census.py

# Deposit archive of the recommended settings and data files (stage2/recommended/README.md), after stage2-inputs
STAGE2_BUNDLE = work/stage2/appendicularia_selenoprotein_annotation
# The bundle carries the release profiles (final profiles plus SELENOK of three Oikopleura species; make
# stage2-release-profiles) under the name the recommended settings use, and the final profiles of the census
# (appendicularia_selenoprofiles_final.tar.gz) for reproducibility.
stage2-bundle:
	@test -s work/stage2/appendicularia_selenoprofiles.tar.gz || { echo "run make stage2-inputs first"; exit 1; }
	@test -s work/stage2/appendicularia_selenoprofiles_release.tar.gz || { echo "run make stage2-release-profiles first"; exit 1; }
	rm -rf $(STAGE2_BUNDLE) && mkdir -p $(STAGE2_BUNDLE)
	cp stage2/recommended/README.md stage2/recommended/params.yaml stage2/recommended/samplesheet.csv $(STAGE2_BUNDLE)/
	cp stage2/run/stage2.config $(STAGE2_BUNDLE)/appendicularia.config
	cp stage2/run/expectation_oikopleura_dioica.csv $(STAGE2_BUNDLE)/
	cp work/stage2/appendicularia_selenoprofiles_release.tar.gz $(STAGE2_BUNDLE)/appendicularia_selenoprofiles.tar.gz
	cp work/stage2/appendicularia_selenoprofiles.tar.gz $(STAGE2_BUNDLE)/appendicularia_selenoprofiles_final.tar.gz
	cp work/stage2/models/oikopleura_secis.cm $(STAGE2_BUNDLE)/oikopleura_secis.cm
	cp work/stage2/profiles/release/members.tsv $(STAGE2_BUNDLE)/profile_members.tsv
	cp work/stage2/profiles/final/members.tsv $(STAGE2_BUNDLE)/profile_members_final.tsv
	tar -czf $(STAGE2_BUNDLE).tar.gz -C work/stage2 appendicularia_selenoprotein_annotation
	@echo "$(STAGE2_BUNDLE).tar.gz"

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
