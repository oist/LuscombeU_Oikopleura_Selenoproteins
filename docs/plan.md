# Plan of the analysis

## Stage 1 — published tools and parameters (justifies stage 2)

For the nine appendicularian genomes:
- Selenoprofiles 4.5.8 with the default metazoan profiles; Secmarker for tRNA-Sec.
- Generic protein survey: miniprot alignments of profile sequences, Sec codon = in-frame TGA matching the reference
  Sec context.
- SECIS: Rfam RF00031 (Infernal) on 3′ regions of the Sec genes.
- Public annotations of the three *O. dioica* assemblies compared with the Sec codons (Ticó et al. 2026 classes).
- Transcript evidence: OKI2018 transcriptome, OSKA2016 transcripts, OdB3 ESTs (splice-aware alignment; Sec codon in
  an exon; screens for Sec-extended and C-terminal Sec reading frames with a SECIS).
Outputs that justify stage 2: genes missed by every tool; RF00031 sensitivity on appendicularian 3′ regions;
annotation error classes.

## Curation — baseline from the three *O. dioica* assemblies

Tables in `curation/`: one row per Sec gene and assembly, with the Sec codon, the evidence (tool, transcript, cross-
assembly comparison) and notes on assembly errors (recorded, not fixed).

## Stage 2 — final parameters (delivers the results)

- *Oikopleura* SECIS model built from baseline 3′ regions (Oikopleuridae seeds, OSKA2016 held out as control);
  held-out, shuffled, non-target (other tunicates) and covariation tests; threshold choice.
- Lineage reference proteins (curated *O. dioica* sequences) and profile HMMs for genes missed in stage 1
  (SELENOK, PSTK), each tested with the target genome left out.
- Re-calling on the nine genomes with the final parameters; census per genome; machinery; Cys homologs.
- Operons in the three *O. dioica* assemblies (collinear genes at most 500 nt apart; Plessy et al. 2024).

## Figures

Built from `results/tables/` only.

## Notes for the manuscript

- State that `OKI2018_I68` (NCBI GCF_907165135.1) is the `OKI2018_I69` assembly of Bliznina et al. 2021 (I69: the
  sequenced individual; I68 is a typo in the NCBI record), and that the NCBI name is used for consistency.

## Decisions and open items (2026-10-16)

- Order of the results: stage 1; Sec machinery; complete pool of Sec genes (annotated from proteins and transcripts);
  SECIS model trained on that pool (with same-locus and 90% identity redundancy filters, family balance, and
  leave-one-family-out evaluation).
- tRNA-Sec: standard tools are enough (Secmarker and RF01852 find at least one supported locus in every genome); no
  special-casing.
- Multiple copies of a gene (tRNA-Sec loci, eEFSec, SEPHS2) in some appendicularian assemblies may come from
  incomplete haplotype purging or from contamination: the B. stygius assembly may be contaminated by M.
  erythrocephalus (Pichon et al. 2019, F1000Research, doi:10.12688/f1000research.21551.2; cite where copy numbers or
  B. stygius results are discussed).
- Fragmented assemblies: the pipeline may miss genes whose locus is broken by the assembly (current long-read
  assemblies do not have this problem); the analysis shows such genes as fragments (for example PSTK in
  O. longicauda, B. stygius and M. erythrocephalus) to validate the baseline. At the end, review the assemblies and
  drop those of too low quality that add little phylogenetically.
- Open: which genes form the machinery set (minimum: tRNA-Sec, SEPHS2, PSTK, SecS, eEFSec, SBP2; possible: SECp43,
  SCLY).
- SECIS model, when rebuilt from the complete pool (restoring what the earlier analysis behind the slides had):
  1. 90% identity filter on the seeds (and one seed per genomic locus);
  2. leave-one-species-out evaluation;
  3. non-selenoprotein 3' regions and shuffled genomes as backgrounds;
  4. whole-genome and non-target (Ciona, Styela, Thalia) scans, to set the thresholds;
  5. covariation of the base pairs.
- Machinery (in progress): appendicularian Selenoprofiles profiles (<family>_app) built from the standard profiles
  plus appendicularian proteins (stage 1 predictions, curated transcripts), OdB3 held out; tested by Selenoprofiles
  on the nine genomes against stage 1; then deposited (Zenodo) and passed to the pipeline as extra profiles.

## Manuscript style

- Write "selenocysteine" when the amino acid is meant alone; keep "Sec" only inside names (tRNA-Sec, eEFSec, SecS)
  and in the first definition, "selenocysteine (Sec)".
- Give the gene symbol (human HGNC symbol) at the first mention of each gene, and the full name at the first mention
  of a symbol (for example "SECIS-binding protein 2 (SECISBP2)", "selenoprotein T (SELENOT)"); expand every acronym at
  its first use (for example "selenocysteine insertion sequence (SECIS)").
- Orthology of the machinery genes (to establish before the paper; outside the pipeline, in stage2/orthology/):
  1. first pass for every predicted gene: reciprocal best hit against a reference proteome that contains the related
     families (human RefSeq proteome, from the control genomes), recording the best and second-best hits;
  2. for the ambiguous cases of step 1 (decided 2026-10-17: reciprocal best hits first, phylogeny only where they
     are ambiguous), a maximum-likelihood tree (IQ-TREE) of the appendicularian proteins with orthologs from
     other tunicates and vertebrates and with the closest related family as outgroup (PSTK with KTI12; eEFSec with
     eEF1A and other translation GTPases; SEPHS2 with the SPS1-type paralogs of Mariotti et al. 2015; SECISBP2 with
     SECISBP2L; SECp43 with other RRM proteins such as TIA1); fast-evolving appendicularian sequences need models that
     resist long-branch attraction (for example LG+C20 or site removal);
  3. diagnostic features where they exist: selenocysteine or cysteine at the catalytic site of SEPHS2 (SPS2-type),
     the domain IV extension of eEFSec, the domain architecture (Pfam);
  synteny is of little use in O. dioica (scrambled genomes), except between its cryptic species.
- Discussion: a selenoprotein gene can be missed for three reasons that a survey must distinguish: the prediction
  method fails, the annotation pipeline fails, or the assembly is wrong at that locus (moved out of the introduction).
- Introduction TODO: confirm by a full-text screen of the large surveys that SEPHS2 (Mariotti et al. 2015) is the only
  appendicularian selenoprotein gene reported.
- Introduction, to revisit once the results are complete: the first paragraph (machinery; consider a table of the
  machinery genes with their roles and symbols); the importance of the O. dioica selenoproteome (placeholder in the
  text); the last paragraph (concrete contributions).
- Gene symbols: TRU-TCA1-1 (human tRNA-Sec gene) appears once, in the Methods, for text searches.
- Background noted: community metabolomics of marine plankton (Llewellyn et al. 2015 Prog Oceanogr, keywords
  "meta-metabolomics"; Heal et al. 2021 mSystems) and meta-omics models of plankton ecosystems (Stec et al. 2017);
  zooplankton-derived dissolved organic matter linked to zooplankton meta-transcriptomes (bioRxiv
  10.1101/2024.10.27.620474). Our genome annotations of zooplankton feed such approaches.
- Definition (agreed 2026-10-17): "the machinery" is the six components dedicated to selenocysteine and required for
  it: tRNA-Sec, SEPHS2, PSTK, SEPSECS (synthesis), EEFSEC, SECISBP2 (insertion). SECp43 (TRNAU1AP) is reported as an
  associated factor: not required in mouse liver and neurons (Mahdi et al. 2015), homologs in lineages without
  selenocysteine (yeasts, land plants in the Selenoprofiles profile), but reported by Selenoprofiles by default and
  absent from all oikopleurids (discussion). Other borderline genes (SEPHS1, SECISBP2L, SCLY, SARS1, ALKBH8, TRIT1,
  RPL30, nucleolin, EIF4A3) are excluded by the specificity criterion; SECISBP2 versus SECISBP2L is an orthology
  question for the reciprocal best hits.
- Reciprocal best hits (decided 2026-10-17): reference proteome = human RefSeq proteins (GCF_000001405.40, fetched
  with the genome; current, with gene symbols and the related families such as KTI12 and SECISBP2L); the RefSeq
  proteomes of OKI2018 and of the other annotated genomes serve for the reverse searches. Optional cross-check: the
  OIST Swiss-Prot BLAST/DIAMOND databases (2020_06, BioinfoUgrp modules), which are older and site-specific.
- Reciprocal best hits, first results (2026-10-18, human RefSeq proteome):
  - SecS (SEPSECS), eEFSec (EEFSEC) and SECp43 in F. borealis: reciprocal best hits; extra eEFSec copies only forward.
  - SECISBP2: in every tunicate, appendicularian or not, the single SBP2-family gene is closest to human SECISBP2L
    (margins -2 to -38 bits; H. roretzi +1), whereas vertebrates have both paralogs, cleanly separated. The tunicate
    gene is therefore a co-ortholog of the vertebrate SECISBP2 and SECISBP2L (duplication in vertebrates, to check
    in the literature, for example Donovan and Copeland 2009); the manuscript must name it accordingly.
  - SEPHS2: the selenocysteine-containing SEPHS of tunicates (and of zebrafish, weakly) are closer to human SEPHS1
    or only marginally to SEPHS2: best hits cannot separate SPS1 and SPS2 outside mammals (Mariotti et al. 2015).
    Use the residue at the catalytic site (selenocysteine or cysteine: SPS2 type) as the criterion.
  - PSTK: other tunicates are reciprocal best hits (55-108 bits); the O. dioica, O. albicans and B. stygius PSTK have
    no blastp hit in the human proteome at E <= 1e-3, O. vanhoeffeni 48 bits: too divergent for best hits against
    human; the absence test (tunicate references, KTI12 outgroup) and a phylogeny are the evidence.
  - Phylogenies needed: PSTK (with KTI12), SEPHS (SPS1/SPS2, with the Mariotti et al. 2015 sequences), SBP2 family
    (placement before the vertebrate duplication).
- Phylogenies, provisional (2026-10-18; sandbox, one CPU, LG+G4/R4 only; to redo on Deigo with the full model
  search and a site-heterogeneous model):
  - PSTK (154 sequences, rooted on KTI12): appendicularian PSTK form a clade (100) inside a tunicate group (74), not
    with KTI12; the PSTK/KTI12 split has only 54; the appendicularian clade is sister to insect-gut fungi (83), a
    long-branch artefact.
  - SEPHS (131, midpoint): vertebrate SEPHS1/SEPHS2 separate (100 each); tunicate selenocysteine copies (89) and
    glycine forms (97) separate; appendicularians have a second, non-selenocysteine SEPHS (Ala/Ser; Olon, Oalb, Ovan,
    OKI2018; 100) sister to the ascidian glycine form (95); the appendicularian SEPHS2 with selenocysteine (all nine)
    sit on a long branch outside the animal group. To check: the non-selenocysteine copy in OSKA2016 and OdB3; the
    Mariotti et al. 2015 SPS1/SPS2 sequence set.
  - SECISBP2 family (118, rooted on non-animal sequences; ClipKIT gappy 0.5, 304 columns; length filter 0.2 so that
    the domain-only tunicate predictions enter): oikopleurid proteins (99) nest inside the tunicate clade (87), sister
    to Thalia; tunicates group with protostome and echinoderm SECIS-binding proteins (39), and the vertebrate
    SECISBP2 and SECISBP2L clades (100 each) are not resolved against them (34-41): the tunicate gene is a single
    SECIS-binding protein that cannot be assigned to either vertebrate paralog. F. borealis is misplaced (13).
- Phylogenies, final (2026-10-19, Deigo, 8 CPUs; full ModelFinder: PSTK JTT+F+R6, SEPHS LG+R6, SECISBP2 family
  LG+I+R5; plus LG+C20+F+G for PSTK and SEPHS in results/trees/C20):
  - PSTK: appendicularian PSTK monophyletic (100) and inside PSTK, not KTI12 (PSTK/KTI12 split 94 with JTT+F+R6, 52
    with C20); their position inside PSTK is unstable, a long-branch effect (next to birds with JTT+F+R6, next to
    insect-gut fungi with C20, 81). Conclusion: PSTK, position within the family unresolved.
    CORRECTION (2026-10-19): the PSTK/KTI12 split has 100 with both models. The values 94 and 52 above, and the
    earlier provisional 54, came from a display error: Bio.Phylo moves supports by one branch along the rerooting
    path. Supports are now read from the unrooted bipartitions, and the figure scripts put them back on their
    bipartitions after rerooting (reroot_keeping_supports). Other values checked from the unrooted trees: PSTK
    appendicularian clade 100; SEPHS appendicularian selenocysteine clade 100 (C20 60), non-selenocysteine clade 100,
    with the tunicate glycine forms 94 (C20 95), tunicate selenocysteine copies 91 (C20 56); SECIS-binding proteins:
    oikopleurid clade 100, with Nematostella 70, other tunicates 90.
  - SEPHS: the non-selenocysteine appendicularian copies (Olon, Oalb, Ovan, OKI2018) group with the tunicate glycine
    forms (94; C20 95). The appendicularian selenocysteine SEPHS2 form their own clade (100, with F. borealis) at the
    base of the eukaryotic SEPHS (90 is the support of the whole eukaryotic group, not of a grouping with the glycine
    forms: corrected after drawing the collapsed tree); their position is unresolved (long branch). The
    selenocysteine copies of Styela, Halocynthia and Thalia group together (91; C20 56).
  - SECISBP2 family: the oikopleurid proteins form a clade (100, with the O. dioica profile sequence), sister to
    Nematostella (70), not with the other tunicates (90, among protostome and echinoderm proteins, 18): unlike the
    provisional tree, no support for grouping appendicularians with other tunicates. Conclusion: a member of the
    SECIS-binding protein family, orthology to the tunicate gene and to the vertebrate paralogs unresolved.
- Display trees (figures/phylogeny_<family>_collapsed.*, figures/phylogeny_collapsed.py): nodes with ultrafast
  bootstrap < 50 collapsed into polytomies; clades of a single taxonomic group (config/taxon_groups.tsv, assigned
  by hand from genus names) drawn as one tip; appendicularian and other tunicate sequences always shown.
- Long-branch analyses (stage2/machinery_lba.py, make machinery-lba; to run on Deigo): LG+C60+F+G with PMSF,
  slow-fast removal of the fastest sites (10 to 40 %), removal of non-target long branches (insect-gut fungi in
  PSTK), and an AU test of the grouping of appendicularians with the other tunicates. Background: O. dioica proteins
  evolve about twice as fast as those of Ciona, with homogeneous divergence across amino acids and a fast-evolving
  cysteine (Berna et al. 2012, doi:10.1016/j.ympev.2011.11.013; Berna and Alvarez-Valin 2015,
  doi:10.1016/j.margen.2015.07.007); the cysteine result also matters for the selenocysteine/cysteine homologs.
- Taxon sampling of the machinery trees (2026-10-19): added four tunicates (Botryllus schlosseri GCF_051294905.1,
  Phallusia mammillata GCF_965637545.1, the salps Pegea sp. GCA_034620365.1 and Salpa thompsoni GCA_047495815.1)
  to break the appendicularian long branch, and two cephalochordates (Branchiostoma floridae GCF_000003815.2,
  B. lanceolatum GCF_035083965.1) as chordate anchors. No echinoderm or hemichordate: the deuterostome clade itself
  may be an artefact of unequal branch lengths and compositional heterogeneity (Serra Silva et al. 2025, Current
  Biology, doi:10.1016/j.cub.2025.06.045), whereas chordates are well supported. SEPHS tree: the SPS1/SPS2-labelled
  sequences of Mariotti et al. 2015 (Supplemental Material S8; curation/mariotti2015_sps.faa). PSTK and
  SECIS-binding proteins: RefSeq proteins replace the partial Selenoprofiles predictions where a genome has them.
  Profile and reference sequences capped at three per taxonomic group (and SPS type).
- Machinery trees with the extended taxon sampling (2026-10-19, Deigo, commit 149cf6d; collapsed figures redrawn):
  - PSTK (86 sequences, 691 sites, LG+F+R6): PSTK/KTI12 split 100. Appendicularian clade 100 (OKI2018, OSKA2016,
    OdB3, O. albicans, O. vanhoeffeni; the O. longicauda, B. stygius and M. erythrocephalus fragments are below the
    length filter), sister to all other PSTK (91; C20 91). Other tunicates (ascidians and salps) form a clade (99)
    elsewhere in PSTK. Long-branch analyses: appendicularian clade 100 and PSTK/KTI12 100 in every analysis; never
    with the other tunicates; sister to the other PSTK with PMSF LG+C60 (90), 10-30 % of fast sites removed (86-92)
    and without the insect-gut fungi (89); with 40 % removed, next to the fungi (35). AU test: the tree with
    appendicularians and other tunicates together is not rejected (p-AU 0.275, delta logL 5.7). Conclusion: PSTK;
    position within the family unresolved; grouping with the other tunicates not rejected.
  - SEPHS (129 sequences, 976 sites, LG+R6): appendicularian selenocysteine SEPHS2 clade of 11 (all nine genomes
    except B. stygius, with F. borealis, M. erythrocephalus and the O. dioica SPS2 of Mariotti et al. 2015) 83 (C20
    89), position unresolved (32; C20 56). Non-selenocysteine copies (OKI2018, O. longicauda x2, O. albicans,
    O. vanhoeffeni) 100, inside the tunicate glycine forms (95; C20 100), next to Phallusia (78; C20 73). Ascidian
    selenocysteine copies 100; with the salps and Thalia 93 with C20 only. All tunicate SEPHS except the
    appendicularian selenocysteine clade: 78 (C20 92). Reading: the tunicate selenocysteine/glycine duplication
    precedes the divergence of appendicularians, so the appendicularian selenocysteine SEPHS2 is expected next to the
    ascidian and salp selenocysteine copies; it sits outside them on a long branch. Long-branch analyses pending.
  - SECIS-binding protein family (90 sequences, 391 sites, LG+R5): Oikopleuridae clade (with B. stygius and
    M. erythrocephalus and the O. dioica profile sequence CBY19450.1) 100, next to Nematostella (48), outside a group
    of the other tunicates, vertebrates and other bilaterians (71); F. borealis unresolved (34). Other tunicates 100,
    branching before the vertebrate SECISBP2/SECISBP2L duplication (vertebrate clade 86). Long-branch analyses
    pending.
- stage2/machinery_lba.py corrected (2026-10-19): (1) the p-AU value was read from a line of the rate table (0.1646
  for PSTK) instead of the AU table (0.275); (2) leaves were assigned from the groups table, which left out the
  reference and profile sequences of appendicularians and tunicates (O. dioica SPS2 and CBY19450.1, Botryllus SPS1
  and SPS2, Ciona XP_002129526) and put both appendicularian SEPHS paralogs in one set; the AU constraint for SEPHS
  and the SECIS-binding proteins therefore forced these sequences out of their clades. Now: leaves assigned by taxon
  (config/taxon_groups.tsv), the appendicularian clade is the largest appendicularian clade of the ML tree, the
  other tunicates are those of the orthologous form (SEPHS glycine forms left out), and a constrained tree that does
  not satisfy the current constraint is recomputed with its AU test. PSTK summary regenerated (trees unchanged).
- Long-branch analyses of SEPHS and the SECIS-binding proteins (2026-10-20, Deigo, commit 72da760;
  results/trees/lba/<family>.summary.tsv):
  - SEPHS: the appendicularian selenocysteine SEPHS2 clade holds with PMSF LG+C60 (79) and with 30 % of fast sites
    removed (84); with 10, 20 and 40 % removed, F. borealis leaves it (the other ten stay together, 100). Its position
    is unstable: at the base of the animal SEPHS in the ML, C20 and PMSF trees (32-56), and outside a group that
    includes kinetoplastid, stramenopile and amoebozoan sequences in the slow-fast trees, a position not possible
    for an animal gene. The tunicate SEPHS clade without it is supported in every analysis (ML 78, C20 92, PMSF 97),
    with the appendicularian non-selenocysteine copy inside the glycine forms (95-100; next to Phallusia 73-98)
    and the ascidian and salp selenocysteine copies together with C20 (93) and 40 % removal (83). AU test: the tree
    in which the appendicularian selenocysteine SEPHS2 groups with the ascidian and salp selenocysteine copies is
    rejected (p-AU 0.014, delta logL 52; LG+R6); in that constrained tree the glycine forms move next to the
    vertebrates. Reading: the tunicate duplication into selenocysteine and glycine forms precedes the divergence of
    appendicularians (their non-selenocysteine copy is a glycine-form ortholog); the appendicularian
    selenocysteine SEPHS2 is too divergent for its position to be inferred, and the AU rejection is made under a
    site-homogeneous model, exposed to the same long-branch effect. Its assignment rests on the selenocysteine at
    the catalytic site (SPS2 type; Mariotti et al. 2015) and on the presence of the separate glycine-form ortholog.
  - SECIS-binding protein family: the Oikopleuridae clade holds in every analysis (99-100). With PMSF LG+C60 it
    groups with the other tunicates (81), and with 40 % of fast sites removed too (71, F. borealis joining it,
    71); in the ML tree and with 10-30 % removed it sits next to Nematostella (below 50). AU test: grouping with the
    other tunicates not rejected (p-AU 0.503, delta logL 0.1). The vertebrate SECISBP2 and SECISBP2L form one clade
    that excludes the tunicates in every tree (86-100). Reading: the appendicularian protein is the ortholog of the
    single tunicate SECIS-binding protein gene, which branches before the vertebrate duplication (a pro-ortholog of
    both SECISBP2 and SECISBP2L); support comes from the site-heterogeneous model, not from the ML tree.
  - Summary of the three families: PSTK and SECIS-binding protein: appendicularian placement with the other
    tunicates not rejected (p-AU 0.275 and 0.503), supported for the SECIS-binding protein by PMSF C60 (81).
    SEPHS2: placement unresolved; the non-selenocysteine copy settles the timing of the tunicate duplication.
- SEPHS orthology AU tests (2026-10-20, Deigo, commit 559d9c2; results/trees/lba/SEPHS2.{ml,pmsf}.*): both positions
  of the appendicularian selenocysteine SEPHS2 clade are rejected under both models. Inside the clade of all tunicate
  SEPHS: p-AU 0.0067 (LG+R6, delta logL 43) and 0.019 (LG+C60+F+G PMSF, delta logL 34); next to the ascidian and salp
  selenocysteine copies, with the glycine forms and the second appendicularian copy as their sister: p-AU 0.012
  (delta logL 51) and 0.0082 (delta logL 39). Checks: no appendicularian selenocysteine SEPHS2 fails the IQ-TREE
  composition test (median p 74 %); its path from the common node of the tunicate and vertebrate SEPHS (2.2
  substitutions per site) is shorter than that of the second, non-selenocysteine appendicularian copy (3.0), which is
  placed inside the tunicate glycine forms in every analysis. Reading: under both models the data reject the
  appendicularian selenocysteine SEPHS2 as a member of the tunicate SEPHS clade, and neither composition nor the
  longest branch explains it. Two readings remain: a divergent ortholog of the tunicate selenocysteine copy whose
  signal these models do not recover, or a separate SEPHS2 lineage (with the ortholog of the tunicate selenocysteine
  copy lost in appendicularians). The manuscript states the rejection and leaves the origin open; the protein is a
  selenocysteine SEPHS2 by its catalytic residue in any case.
- Curation of the machinery tree sets (2026-10-21, after review of the PSTK alignment by eye: a truncated lancelet
  PSTK, rare insertions in mammals): stage2/machinery_phylogeny.py now (1) leaves out the human stage 1 predictions
  (the human RefSeq proteins represent human; the PSTK prediction carried an N-terminal extension and the SEPHS
  predictions included pseudogenes); (2) keeps identical sequences once (genome or RefSeq names preferred); (3) after a
  first alignment, removes sequences covering less than half of the core columns, except appendicularian, curated and
  tunicate sequences, then realigns; (4) trims with ClipKIT gappy 0.9 instead of smart-gap (smart-gap had kept 224 of
  691 PSTK columns and 426 of 976 SEPHS columns occupied by at most 6 sequences, for example the 75-residue insertion
  of mammalian KTI12, which is in the canonical human and mouse proteins). Removals are listed in
  results/trees/<family>.removed.tsv. Preview on the previous sets (local MAFFT 7.505): PSTK 83 sequences, 413
  columns (removed: human prediction, duplicate human KTI12, Branchiostoma floridae PSTK fragment of 106 residues);
  SEPHS 117, 448 (6 human predictions, 6 duplicates); SECIS-binding proteins 82, 427 (2 human predictions; fragments
  of Nematostella (78 residues, the ML sister of the oikopleurid clade), Helobdella, Trichoplax, Daphnia, Bodo and
  Sphaeroforma). Trees, long-branch analyses, figures, tables and manuscript numbers to recompute.
- Machinery trees after curation (2026-10-22, Deigo, commit 806e317; alignments PSTK 83 x 419, LG+F+I+R5; SEPHS
  117 x 454, LG+R6; SECIS-binding proteins 81 x 444, LG+I+R6). Removed: duplicate human KTI12 and the B. floridae
  PSTK fragment; 6 SEPHS duplicates; 7 partial SECIS-binding proteins (KOO29690.1, Daphnia, Bodo, Sphaeroforma,
  Trichoplax, Helobdella, Nematostella); human predictions (not listed in this run's removed.tsv; recorded from now on).
  - PSTK: PSTK/KTI12 100 and appendicularian clade 100 everywhere; other tunicates 97-100; the appendicularian clade
    now groups with the PSTK of the insect-gut fungi (56 ML, 74 C20, 84 PMSF, 60-77 site removal), and without the
    fungi it is sister to all other PSTK (94): long-branch attraction between two long branches. AU p 0.268.
  - SEPHS: oikopleurid selenocysteine SEPHS2 (10 sequences, 7 genomes) 100 everywhere; F. borealis joins under C20
    and PMSF only (83, 76). Placement: base of animal SEPHS under C20 and PMSF (39, 56), outside a group with
    kinetoplastid, stramenopile and amoebozoan proteins in ML and every site-removal tree. Tunicate SEPHS without it
    82-95; glycine forms with the second appendicularian copy 75-99; Phallusia sister 79 (ML), 95 (PMSF). AU: all
    rejected, more strongly than before (inside tunicate SEPHS 0.011 ML, 0.0045 PMSF; orthologous position 0.0047,
    0.0009; generic 0.0057). No SEPHS fails the composition test; paths from the tunicate-vertebrate common node:
    appendicularian selenocysteine copy 1.3, second copy 2.5 substitutions per site.
  - SECIS-binding proteins: Oikopleuridae clade 99-100; F. borealis joins it in all but the 40 % site-removal tree
    (47-77). Appendicularians with the other tunicates: 70 ML, 93 PMSF, 53-68 site removal (oikopleurids with
    tunicates 76 at 40 %); AU p 0.586. Vertebrate SECISBP2+SECISBP2L clade excluding tunicates and amphioxus 79-98.
    The earlier ML placement next to Nematostella was due to a 78-residue fragment.
  - Provisional Results, Table MACHINERY and Methods updated with these numbers; collapsed figures redrawn (title on
    two lines, legend below the axis label). The C20 run now uses --max-per-group 3 for its sequence collection too
    (its trees always used the copied trimmed alignments; only its unused faa and groups files had all sequences).
- Repository brought in line with the curated machinery analysis (2026-10-22): Table MACHINERY is now generated by
  stage2/machinery_table.py (make machinery-table, also run at the end of make machinery-lba) from the census, the
  reciprocal best hits, the long-branch table and the SEPHS trees, with the orthology column and the second SEPHS row
  computed from them; the manuscript table is pasted from results/tables/stage2_machinery_table.md. Manuscript: the
  machinery paragraph calls the appendicularian protein a SECIS-binding protein co-orthologous to SECISBP2 and
  SECISBP2L. READMEs: top-level steps (actual make targets), data sources (all 24 genomes, Mariotti et al. 2015
  references) and layout; stage2/README.md (curation, trees, long-branch analyses, table, run times); figures and
  results/tables indexes. Still open: results/tables/stage1_machinery_rbh.tsv (make machinery-rbh on Deigo); the
  second SEPHS copy in OSKA2016 and OdB3.
- Domain analysis of the machinery trees (2026-10-22, local run, 2 CPUs; stage2/machinery_domains.py, make
  machinery-domains): each trimmed alignment split at a domain boundary of the human protein (InterPro entries of
  Q8IV42, Q99611, Q96T21; curation/machinery_domains.tsv), the parts analysed on the same sequences. PSTK kinase
  1-263 (294 columns) / C-terminal (125), 82 sequences; SEPHS N-terminal 1-238 (234) / C-terminal (220), 116; SECIS-
  binding protein upstream of L7Ae 1-631 (253) / L7Ae 632-854 (191), 76.
  - In every family, each part rejects the tree of the other (AU p <= 0.0023), but no clade with >= 95 % support in
    one domain tree is incompatible with one so supported in the other, except arrangements within taxonomic groups
    and among prokaryotic and alveolate SEPHS.
  - PSTK: appendicularian clade 100 (kinase) and 97 (C-terminal), without the fungi; sister to all other PSTK (57) or
    next to dipteran PSTK (54); appendicularians with tunicates not rejected (0.44, 0.32). C-terminal data reject the
    full-length tree (0.0001), kinase data do not (0.19).
  - SECIS-binding protein: L7Ae tree as the full-length tree (oikopleurids with tunicates 81; vertebrate SECISBP2 +
    SECISBP2L 95); upstream region: oikopleurids next to vertebrate SECISBP2, with tunicates (67, 67). With tunicates
    not rejected (0.38, 0.49). F. borealis (315 aa) elsewhere in both.
  - SEPHS: N-terminal data reject no tested position (with tunicate Sec copies 0.45, inside all tunicate SEPHS 0.33,
    orthologous position 0.19); C-terminal data reject the orthologous position (0.0034; others 0.08). The N-terminal
    tree groups both appendicularian copies with dipteran and louse SPS2 (84; second copy with them 94), against the
    C-terminal tree (second copy in tunicate glycine forms 77-83; insect SPS2 with Rhodnius SPS2 100): the only
    conflict >= 90 % involving an appendicularian protein. The rejection of SEPHS2 orthology rests on the C-terminal
    half.
  - Provisional Results and Methods updated; display trees figures/phylogeny_<family>_collapsed_domain_<domain>.pdf;
    figures/phylogeny_collapsed.py got --formats.
- Schematic of the selenocysteine machinery (2026-10-22): figures/machinery_scheme.py (make machinery-scheme), drawn
  from scratch with SBGN-style glyphs (metabolite, RNA, protein, complex; conversion, catalysis, binding); PDF and
  editable SVG; draft legend in figures/README.md. Biochemistry checked against the abstracts of Xu2006 and
  SmallHoward2006. Note: SmallHoward2006 reports SECp43 complexes with SLA/LP (SEPSECS) and SPS1, not SPS2; the
  manuscript sentence "in complexes with SecS and SEPHS2" needs checking.
