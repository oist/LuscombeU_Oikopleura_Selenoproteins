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
