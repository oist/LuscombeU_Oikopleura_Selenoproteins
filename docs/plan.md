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
