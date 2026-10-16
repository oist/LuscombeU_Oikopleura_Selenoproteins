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
  incomplete haplotype purging or from contamination. TODO: ask Charles for the reference documenting this, and cite
  it where copy numbers are discussed.
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
