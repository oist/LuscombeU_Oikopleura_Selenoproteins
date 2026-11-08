# EEFSEC: text removed on 2026-11-08

The EEFSEC sentences of the gold-standard section ("Five models have splice sites outside the intron classes ...")
were replaced by a shorter text in the new subsection "SEPSECS and EEFSEC" of the machinery section, derived from a
replacement paragraph by Charles (see the commit that created this file); the gold-standard section keeps a pointer.
Figure EEFSEC moved with them and its caption was rewritten. The removed text below is verbatim.

## What the removed text holds, and how to regenerate it

- The retrotransposon of OdB3 EEFSEC (make eefsec-insertion, curation/insertion.py, results/tables/eefsec_insertion.tsv):
  intron of 8,982 nt, insertion of 8,977 nt, its duplicated target site giving the intron its equivalent placements;
  longest open reading frame of 1,244 codons with a reverse transcriptase motif (YLDD) followed by an integrase-like
  HHCC motif (orf1 of the table; three more ORFs of 241 to 651 codons).
- The OdB3 EST that joins the two exons: from an allele without the insertion or from a transcript that splices it
  out; the protein is the same.
- Copies of the element (results/tables/eefsec_insertion_copies.tsv): all loci, coverage and identity, in the six O.
  dioica assemblies.
- OKI2018 EEFSEC: the 73-nt intron of the prediction around the frameshift, the 50-nt first intron of the transcripts
  at the position of the first intron of OSKA2016 and OdB3, the model taking the protein of a transcript (523
  residues, the length of the OSKA2016 and OdB3 proteins).
- Divergence of the coding sequences between the two individuals of each species for the five machinery genes (make
  cds-divergence, results/tables/cds_divergence.tsv): the values of the two other species (EEFSEC 6% and 5% of third
  positions; other genes 0.3 to 4.5%), and the caveat that with two individuals per species a single gene can diverge
  this much by chance, so that this does not by itself point to selection or to a rearrangement.

## Superseded statements (do not reuse as they stand)

None.

## Removed text (gold-standard section)

In OdB3 EEFSEC (Figure EEFSEC), the unusual intron (8,982 nucleotides)
is an insertion inside an exon that is continuous in OKI2018, OSKA2016 and in another assembly of the same species
(Bar2_p4): a long terminal repeat (LTR) retrotransposon of 8,977 nucleotides, whose two LTRs of 184 nucleotides are
identical, flanked by a duplication of five bases of the exon (GATAA), which gives the intron its equivalent placements;
its longest open reading frame (1,244 codons) has a reverse transcriptase motif (YLDD) followed by an integrase-like
HHCC motif. The EST that joins the two exons can come from an allele without the insertion as well as from a
transcript that splices it out; the protein is the same. Identical LTRs point to a recent insertion; the OdB3 assembly
holds another copy over 92% of its length (99.7% identity), two partial copies at contig ends (99.6 and 98.8%) and
fragments at five loci (90 to 94%); Bar2_p4 holds a complete copy (97.0%) and five fragments, all on its Y-specific
region (sequence YSR); the four assemblies of the two other species (OKI2018, KUM-M3-7f, OSKA2016, AOM-5-5f) have
none (BLASTN, E at most 10^-20^; `make eefsec-insertion`, `curation/insertion.py`). In OKI2018 EEFSEC, the assembly
has a run of five T where the four transcripts and the KUM-M3-7f assembly have six; the prediction splices around the
frameshift with a 73-nucleotide intron and reads 5 residues where the OSKA2016 and OdB3 proteins have 13, whereas the
transcripts splice a 50-nucleotide intron at the position of the first intron of OSKA2016 and OdB3; the model takes the
protein of a transcript (523 residues, the length of the OSKA2016 and OdB3 proteins). These transcripts differ from the
assembly at 17 to 52 positions per kb, and from each other, against at most 14 for the OKI2018 transcripts of the other
machinery genes: of five machinery genes, EEFSEC has the most divergent coding sequences between the two individuals of
each species, at third codon positions 16% between OKI2018 and KUM-M3-7f (0.3 to 5% for the other four genes), 6% and
5% within the two other species (0.3 to 4.5%), with few amino-acid changes (1 to 2%; `make cds-divergence`). The
transcripts thus carry other haplotypes of a gene whose alleles diverge more than usual in the Okinawa species; with
two individuals per species, a single gene can diverge this much by chance, so this does not by itself point to
selection or to a rearrangement. <!-- results/tables/cds_divergence.tsv (curation/cds_divergence.py, 2026-11-06):
codons through the gold-standard protein of each species aligned with exonerate; between North Pacific and North
Atlantic, EEFSEC third positions 31%, second of the five genes (13 to 39%). Mismatch rates per kb from the transcript
rows of figures/inputs/reannotation_<gene>.tsv. Tests noted in docs/plan.md: divergence in windows along
NC_137291.1 (OKI2018 against KUM-M3-7f), population data. -->

## Removed text (Figure EEFSEC and its comment, from the gold-standard section)

<!-- Figures EEFSEC, SBP2 and SEPHS_b: provisional (2026-11-07), placed here for Charles's review of the machinery
genes; they may be removed or aggregated later. Panels drawn by make machinery-reannotation (figures/reannotation_page.py
on figures/inputs/reannotation_<gene>.tsv): EEFSEC from run 5829a72, SBP2 and SEPHS_b from run 0d19c70; the later
runs changed only models of other species, which the panels do not draw. -->

![Figure EEFSEC](../figures/machinery_EEFSEC.png)

Figure EEFSEC (provisional; `figures/machinery_EEFSEC.pdf`, `make machinery-reannotation`). EEFSEC in the three *O.
dioica* assemblies, drawn as Figure PSTK A and B. Top, public models, gold-standard models and supporting transcripts
with their splice sites and their differences to the assembly. OKI2018: the gold-standard model follows the
transcripts (first intron of 50 nucleotides, at the position of the first intron of OSKA2016 and OdB3), whose many
mismatches come from other haplotypes (see text). OdB3: the long intron is the LTR retrotransposon inserted in an exon;
the two ESTs that reach the start codon have their 5′ ends masked in the EST collection used, so their spliced leader
is not tested. Bottom,
reference and gold-standard proteins.
