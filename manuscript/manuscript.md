---
title: "The selenoproteome of *Oikopleura dioica* and other appendicularians"
bibliography: references.bib
---

<!-- Draft. Citations use Pandoc keys from references.bib; [REF: …] marks references still to add. -->

# Introduction

<!-- Part 1: justification -->

Selenoproteins carry selenocysteine (Sec), an amino acid encoded by UGA, which is otherwise a stop codon. Because
gene prediction programs read UGA as a stop, selenoprotein genes are often truncated, split or missing in genome
annotations [@Mariotti2015]. In vertebrates, a systematic comparison found that only 5% of selenoprotein genes are
correctly annotated in GenBank and 11% in Ensembl; RefSeq, which applies a dedicated selenoprotein pipeline, reaches
77%, with most remaining errors in families whose Sec is near the C-terminus [@Tico2026]. Outside vertebrates, and for
genomes without such a pipeline, the situation is largely unmeasured.

The appendicularian *Oikopleura dioica* is a planktonic tunicate, part of the sister group of vertebrates. Its short
life cycle, transparent body and invariant embryonic cell lineages make it a model for developmental biology
[@Wang2020], and its very small and compact genome makes it a model for genome evolution [@Denoeud2010;
@Naville2019]. Three chromosome-scale or scaffold-level assemblies are public: OdB3 from a Norwegian population
[@Denoeud2010], OSKA2016 from a Japanese population [@Wang2020] and OKI2018 from an individual collected in Okinawa
[@Bliznina2021]. These populations differ so much in genome structure that they have been proposed to be cryptic
species [@Plessy2024]. Each assembly has its own gene annotation, produced by a different pipeline. The RefSeq
annotation of OKI2018 comes from the NCBI pipeline, which handles selenoprotein genes well in vertebrates
[@Tico2026]; how the three annotations represent selenoprotein genes in *O. dioica* has not been examined.

A census of selenoprotein genes is part of the basic description of a model organism. Selenoproteins include
glutathione peroxidases, thioredoxin reductases and iodothyronine deiodinases, enzymes of redox homeostasis and of
thyroid hormone metabolism; studies of these functions in *O. dioica* start from the list of genes and their correct
protein sequences. The census also matters beyond *O. dioica*: selenoprotein sets differ between animal lineages
through gene losses and replacements of Sec by cysteine [@Mariotti2012; @Jiang2012], and appendicularians, with
their fast-evolving genomes, have not been surveyed. Previous studies included *O. dioica* only for single families:
its selenophosphate synthetase SEPHS2 is a Sec gene with a SECIS element [@Mariotti2015].

<!-- Part 2: background for the reader -->

Sec insertion requires a dedicated machinery. In eukaryotes, a stem–loop structure in the 3′ untranslated region of
selenoprotein mRNAs, the SECIS element, is bound by SECIS-binding protein 2, which recruits the Sec-specific
elongation factor eEFSec and the Sec-charged tRNA-Sec to the UGA codon [REF: review of eukaryotic Sec insertion].
Sec is synthesised on its tRNA: seryl-tRNA-Sec is phosphorylated by PSTK, then converted to selenocysteinyl-tRNA-Sec
by SecS, using selenophosphate produced by SEPHS2 [REF: Sec biosynthesis on tRNA]. Additional RNA structures next to
the UGA, the Sec redefinition elements, can stimulate UGA recoding [@Howard2005]. A genome that encodes
selenoproteins is therefore expected to also encode tRNA-Sec, SEPHS2, PSTK, SecS, eEFSec and SECIS-binding protein 2.

Selenoprotein genes are predicted by homology rather than by ab initio gene finding. Selenoprofiles aligns curated
protein profiles of each selenoprotein family to a genome and labels the predictions with Sec, cysteine or other
residues at the Sec position [@Mariotti2010; @Tico2026]. Selenoprofiles does not require a SECIS element; SECIS
elements are searched separately, for example with covariance models such as the Rfam SECIS family, whose seed
derives from a survey of metazoan SECIS elements [@Lambert2002; REF: Rfam; REF: Infernal], or with SECISearch3
[REF: SECISearch3]. tRNA-Sec genes are found with dedicated covariance models [@Santesmasses2017]. The Rfam SECIS
seed contains no tunicate sequence, and how these methods perform on divergent genomes is not known in advance.

Tunicate selenoproteins are known mostly from ascidians. The selenoproteome of *Ciona intestinalis* comprises 19
genes [@Jiang2010], and SelenoDB 2.0 provides automatic annotations of selenoprotein genes and of their cysteine
homologs for 58 animal genomes, including *Ciona* [@Romagne2014]. A Sec-containing deiodinase with a 3′ UTR SECIS
element was characterised in *Halocynthia roretzi* [@Shepherdley2004], and in ascidians the selenophosphate
synthetase gene evolved through alternative isoforms and a retroposed copy [@Mariotti2015]. A survey of the thyroid hormone toolkit
across chordates reported one deiodinase per non-vertebrate chordate and noted that tunicate deiodinase sequences in
public databases often lack parts of the protein, calling for better annotation [@Esposito2021].

*O. dioica* adds difficulties to selenoprotein annotation. Its genome is compact, and many of its genes are
transcribed as operons whose polycistronic transcripts are resolved by trans-splicing of a spliced leader to the
downstream cistrons; the same leader is also added to the 5′ end of many monocistronic transcripts [@Ganot2004;
@Wang2015]. Part of its introns use non-canonical GA–AG splice sites [@Denoeud2010], which gene prediction programs
built for canonical GT–AG introns model poorly, so that de novo gene models are often wrong at these introns. Gene
models from the three assemblies differ, and assemblies can contain gaps, duplicated haplotypes and local errors. A selenoprotein gene can therefore be missed for three reasons that a survey
must distinguish: the prediction method fails, the annotation pipeline fails, or the assembly is wrong at that
locus.

<!-- Part 3: rationale for stage 1 -->

Here we establish the selenoproteome of *O. dioica*, and we use the genomes of six other appendicularian species,
together with the selenoproteomes of other tunicates, to extrapolate a baseline for appendicularians. We also
document what is needed to annotate these genes reliably. We proceed in two stages. In stage 1, we apply current
tools with their default parameters. For this, we wrote a pipeline that combines Selenoprofiles, protein alignments
with miniprot [REF: miniprot], the Rfam SECIS model and tRNA-Sec searches, and that compares the predicted Sec codons
with the public annotations and with transcript evidence. Stage 1 answers two questions: how many selenoprotein genes do
default tools find, and which ones do they miss? To answer them, we built a curated baseline from the three
*O. dioica* assemblies. Because they come from closely related cryptic species, a gene present in two assemblies and
absent from the third is first checked for an assembly or annotation problem, and transcripts from each population
test the Sec codons directly. The genes and SECIS elements that stage 1 misses define the changes tested
in stage 2: a SECIS covariance model built from appendicularian elements, lineage-specific reference proteins and
profile searches for divergent genes. Stage 2 then gives the selenoproteome of *O. dioica* and, from the genomes of
six other species [@Naville2019], whose assemblies are often too fragmented for a complete census, a baseline for
appendicularians. The *Ciona* selenoproteome [@Jiang2010; @Romagne2014] serves as the tunicate reference to which
this baseline is compared.

<!--
Notes for the next pieces
- Data: the NCBI assembly OKI2018_I68 (GCF_907165135.1) is the OKI2018_I69 assembly of Bliznina et al. 2021 (I69 is
  the sequenced individual; I68 is a typo in the NCBI record); we use the NCBI name.
- References to add: review of Sec insertion; Sec biosynthesis (PSTK, SecS); Rfam; Infernal; SECISearch3; miniprot;
  minimap2; Nextflow; nf-core. Selenoprofiles 4 is described in Tico2026.
- SelenoDB 2.0 and the Selenoprofiles expectation table were used for the chordate comparison (Ciona gene lists;
  vertebrates); say so in the Methods.
-->
