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

<!-- Part 3: aims and approach -->

Here we establish the selenoproteome of *O. dioica*, and we use the genomes of six other appendicularian species,
together with the selenoproteomes of other tunicates, to extrapolate a baseline for appendicularians. Because no
tool had been tested on these genomes, we proceeded in two stages. We first measured what current methods find with
their default parameters, against a baseline curated from the three *O. dioica* assemblies; the genes and SECIS
elements they missed defined the adaptations tested in the second stage.

# Results

<!-- Stage 1 gene recovery (baseline_vs_stage1) to be written before this subsection. Figures to be numbered. -->

## The Rfam SECIS model does not separate appendicularian SECIS elements from genomic background

In stage 1, no SECIS element reached the Rfam gathering threshold (35 bits) downstream of any selenoprotein gene
predicted in the nine genomes. To see whether the threshold alone was too strict, we scored with the same model,
without threshold, the 3′ sequences of all Selenoprofiles predictions and, as background, 200 random genomic windows
of the same length per genome (Figure SECIS-1). Predictions fall into two groups: 143 with a Sec codon, and 578
others. The others are mostly homologs with cysteine at the Sec position (375, of which 235 in the EhSEP2 family),
predictions whose alignment does not cover the Sec position of their profile (122, labelled "unaligned"), homologs of
Sec machinery proteins that are not selenoproteins (40, such as eEFSec and SECIS-binding protein 2), and predictions
with another amino acid or a stop codon at the Sec position (41).

Of the 143 3′ sequences that follow predictions with a Sec codon, 96 (67%) had no hit at all, even when every hit
of at least 0 bits was reported. This was not caused by fragmented assemblies separating coding sequences from their UTRs: in the
chromosome-scale assembly OKI2018, 5 of 12 had no hit. The remaining 47 scored from 12 to
31 bits; 18 reached 20 bits, and none reached the threshold. Random windows rarely had a hit (45 of 1,800), but those
hits reached 22.6 bits. Lowering the threshold to 20 bits would therefore have recovered at most four SECIS elements
per genome while accepting background hits, and most selenoprotein 3′ sequences would still have no hit. The Rfam
model does not recognise most appendicularian SECIS elements, rather than recognising them below its threshold.

## An appendicularian SECIS model built from the stage 1 predictions

We built a new covariance model from the stage 1 outputs alone, by iterative bootstrap (Methods). Training used the
110 3′ sequences that follow predictions with a Sec codon in seven oikopleurid genomes: OKI2018 and OdB3 for
*O. dioica*, *O. longicauda*, *O. albicans*, *O. vanhoeffeni*, *B. stygius* and *M. erythrocephalus*. The OSKA2016
assembly of *O. dioica* and the genome of *F. borealis*, a fritillariid, were kept out of training as controls. The
initial seeds were the 26 training sequences whose best Rfam hit reached 15 bits. After five rounds, in which each
new model was used to collect the training hits that score above every background window, the set of seeds no
longer changed: 73 hits, in 73 of the 110 training sequences. The seeds lie close to the predicted stop codons (median
152 nucleotides downstream; quartiles 74 and 388) and are 60 to 75 nucleotides long.

We evaluated the final model on new random windows and on the curated *O. dioica* Sec genes (Figure SECIS-2). Of
1,800 random windows, 8 had a hit of at least 20 bits and 1 of at least 25 bits (27.3 bits, in *F. borealis*); in
the three *O. dioica* assemblies, no random window exceeded 20.8 bits. Hits downstream of Sec genes reached 62 bits.
In the held-out OSKA2016 assembly, 12 of the 14 Sec predictions and 15 of the 18 curated Sec genes had a hit of at
least 25 bits, against a best background window of 18.1 bits. In the training genomes, 12 of 12 Sec predictions of
OKI2018 and 18 of 21 of OdB3 scored at least 20 bits, but these sequences were used for training and overstate the
performance. In the held-out *F. borealis*, only 8 of 19 Sec predictions reached 20 bits, and 3 reached 25 bits.

Across the three *O. dioica* assemblies, 42 of 53 curated Sec genes scored at least 25 bits. Below 25 bits were
SELENOT in all three assemblies, SEPHS2 in OKI2018 and OdB3, SELENOO in OKI2018 and OSKA2016, SELENOK in OKI2018 and OdB3 (23.4 and
24.9 bits), and SELENOF in OSKA2016 and OdB3. For SEPHS2 and SELENOT, the Sec codon lies far upstream of the stop
codon, so that the 6,000 nucleotides scored after the Sec codon may not reach the SECIS element. Outside the Sec
predictions, 19 of the 578 other Selenoprofiles predictions had a hit of at least 25 bits. Seventeen of them belong to
selenoprotein families and either do not cover the Sec position ("unaligned": GPX, DIO, SEPHS2, SELENOT and others) or
carry another residue there (SELENOW with isoleucine, tyrosine or leucine), as expected if some are Sec genes that
Selenoprofiles did not recognise as such, including SEPHS2 of OSKA2016 (labelled "pseudo"). Only two are cysteine
homologs, both in *F. borealis* and at 25.5 bits, the level of its best random window; none of the cysteine homologs
of the eight other genomes reached 25 bits.

<!-- Next: leave-one-genome-out models for the training genomes; genome-wide scans of non-target tunicates and
shuffled sequences; threshold; the 19 other predictions; the curated genes below 25 bits; F. borealis. -->

# Methods

<!-- Stage 2 subsections are placeholders until stage 2 is run. -->

## Genome assemblies, annotations and transcripts

We analysed three assemblies of *O. dioica*: OKI2018 (RefSeq GCF_907165135.1, with its RefSeq annotation)
[@Bliznina2021], OSKA2016 (genome and transcript models from the ANISEED download archives
`Oidioi_genome_fasta_gff3.zip` and `Oidioi_transcripts_fasta_gff3.zip`) [@Wang2020] and OdB3 (GenBank
GCA_000209535.1, with its GenBank annotation) [@Denoeud2010]. At NCBI, the OKI2018 assembly is named OKI2018_I68;
it is the OKI2018_I69 assembly of Bliznina et al. [-@Bliznina2021], I69 being the identifier of the sequenced
individual, and we use the NCBI name for consistency with the database. The six other appendicularian genomes are
those of *Oikopleura longicauda* (GCA_004367895.1), *Oikopleura albicans* (GCA_004367875.1), *Oikopleura vanhoeffeni*
(GCA_004367855.1), *Bathochordaeus stygius* (GCA_004367955.1), *Mesochordaeus erythrocephalus* (GCA_004367975.1) and
*Fritillaria borealis* (GCA_004368075.1) [@Naville2019]; they have no public gene annotation. As non-target genomes
for SECIS models, we used four other tunicates: *Ciona intestinalis* (GCF_018327825.1), *Ciona savignyi*
(GCA_054131345.1), *Styela clava* (GCF_964204865.1) and *Thalia democratica* (GCF_965202585.1).

Transcript evidence came from one source per *O. dioica* assembly: for OKI2018, the Trinity transcriptome assembled
by Bliznina et al. [-@Bliznina2021] (109,877 sequences; to be added to its Zenodo record [REF: Zenodo record 4604144]); for OSKA2016,
the Trinity assembly of Wang et al. [-@Wang2015] (GEO GSE64421, 86,898 sequences); for OdB3, the cDNA collection of
the genome project (Genoscope, 207,355 sequences: 7,862 clone contigs and 199,493 single reads) [@Denoeud2010],
also distributed by the now offline OikoBase. All files were retrieved by a script that checks their MD5 sums.

## Selenoprotein gene pipeline

We wrote a Nextflow [REF: Nextflow] pipeline, based on the nf-core template [REF: nf-core], that runs the following
steps on each genome; all tools run from container images, and stage 1 used release 0.1.2 with default parameters.

*Selenoprotein genes.* Selenoprofiles 4.5.8 [@Mariotti2010; @Tico2026] searched the profiles of eukaryotic
selenoprotein families and of the Sec machinery (profile sets `eukarya` and `machinery`), with the Selenoprofiles data
set (profiles, Gene Ontology tables and UniRef50 tag database) downloaded from its server; the MD5 sums of the
archives used are recorded. Independently, miniprot 0.11 [REF: miniprot] aligned to the genome the Sec-containing
sequences of the same profiles. A TGA codon was called a Sec codon when miniprot aligned it to a query Sec, or
when the query Sec was the first residue after the end of the alignment and the next genomic codon was TGA, since
miniprot tends to stop before a C-terminal Sec. Sec codons found by both methods at the same position were merged.

*SECIS elements.* Infernal 1.1.5 [REF: Infernal] `cmsearch` searched, with the Rfam SECIS model RF00031
[REF: Rfam], the 6,000 nucleotides downstream of each Selenoprofiles prediction, at the Rfam gathering threshold.
A SECIS element counted as support for a prediction, not as a filter.

*tRNA-Sec.* Secmarker 1.0.0 [@Santesmasses2017] and `cmsearch` with the Rfam tRNA-Sec model RF01852 (gathering
threshold) searched each genome. Predictions overlapping on the same strand were merged into loci; a locus was
considered supported when found by both tools and, where Secmarker reported it, carrying the discriminator base G73.

*Annotation check.* For the genomes with an annotation, each Sec codon predicted by Selenoprofiles was classified
according to the annotated coding sequences overlapping it, following the classes of Ticó et al. [-@Tico2026]: well
annotated (Sec codon inside the CDS, in frame, with or without Sec declared), stop codon (CDS ending at the Sec
codon), skipped (CDS on both sides of the codon but not over it), upstream or downstream (CDS entirely on one side),
out of frame, spliced (CDS boundary inside the codon) or absent. The same predictions were also assessed with
`selenoprofiles assess`, identifying annotated transcripts by the `Parent` of their CDS features.

*Transcript support.* Transcripts were aligned to their genome with minimap2 2.30 [REF: minimap2] (`-x splice
--secondary=no`). A transcript supported a Sec codon when the three nucleotides of the codon lay inside one aligned
exon of its primary alignment; for each codon we recorded the number of supporting transcripts and the extension of
the longest one 3′ of the codon.

Results of all genomes were summarised with MultiQC [REF: MultiQC].

## Baseline curation

<!-- To complete: reference proteins, Sec context, criteria; tables in curation/ of the analysis repository. -->

The baseline lists, for each of the three *O. dioica* assemblies, the Sec codons of the selenoprotein genes, placed
by aligning curated reference proteins and checking the residues around the Sec, and supported by transcripts and
by comparison between assemblies. Assembly errors at a gene were recorded, not corrected.

## SECIS scores against genomic background

For each prediction, Selenoprofiles writes the sequence downstream of the predicted coding sequence, up to 6,000
nucleotides or the end of the contig. We scored these 3′ sequences with Infernal 1.1.5 `cmsearch` on their top
strand, with the default filters, reporting all hits of at least 0 bits, and kept the best score of each sequence.
Predictions were grouped into those with a Sec codon (label "selenocysteine") and the others (cysteine homologs,
predictions not aligned at the Sec position, homologs of the Sec machinery). As background, we drew from each genome
200 windows of 6,000 nucleotides, at positions chosen uniformly over the contigs of at least that length and on a
random strand, and scored them the same way. For the curated *O. dioica* Sec genes, we scored the 6,000 nucleotides
that follow each Sec codon on its strand.

## Appendicularian SECIS model

Training sequences were the 3′ sequences that follow Selenoprofiles predictions labelled "selenocysteine" in seven
oikopleurid genomes (OKI2018, OdB3, *O. longicauda*, *O. albicans*, *O. vanhoeffeni*, *B. stygius* and
*M. erythrocephalus*); OSKA2016 and *F. borealis* were held out. The background for training was 200 random windows
per training genome, drawn as above with a different random seed from the windows used for evaluation. Round 0
selected, as seeds, the best RF00031 hit of each training sequence scoring at least 15 bits. Each round then aligned
the seeds to the current model with `cmalign` (in the first round, RF00031, from which the consensus secondary
structure is inherited), built a model with `cmbuild` (default options), calibrated it with `cmcalibrate -L 0.2`, and
searched the training sequences and the background with `cmsearch` as above. The next seeds were the best hits of
the training sequences scoring above the best background hit of that round. Rounds stopped when the set of seed
sequences did not change, which happened at round 5; the best background hit was 29.3, 23.7, 21.6, 21.0 and 24.8 bits
in rounds 1 to 5. The final model was built from the last seeds and calibrated with default `cmcalibrate` settings.
The procedure is implemented in `stage2/build_secis_model.py` of the analysis repository.

## Stage 2

<!-- To write after stage 2: leave-one-genome-out evaluation; non-target genomes; threshold; lineage reference
proteins and profile HMMs; re-calling. -->

## Data and code availability

The pipeline is available at <https://github.com/oist/LuscombeU_Selenoproteins_Pipeline> (release 0.1.2 for
stage 1). The analysis, from data download to figures, is at <https://github.com/oist/LuscombeU_Oikopleura_Selenoproteins>.
[Zenodo DOIs to add.]

<!--
Notes for the next pieces
- References to add: review of Sec insertion; Sec biosynthesis (PSTK, SecS); Rfam; Infernal; SECISearch3; miniprot;
  minimap2; Nextflow; nf-core; MultiQC; Zenodo record of the OKI2018 transcriptome.
- SelenoDB 2.0 and the Selenoprofiles expectation table were used for the chordate comparison (Ciona gene lists;
  vertebrates); say so in the Methods.
-->
