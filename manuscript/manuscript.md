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

## Stage 2

<!-- To write after stage 2: SECIS covariance model (seeds, held-out tests, non-target genomes, thresholds,
covariation); lineage reference proteins and profile HMMs; re-calling. -->

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
