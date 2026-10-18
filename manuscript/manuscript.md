---
title: "The selenoproteome of *Oikopleura dioica* and other appendicularians"
bibliography: references.bib
---

<!-- Draft. Citations use Pandoc keys from references.bib; [REF: …] marks references still to add. -->

# Introduction

<!-- Revised after the comments of 2026-10-16. [REF: …] marks references still to add. -->

Selenoproteins carry selenocysteine (Sec), an amino acid encoded by UGA, which is otherwise a stop codon. Selenium is
an essential trace element, with a narrow range between deficiency and toxicity [REF: selenium requirement and
toxicity]. Selenocysteine is a more reactive analogue of cysteine, and most selenoproteins are oxidoreductases in
which it forms the catalytic site: glutathione peroxidases (GPX) and thioredoxin reductases (TXNRD) maintain redox
homeostasis, and iodothyronine deiodinases (DIO) activate and inactivate thyroid hormones [REF: review of
selenoprotein functions]. Selenoprotein sets differ between lineages, which have lost selenoprotein genes or replaced
their selenocysteine by cysteine [@Mariotti2012; @Jiang2012]. Selenocysteine is inserted by a dedicated machinery. In
eukaryotes, a stem–loop structure in the 3′ untranslated region of selenoprotein mRNAs, the selenocysteine insertion
sequence (SECIS) element, is bound by SECIS-binding protein 2 (SECISBP2), which recruits the selenocysteine-specific
elongation factor eEFSec (EEFSEC) and the selenocysteine-charged tRNA-Sec to the UGA codon [REF: review of eukaryotic
selenocysteine insertion]. Selenocysteine is synthesised on its tRNA: seryl-tRNA-Sec is phosphorylated by
O-phosphoseryl-tRNA-Sec kinase (PSTK), then converted to selenocysteinyl-tRNA-Sec by selenocysteine synthase (SecS,
encoded by SEPSECS), using selenophosphate made by selenophosphate synthetase 2 (SEPHS2) [REF: selenocysteine
biosynthesis on tRNA]. This machinery is expected in any genome that encodes selenoproteins.

Because gene prediction programs read UGA as a stop codon, selenoprotein genes are often truncated, split or missing
in genome annotations [@Mariotti2015]. A gene missing from an annotation is invisible to the analyses that start from
annotated genes, such as proteomic searches, transcriptome quantification and orthology inference. In vertebrate
genomes, most selenoprotein genes are misannotated in public databases unless a dedicated pipeline is used, and even
then about one in five remains wrong [@Tico2026]. Outside vertebrates, and for genomes without such a pipeline, the
situation is largely unmeasured.

Programs dedicated to selenoproteins find their genes either by homology to known selenoprotein families, such as
Selenoprofiles [@Mariotti2010; @Tico2026], or from SECIS elements, such as Seblastian, which predicts selenoprotein
genes upstream of the SECIS elements found by SECISearch3 [@Mariotti2013]; tRNA-Sec genes have their own detectors
[@Santesmasses2017]. These methods were built from well-studied lineages, and the SECIS models in particular from
vertebrate and other model genomes: the seed of the Rfam SECIS family, for example, derives from a survey of metazoan
SECIS elements and contains no tunicate sequence [@Lambert2002; REF: Rfam]. How they perform on divergent genomes is
not known in advance.

The appendicularian *Oikopleura dioica* is a planktonic tunicate found in the oceans worldwide [@Masunaga2020].
Appendicularians are among the most abundant zooplankton, and the filter-feeding houses that they discard contribute
to the vertical flux of carbon in the ocean [@Masunaga2022; REF: primary sources, for example Alldredge 2005 and
Robison et al. 2005]. Meta-omic studies now aim at the metabolism of whole plankton communities [@Stec2017;
@Heal2021], and depend on the annotated genomes of their members. The short life cycle, transparent body and
invariant embryonic cell lineages of *O. dioica* make it a model for developmental biology [@Wang2020], and its small
and compact genome, which has lost many genes conserved in other chordates, makes it a model for evolution by gene
loss [@Denoeud2010; @Albalat2016]. *O. dioica* has even lost its peroxisomes [@Zarsky2015; @Kienle2016;
@Bliznina2021], which raises the question of how it detoxifies reactive oxygen species, a function in which several
selenoproteins take part.

*O. dioica* comprises at least three cryptic species, not yet formally described: animals found in the main islands
of Japan and in Okinawa do not interbreed [@Masunaga2022], and the genomes of the three species are extensively
rearranged relative to each other [@Plessy2024]. A genome assembly is public for each: OdB3 from Norway
[@Denoeud2010], OSKA2016 from near Osaka [@Wang2020] and OKI2018 from Okinawa [@Bliznina2021]. How the annotations of
the three assemblies represent selenoprotein genes has not been examined.

Tunicate selenoproteins are known mostly from ascidians, where they show a dynamic evolution. *Ciona intestinalis*
has 19 selenoprotein genes [@Jiang2010] and a selenocysteine-containing deiodinase was characterised in *Halocynthia
roretzi* [@Shepherdley2004]. A survey of the thyroid hormone toolkit across chordates found that tunicate deiodinase
sequences in public databases often lack parts of the protein [@Esposito2021], an example of how annotation errors
limit comparative studies. No appendicularian selenoproteome has been described; to our knowledge, the only
appendicularian selenoprotein gene reported is the selenophosphate synthetase SEPHS2 of *O. dioica*, with its SECIS
element [@Mariotti2015]. <!-- TODO: confirm by a full-text screen of the large surveys (SelenoDB 2.0, Secmarker,
Baclaocos 2019, Trenz 2021, Esposito 2021). -->

*O. dioica* adds difficulties to the annotation of these genes. Like other tunicates, it transcribes many genes as
operons, whose transcripts are resolved into individual mRNAs by trans-splicing [@Ganot2004; REF: operons in
ascidians]; annotation pipelines can merge the genes of an operon into a single gene model [REF: example]. In
addition, part of its introns use non-canonical GA–AG splice sites, which gene predictors model poorly
[@Denoeud2010]. Selenoprotein genes can also take unexpected forms: in ascidians, the selenophosphate synthetase gene
produces an alternative isoform that does not carry selenocysteine [@Mariotti2015].

<!-- TODO (once the results are complete): why the selenoproteome of O. dioica matters, for example whether a lineage
that lost many genes kept the full selenocysteine system or reduced it, as insects and nematodes did; and the
selenoproteome size of marine organisms (Lobanov et al. 2007). -->

<!-- To revise once the results are complete: list the concrete contributions. -->
Here we establish the selenoproteome of *O. dioica*, and we use the genomes of six other appendicularian species,
together with the selenoproteomes of other tunicates, to extrapolate a baseline for appendicularians. We built a
pipeline that combines the detection of selenoprotein genes, of their selenocysteine machinery and of SECIS elements
with annotation and transcript evidence. With their default parameters, current tools did not recover this baseline
in appendicularian genomes; we provide lineage-specific models and settings that do, and the data to evaluate them.
We hope that these data will help build defaults that apply to more lineages.

# Results

<!-- Figures to be numbered. -->

## A pipeline for selenoprotein genes, their machinery and SECIS elements

We wrote a Nextflow pipeline that brings together, for each genome, the evidence needed to annotate selenoprotein
genes (Figure PIPELINE, tube map). It predicts selenoprotein genes and the genes of the selenocysteine machinery with
Selenoprofiles; aligns independently the selenocysteine-containing sequences of the same profiles with miniprot and
calls the selenocysteine codons that they align to; searches SECIS elements downstream of each prediction; and
searches tRNA-Sec genes with two methods. When a gene annotation is given, it classifies each predicted
selenocysteine codon according to the annotated coding sequence that overlaps it, and when transcripts are given, it
counts those that support the codon. A MultiQC report summarises the results of all genomes. By default, the pipeline
uses the settings of the underlying tools: the eukaryotic and machinery profile sets of Selenoprofiles for eukaryotes,
the Rfam SECIS model at its gathering threshold, and the Rfam tRNA-Sec model together with Secmarker.

The pipeline's test runs every process on three small genomes. On human chromosome 22, it predicted selenoprotein O
(SELENOO) with its selenocysteine codon, found by both Selenoprofiles and miniprot, supported by a transcript and
correctly annotated in a test annotation; in the yeast genome, which has no selenoprotein, it predicted none.
Chromosome 22 carries no tRNA-Sec gene. The SECIS element of human SELENOO scored only 16 bits with the Rfam model,
below its gathering threshold of 35 bits; the test therefore lowers the threshold to 15 bits so that the SECIS search
is exercised.

## Default settings miss part of the selenoproteome of *O. dioica*

We ran the pipeline with its default settings on the three *O. dioica* assemblies and the genomes of six other
appendicularians, and compared its predictions with a baseline of 18 selenoprotein genes curated in each
*O. dioica* assembly (Methods). The pipeline found 13 of the 18 selenocysteine codons in OKI2018, 14 in OSKA2016 and 13
in OdB3, nearly all through Selenoprofiles; miniprot found 6, 6 and 3. Transcripts supported almost all the codons
found (13 of 13, 14 of 14 and 12 of 13). Two genes were missed in all three assemblies: selenoprotein K (SELENOK), and
SEPHS2, which Selenoprofiles predicted but without calling its selenocysteine codon. Other genes were missed in one or
two assemblies, for reasons that include assembly gaps. No SECIS element reached the gathering threshold in any of the
nine genomes (see below), whereas tRNA-Sec was found, by both methods, in every genome.

The public annotations represent these genes poorly. Of the selenocysteine codons predicted by Selenoprofiles, 5 of 12
fall in a correctly annotated coding sequence in the RefSeq annotation of OKI2018, none of 14 in the OSKA2016
annotation and none of 21 in the OdB3 annotation; in the other cases, the annotated coding sequence ends at the codon,
skips it, lies on one side of it or is absent.

With the same settings, the pipeline performed as expected on vertebrates, from which its tools were mostly built. It
predicted 30 selenoprotein genes in human and in mouse and 47 in zebrafish, found the complete machinery in all
three, and 27 of 34, 24 of 36 and 47 of 50 of the selenocysteine codons that it checked fell in correctly annotated
coding sequences of their RefSeq annotations. In other tunicates, it predicted 16 to 25 selenoprotein genes per
genome, including 19 in *Ciona robusta*, the number reported for this genome [@Jiang2010], and the machinery genes,
including PSTK and SECp43. In the three *Ciona* species, SEPHS2 was predicted with a residue other than
selenocysteine, consistent with its alternative isoform that does not carry selenocysteine [@Mariotti2015]. The Rfam
SECIS model was strict even in vertebrates: it supported 12 of the 30 predictions in human and in mouse and 10 of 47
in zebrafish, against 0 to 2 per genome in the other tunicates and none in appendicularians. Results on the
appendicularians were identical between two releases of the pipeline (0.1.2 and 0.2.0).

## Appendicularian profiles recover the selenocysteine machinery

With the standard profiles, Selenoprofiles did not find PSTK in any of the three *O. dioica* assemblies, although its
transcripts are expressed, and did not call the selenocysteine codon of SEPHS2. We therefore extended each profile of
the machinery set with appendicularian proteins: the predictions of the standard profiles in the appendicularian
genomes, and the PSTK and SEPHS2 proteins reconstructed from *O. dioica* transcripts (Methods). The OdB3 assembly
contributed no sequence and served as a test. With these profiles, the machinery was found in every *O. dioica*
assembly, including PSTK, whose predicted gene models in OKI2018 and OdB3 match the curated ones, and SEPHS2 with its
selenocysteine codon (Table MACHINERY). In the fragmented assemblies of *O. longicauda*, *B. stygius* and *M.
erythrocephalus*, PSTK was predicted only as fragments, which Selenoprofiles discards but which match the PSTK of
other appendicularians more closely than its relative KTI12. A search of each genome with the PSTK proteins of the
other genomes (tblastn), followed by a comparison with KTI12, supported the loci of *O. longicauda* and *M.
erythrocephalus* (38 and 50 bits closer to PSTK than to KTI12) and the PSTK transcripts of the three *O. dioica*
transcriptomes (66 to 84 bits), and left that of *B. stygius* ambiguous (24 bits), like the divergent PSTK of *O.
vanhoeffeni* (13 bits); in *F. borealis* it found no PSTK-like sequence at all, so that the gene may be missing from
that assembly. SEPHS2 was predicted with its selenocysteine codon in seven of the nine genomes; it was not found in
*B. stygius* and only partially in *M. erythrocephalus*. SecS (SEPSECS), eEFSec (EEFSEC), SECIS-binding protein 2
(SECISBP2) and tRNA-Sec were found in every genome. Selenoprofiles also searches for SECp43 (TRNAU1AP), which binds
tRNA-Sec in complexes with SecS and SEPHS2 [REF: SECp43] but is an associated factor whose role is debated:
selenoprotein synthesis did not depend on it in mouse liver and neurons [@Mahdi2015], whereas it enhanced the
translation of several selenoproteins in glioblastoma stem cells [@Zhang2026]. We therefore report it separately.
SECp43 was found only in *F. borealis*, whereas the standard profile found it in all other tunicates and in the
vertebrates. The same search, applied to the eight oikopleurid genomes and the three *O. dioica* transcriptomes,
found only proteins at most 6 bits closer to SECp43 than to TIA1-like RNA-binding proteins, whereas the SECp43
proteins of the other genomes were 65 to 434 bits closer, and it recovered the SECp43 of *F. borealis* (Methods).
SECp43 is therefore absent from the Oikopleuridae, and was most parsimoniously lost in their ancestor.

| Gene | *O. dioica* OKI2018 | *O. dioica* OSKA2016 | *O. dioica* OdB3 | *O. longicauda* | *O. albicans* | *O. vanhoeffeni* | *B. stygius* | *M. erythrocephalus* | *F. borealis* |
|---|---|---|---|---|---|---|---|---|---|
| tRNA-Sec | found | found | found | found (4 loci) | found (2 loci) | found (2 loci) | found | found (2 loci) | found |
| SEPHS2 | selenocysteine | selenocysteine | selenocysteine | selenocysteine | selenocysteine | selenocysteine | not found | partial | selenocysteine |
| PSTK | found | found | found | fragment | found | found | fragment | fragment | not found |
| SEPSECS | found | found | found | found | found | found | found | found | found |
| EEFSEC | found | found | found | found | found | found | found | found | found |
| SECISBP2 | found | found | found | found | found | found | found | found | found |
| TRNAU1AP (associated factor) | not found | not found | not found | not found | not found | not found | not found | not found | found |

<!-- Table generated by stage2/machinery_census.py (results/tables/stage2_machinery_census.md). -->

Table MACHINERY. Selenocysteine machinery, and the associated factor SECp43, found with the appendicularian profiles. "Selenocysteine": SEPHS2 predicted
with its selenocysteine codon; "found": predicted; "partial": predicted only with an in-frame stop codon or frameshift,
or without alignment at the selenocysteine position; "fragment": discarded by Selenoprofiles but supported by
similarity (Methods); for tRNA-Sec, loci found by both methods. Several copies in an assembly may reflect duplicated
haplotypes or contamination [@Pichon2019].

## The Rfam SECIS model does not separate appendicularian SECIS elements from genomic background

In stage 1, no SECIS element reached the Rfam gathering threshold (35 bits) downstream of any selenoprotein gene
predicted in the nine genomes. To see whether the threshold alone was too strict, we scored with the same model,
without threshold, the 3′ sequences of all Selenoprofiles predictions and, as background, 200 random genomic windows
of the same length per genome (Figure SECIS-1). Predictions fall into two groups: 143 with a selenocysteine codon (below, selenocysteine
predictions), and 578 others. The others are mostly homologs with cysteine at the selenocysteine position (375, of which 235 in the EhSEP2 family),
predictions whose alignment does not cover the selenocysteine position of their profile (122, labelled "unaligned"), homologs of
selenocysteine machinery proteins that are not selenoproteins (40, such as eEFSec and SECIS-binding protein 2), and predictions
with another amino acid or a stop codon at the selenocysteine position (41).

Of the 143 3′ sequences that follow predictions with a selenocysteine codon, 96 (67%) had no hit at all, even when every hit
of at least 0 bits was reported. This was not caused by fragmented assemblies separating coding sequences from their UTRs: in the
chromosome-scale assembly OKI2018, 5 of 12 had no hit. The remaining 47 scored from 12 to
31 bits; 18 reached 20 bits, and none reached the threshold. Random windows rarely had a hit (45 of 1,800), but those
hits reached 22.6 bits. Lowering the threshold to 20 bits would therefore have recovered at most four SECIS elements
per genome while accepting background hits, and most selenoprotein 3′ sequences would still have no hit. The Rfam
model does not recognise most appendicularian SECIS elements, rather than recognising them below its threshold.

## An appendicularian SECIS model built from the stage 1 predictions

We built a new covariance model from the stage 1 outputs alone, by iterative bootstrap (Methods). Training used the
110 3′ sequences that follow predictions with a selenocysteine codon in seven oikopleurid genomes: OKI2018 and OdB3 for
*O. dioica*, *O. longicauda*, *O. albicans*, *O. vanhoeffeni*, *B. stygius* and *M. erythrocephalus*. The OSKA2016
assembly of *O. dioica* and the genome of *F. borealis*, a fritillariid, were kept out of training as controls. The
initial seeds were the 26 training sequences whose best Rfam hit reached 15 bits. After five rounds, in which each
new model was used to collect the training hits that score above every background window, the set of seeds no
longer changed: 73 hits, in 73 of the 110 training sequences. The seeds lie close to the predicted stop codons (median
152 nucleotides downstream; quartiles 74 and 388) and are 60 to 75 nucleotides long.

We evaluated the final model on new random windows and on the curated *O. dioica* selenoprotein genes (Figure SECIS-2). Of
1,800 random windows, 8 had a hit of at least 20 bits and 1 of at least 25 bits (27.3 bits, in *F. borealis*); in
the three *O. dioica* assemblies, no random window exceeded 20.8 bits. Hits downstream of selenoprotein genes reached 62 bits.
In the held-out OSKA2016 assembly, 12 of the 14 selenocysteine predictions and 15 of the 18 curated selenoprotein genes had a hit of at
least 25 bits, against a best background window of 18.1 bits. In the training genomes, 12 of 12 selenocysteine predictions of
OKI2018 and 18 of 21 of OdB3 scored at least 20 bits, but these sequences were used for training and overstate the
performance. In the held-out *F. borealis*, only 8 of 19 selenocysteine predictions reached 20 bits, and 3 reached 25 bits.

Across the three *O. dioica* assemblies, 42 of 53 curated selenoprotein genes scored at least 25 bits. Below 25 bits were
selenoprotein T (SELENOT) in all three assemblies, SEPHS2 in OKI2018 and OdB3, selenoprotein O (SELENOO) in OKI2018 and OSKA2016, selenoprotein K (SELENOK) in OKI2018 and OdB3 (23.4 and
24.9 bits), and selenoprotein F (SELENOF) in OSKA2016 and OdB3. For SEPHS2 and SELENOT, the selenocysteine codon lies far upstream of the stop
codon, so that the 6,000 nucleotides scored after the selenocysteine codon may not reach the SECIS element. Outside the selenocysteine
predictions, 19 of the 578 other Selenoprofiles predictions had a hit of at least 25 bits. Seventeen of them belong to
selenoprotein families and either do not cover the selenocysteine position ("unaligned": GPX, DIO, SEPHS2, SELENOT and others) or
carry another residue there (selenoprotein W (SELENOW) with isoleucine, tyrosine or leucine), as expected if some are selenoprotein genes that
Selenoprofiles did not recognise as such, including SEPHS2 of OSKA2016 (labelled "pseudo"). Only two are cysteine
homologs, both in *F. borealis* and at 25.5 bits, the level of its best random window; none of the cysteine homologs
of the eight other genomes reached 25 bits.

<!-- Next: leave-one-genome-out models for the training genomes; genome-wide scans of non-target tunicates and
shuffled sequences; threshold; the 19 other predictions; the curated genes below 25 bits; F. borealis. -->

# Methods

<!-- Stage 2 subsections are placeholders until stage 2 is run. -->

## Genome assemblies, annotations and transcripts

Genes are named by the symbols of their human orthologs (HGNC), given at their first mention (the human tRNA-Sec gene
is TRU-TCA1-1); appendicularian
orthologs take the same symbols, and selenoprotein families the root of their members' symbols (for example GPX for
glutathione peroxidases).

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
(GCA_054131345.1), *Styela clava* (GCF_964204865.1) and *Thalia democratica* (GCF_965202585.1). These tunicates,
with *Ciona robusta* (KH assembly, GCF_000224145.3, the genome of the published *Ciona* selenoproteome, formerly named
*C. intestinalis*) and *Halocynthia roretzi* (GCA_013436055.1), were also analysed with the default settings, to compare the
appendicularians with other tunicates. As positive controls for the default settings, we used the human
(GRCh38.p14, GCF_000001405.40), mouse (GRCm39, GCF_000001635.27) and zebrafish (GRCz11, GCF_000002035.6) genomes, with
their RefSeq annotations.

Transcript evidence came from one source per *O. dioica* assembly: for OKI2018, the Trinity transcriptome assembled
by Bliznina et al. [-@Bliznina2021] (109,877 sequences; to be added to its Zenodo record [REF: Zenodo record 4604144]); for OSKA2016,
the Trinity assembly of Wang et al. [-@Wang2015] (GEO GSE64421, 86,898 sequences); for OdB3, the cDNA collection of
the genome project (Genoscope, 207,355 sequences: 7,862 clone contigs and 199,493 single reads) [@Denoeud2010],
also distributed by the now offline OikoBase. All files were retrieved by a script that checks their MD5 sums.

## Selenoprotein gene pipeline

We wrote a Nextflow [REF: Nextflow] pipeline, based on the nf-core template [REF: nf-core], that runs the following
steps on each genome; all tools run from container images. Stage 1 analysed all genomes (appendicularians,
other tunicates and vertebrate controls) in one run of release 0.2.0 with default parameters.

*Selenoprotein genes.* Selenoprofiles 4.5.8 [@Mariotti2010; @Tico2026] searched the profiles of eukaryotic
selenoprotein families and of the selenocysteine machinery (profile sets `eukarya` and `machinery`), with the Selenoprofiles data
set (profiles, Gene Ontology tables and UniRef50 tag database) downloaded from its server; the MD5 sums of the
archives used are recorded. Independently, miniprot 0.11 [REF: miniprot] aligned to the genome the selenocysteine-containing
sequences of the same profiles. A TGA codon was called a selenocysteine codon when miniprot aligned it to a query selenocysteine, or
when the query selenocysteine was the first residue after the end of the alignment and the next genomic codon was TGA, since
miniprot tends to stop before a C-terminal selenocysteine. Selenocysteine codons found by both methods at the same position were merged.

*SECIS elements.* Infernal 1.1.5 [REF: Infernal] `cmsearch` searched, with the Rfam SECIS model RF00031
[REF: Rfam], the 6,000 nucleotides downstream of each Selenoprofiles prediction, at the Rfam gathering threshold.
A SECIS element counted as support for a prediction, not as a filter.

*tRNA-Sec.* Secmarker 1.0.0 [@Santesmasses2017] and `cmsearch` with the Rfam tRNA-Sec model RF01852 (gathering
threshold) searched each genome. Predictions overlapping on the same strand were merged into loci; a locus was
considered supported when found by both tools and, where Secmarker reported it, carrying the discriminator base G73.

*Annotation check.* For the genomes with an annotation, each selenocysteine codon predicted by Selenoprofiles was classified
according to the annotated coding sequences overlapping it, following the classes of Ticó et al. [-@Tico2026]: well
annotated (selenocysteine codon inside the CDS, in frame, with or without the selenocysteine declared), stop codon (CDS ending at the selenocysteine
codon), skipped (CDS on both sides of the codon but not over it), upstream or downstream (CDS entirely on one side),
out of frame, spliced (CDS boundary inside the codon) or absent. The same predictions were also assessed with
`selenoprofiles assess`, identifying annotated transcripts by the `Parent` of their CDS features.

*Transcript support.* Transcripts were aligned to their genome with minimap2 2.30 [REF: minimap2] (`-x splice
--secondary=no`). A transcript supported a selenocysteine codon when the three nucleotides of the codon lay inside one aligned
exon of its primary alignment; for each codon we recorded the number of supporting transcripts and the extension of
the longest one 3′ of the codon.

Results of all genomes were summarised with MultiQC [REF: MultiQC]. The pipeline's test profile runs every process
on human chromosome 22 (with an annotation of SELENOO and a transcript covering its selenocysteine codon), the yeast
genome and a fragment of the *Escherichia coli* genome, with the SECIS threshold lowered to 15 bits; modules and
pipeline are also tested with nf-test [REF: nf-test].

## Baseline curation

<!-- To complete: reference proteins, selenocysteine context, criteria; tables in curation/ of the analysis repository. -->

The baseline lists, for each of the three *O. dioica* assemblies, the selenocysteine codons of the selenoprotein genes, placed
by aligning curated reference proteins and checking the residues around the selenocysteine, and supported by transcripts and
by comparison between assemblies. Assembly errors at a gene were recorded, not corrected.

## Appendicularian profiles of the selenocysteine machinery

Each Selenoprofiles profile of the machinery set (SEPHS2, PSTK, SecS, eEFSec, SECIS-binding protein 2 and SECp43) was
extended with appendicularian proteins: the predictions of the standard profiles in the appendicularian genomes other
than OdB3, labelled "homologue" or "selenocysteine" and at least half as long as the median sequence of the profile,
and the proteins reconstructed from *O. dioica* transcripts (PSTK from the OKI2018 and OSKA2016 transcriptomes, SEPHS2
with its selenocysteine from the OSKA2016 transcriptome). The added proteins were aligned to the standard profile with
MAFFT `--add --keeplength` [REF: MAFFT], which keeps the columns of the profile and therefore its configuration, and
the extended alignments were built into profiles with `selenoprofiles build`, using all their sequences as BLAST
queries. The profiles keep the standard names; given to the pipeline with `--selenoprofiles_overrides`, they replace
the standard profiles of the same name.

## Census of the selenocysteine machinery

We ran Selenoprofiles with these profiles on the nine genomes. For each genome without a PSTK prediction, we
reran Selenoprofiles with the PSTK profile, keeping the predictions that its last filter discards, and compared
each of them (blastp [REF: BLAST]) with the PSTK proteins found in the other genomes and with KTI12, a related
kinase. This test needs proteins of a related family to compare with, and was not applied to the other genes. A discarded prediction counted as a fragment of the gene when its best hit was a protein of
the gene with an E-value of at most 10^-10^. tRNA-Sec loci were those of stage 1 found by both methods.

## Absence tests

For a machinery gene or associated factor missing from appendicularian genomes, the stage 1 predictions of the gene
in the other genomes were searched in each appendicularian genome and in the *O. dioica* transcriptomes with tblastn
(E-value at most 10^-3^). Hits less than 20 kb apart on the same strand formed a locus, whose translated hits
(genomes) or whose open reading frame around the best hit (transcripts) gave a candidate protein. Each candidate was
compared by blastp with the proteins of the gene (the profile sequences and the stage 1 predictions of the other
genomes) and with proteins of related families: TIA1-like RNA-binding proteins from the Selenoprofiles UniRef50
database for SECp43, KTI12 for PSTK. The margin, the score of the best hit to the gene minus that of the best hit to
the related family, classified the candidate as supported (at least 30 bits), ambiguous (10 to 30 bits) or a related
protein. The stage 1 predictions of the gene, scored in the same way, served as positive controls.

## SECIS scores against genomic background

For each prediction, Selenoprofiles writes the sequence downstream of the predicted coding sequence, up to 6,000
nucleotides or the end of the contig. We scored these 3′ sequences with Infernal 1.1.5 `cmsearch` on their top
strand, with the default filters, reporting all hits of at least 0 bits, and kept the best score of each sequence.
Predictions were grouped into those with a selenocysteine codon (label "selenocysteine") and the others (cysteine homologs,
predictions not aligned at the selenocysteine position, homologs of the selenocysteine machinery). As background, we drew from each genome
200 windows of 6,000 nucleotides, at positions chosen uniformly over the contigs of at least that length and on a
random strand, and scored them the same way. For the curated *O. dioica* selenoprotein genes, we scored the 6,000 nucleotides
that follow each selenocysteine codon on its strand.

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
- References to add: review of selenocysteine insertion; selenocysteine biosynthesis (PSTK, SecS); Rfam; Infernal; SECISearch3; miniprot;
  minimap2; Nextflow; nf-core; MultiQC; Zenodo record of the OKI2018 transcriptome.
- SelenoDB 2.0 and the Selenoprofiles expectation table were used for the chordate comparison (Ciona gene lists;
  vertebrates); say so in the Methods.
-->
