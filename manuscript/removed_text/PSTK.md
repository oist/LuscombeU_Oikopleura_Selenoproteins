# PSTK subsection: text removed on 2026-11-07

Replaced by a shorter version derived from a whole-section replacement paragraph by Charles (see the commit that
created this file). The removed text below is verbatim, from "PSTK has five coding exons" to the end of the domain
trees, before Figure PSTK.

## What the removed paragraphs hold, and how to regenerate it

- Public models of each assembly and why each is wrong (RefSeq OKI2018, ANISEED OSKA2016, GenBank OdB3), and the
  curator's exon chain for OdB3: figures/inputs/reannotation_PSTK.tsv, curation/machinery_overrides.tsv (make
  gold-standard-machinery, make machinery-reannotation-inputs).
- Similarity of the gold-standard proteins to human, Ciona, OKI2018 and OSKA2016 (34, 37, 77, 93%): Figure PSTK B
  (make machinery-reannotation).
- Operon and 3' UTR: transcript ends, next and previous genes (accessions in the comment), distances to the
  trans-splice site of the next gene compared with the operon statistics (section "Selenoprotein genes in operons";
  results/tables/stage2_operons_utr.tsv, make operons).
- Splice sites of the third and fourth introns: equivalent placements and their dinucleotides (Figure PSTK C,
  figures/intron_placements.py), frequency of repeated end bases among the 280 gold-standard introns (make
  intron-content), frequency of TT-AG and AT-AG introns in human [@Soo2026], G(non-T)/AG introns [@Henriet2019].
- Third intron retained or spliced at another acceptor (OSKA2016 transcripts, OdB3 GenBank ESTs; 36 in-frame codons),
  and the limit of the anchor codon of curation/gold_standard.py for alternative isoforms (docs/plan.md, 2026-11-04).
- Pairing of the ends of the fourth intron (z-scores, ranks among 271 introns): results/tables/intron_content.tsv,
  make intron-content, make intron-pairing.
- Intron positions across chordates: make pstk-introns, make pstk-intron-columns (results/tables/
  pstk_intron_columns.tsv); scheme of all species: make pstk-intron-scheme (figures/pstk_intron_scheme.*).
- Transcript support of the Ciona models (Ghost database [@Satou2022]).
- O. albicans PSTK: curation of the last exon (tblastn), inverted repeat, target-site duplication, copies in the
  genome, pairing of the intron ends (results/tables/intron_content.tsv).
- Presence of PSTK in the fragmented assemblies, scored against KTI12 (bit differences), second loci of O. longicauda
  and O. vanhoeffeni (make machinery-fragments, make pstk-regions).
- Reciprocal best hits with human, speed of PSTK evolution against KTI12 and other machinery genes (make
  machinery-rbh, make pstk-rates).
- Trees: maximum likelihood, C20, PMSF, slow-fast, one leaf per species, without KTI12, more eukaryotes, archaea, AU
  tests, domain trees, with their supports (make machinery-phylogeny, machinery-lba, pstk-one-leaf, pstk-no-outgroup,
  pstk-enriched, machinery-domains(-pmsf), machinery-clade-supports; Tables LBA and DOMAINS).

## Superseded statements (do not reuse as they stand)

- "In O. longicauda, a second locus, on a contig of 604 nucleotides, has no intron over residues 12 to 198 but two stop
  codons in frame ... a pseudogene, or contamination": the two stop codons lie in a 66-nt insertion that splices out
  as a GG-AG intron in phase 2 inside the W codon of SMRHFW, at the position, phase and sites of the second O. dioica
  intron (make pstk-intron-scheme, 2026-11-07; docs/plan.md).
- "in O. longicauda, an intron of about 70 nucleotides, in phase 0, five alignment columns before the second O. dioica
  intron; in M. erythrocephalus and B. stygius, introns near residues 23 and 92 of the OdB3 protein, and in both
  species an intron in phase 1 near residue 136": the ~70-nt intron (GT-AG, phase 0) is at the same position in the
  three species, after OdB3 residue 81 (R81|S82), not near 92; M. erythrocephalus has no HSP beyond residue 132 of the
  O. vanhoeffeni protein, so its intron near 136 is not seen; M. erythrocephalus also has an intron after residue 22
  of the O. vanhoeffeni protein (GT-AG, 85 nt, phase 2).
- Placements of the fourth intron ("AG-AA in OKI2018, TT-GT in OSKA2016 and AT-AG in OdB3") and the pairing of its
  ends (z = -6.1 in OSKA2016): with the AG convention of 2026-11-07, the gold standard places it on the AG acceptor in
  the three species (TT-AG, TT-AG, AT-AG, phase 1), and the OSKA2016 ends pair with z = -5.3 (run f16a589).
- Tree supports: from commits c7865cd and 3662128, before B. floridae and B. belcheri PSTK were added to the trees
  (make pstk-trees and the following targets, 2026-11-07).

## Removed text

PSTK has five coding exons in the three *O. dioica* assemblies, the third of them 38 nucleotides long (Figure PSTK A).
No public annotation gives the gene: the RefSeq model of OKI2018 holds the first three exons and reads on into the third
intron, the ANISEED model of OSKA2016 starts in the fourth exon and continues past the stop codon into seven exons of
the next gene, and the GenBank model of OdB3 has the first three exons, another first site for the fourth exon and no
last exon. In OdB3, where the prediction skipped the 38-nucleotide exon and read 25 nucleotides further before the last
intron, and where no EST covers the whole gene, the curator set the exon chain: the first three exons of the GenBank
model, the acceptor of the third intron from the prediction, and the last intron of the ESTs. The gold-standard
proteins, of 241, 246 and 246 residues, are about 100 residues shorter than human PSTK (348 residues) and 75 shorter
than that of *Ciona intestinalis* (320 residues), and they are distant from both (Figure PSTK B): the OdB3 protein has
34 and 37% similar residues with them, against 77 and 93% with the OKI2018 and OSKA2016 proteins.

PSTK is an internal gene of an operon in the three assemblies. Its transcripts start with the spliced leader in OKI2018
and OSKA2016, and so do three OdB3 ESTs in GenBank (FP752340, FP831546, FP830122; the 5′ ends are masked in the EST
collection used for the models), and all end 3 to 11 nucleotides after the stop codon. The start
codon of the next gene on the same strand begins 34 to 48 nucleotides after the stop codon, one nucleotide after an AG
acceptor, a candidate trans-splice site; the transcript ends lie 27 to 35 nucleotides before this acceptor, the spacing
found between the polyadenylation site of a gene and the trans-splice site of the next gene of its operon (median 31
and 32 nucleotides; section "Selenoprotein genes in operons"). The 3′ UTR of PSTK is thus a few nucleotides long. None
of these transcripts ends with a poly(A) tail (adapter sequence follows the end of the OKI2018 transcripts), so that the
polyadenylation site is inferred from where the transcripts end.
<!-- 2026-11-03: transcript ends from the candidate tables (OKI2018 1941206, three transcripts; OSKA2016 1355815, four;
OdB3 777154, two ESTs), sequences of the transcript ends printed on Deigo (Charles), genome after the stop codon from
results/gene_regions/PSTK; next gene: XP_081641984.1 (OKI2018), the downstream part of the ANISEED model
S10.g09825.01 (OSKA2016), CBY23134.1 (OdB3); previous gene on the same strand, coding sequences 134 to 143 nt
apart: XP_081641987.1, S10.g09826.01, CBY23136.1. -->

The third and fourth introns have splice sites outside the intron classes of *O. dioica* in OKI2018 and OSKA2016, and so
has the fourth intron in OdB3 (Figure PSTK A). The transcripts that splice these introns all splice them at the same
junctions (isoforms of one assembled locus in OKI2018 and in OSKA2016, two ESTs in OdB3), and the assemblies of other
individuals of the three species (KUM-M3-7f, AOM-5-5f and Bar2_p4) have the same exon chains and the same bases at the
ends of both introns (Figure PSTK C): these splice sites are not assembly errors. Their names, however, depend on where
the junction is placed. The first bases of each of these introns repeat the first bases of the next exon, or its last
bases the last bases of the previous exon (one to three bases, shaded in Figure PSTK C), so that the intron can be moved
across them without changing the spliced sequence, and each placement reads other dinucleotides; when none of them is an
intron class, the gold-standard script keeps the placement nearest to the alignment, which for the fourth intron reads
AG-AA in OKI2018, TT-GT in OSKA2016 and AT-AG in OdB3. Bases repeated at the two ends of an intron are the rule rather
than the exception: of the 280 gold-standard introns of OKI2018, OSKA2016 and OdB3, 243 have at least one and 160 have
two or three (2.2 on average in the 213 GT-AG introns; `make intron-content`), because exons often end with AG, as
introns do. What is unusual in PSTK is not this ambiguity, but that no placement within it gives an intron class. Placed
on an AG acceptor, the fourth intron is a TT-AG intron in OKI2018 and OSKA2016 and an AT-AG intron in OdB3, in the same
phase, and the third intron a TG-AG, AA-AG and GA-AG intron. All end with AG without starting with GT, as does the
second class of *O. dioica* introns after GT-AG [@Soo2026], but only the GA-AG intron belongs to the nearly canonical
G(non-T)/AG introns, an important minority in oikopleurid genomes that human cells do not splice [@Henriet2019]; TT-AG
and AT-AG introns are minor classes of human introns (0.19 and 0.42% in the compilation shown by @Soo2026). The third
intron is retained in two of the four OSKA2016 transcripts (63 nucleotides, with a stop codon in the frame of the
protein) and in three of the five OdB3 ESTs of GenBank (58 nucleotides, which shift the frame); the two others
(FP831546, FP830122) splice it to an acceptor 108 nucleotides further, which removes in frame 36 codons
(FVEMRTPAHLCKERNALRSHPVSDHSMAKIIENYE) that the OKI2018 and OSKA2016 transcripts keep (the GenBank model CBY23135.1 lacks
them too). No OdB3 EST uses the acceptor of the model, which rests on the similarity to OKI2018 and OSKA2016. These two
ESTs were not among the transcripts of the gold-standard script, which takes only those that hold the anchor codon
of the model (here inside the 36 codons). In OSKA2016 and OdB3, the first and last 25
nucleotides of the fourth intron can pair: the energy of their duplex is lower than in dinucleotide-shuffled introns (z
= -6.1 in both), the pairing that would bring the two splice sites together being among the four strongest of the 271
gold-standard introns of *O. dioica* scored that are not insertions (see the gold-standard section); the OKI2018 intron, whose interior differs, has no such pairing (z = -0.7),
so that it is not needed to splice the intron.

In vertebrates and amphioxus, PSTK has five introns at the same alignment columns and in the same phases (RefSeq models
of human, mouse, the two zebrafish genes and *Branchiostoma lanceolatum*; `make pstk-introns`,
`results/tables/pstk_intron_columns.tsv`), the structure of the gene in the chordate ancestor, whereas none of the six
other tunicates with a RefSeq model has an intron in PSTK (*C. intestinalis*, *C. robusta*, *Phallusia mammillata*,
*Styela clava*, *Botryllus schlosseri* and *Thalia democratica*); the *Ciona* models are supported by transcripts over
their whole length (RNA-seq for *C. intestinalis*, ESTs for *C. robusta* in the Ghost database [@Satou2022]). No
appendicularian intron lies at the column of a chordate intron: the nearest, the second and third *O. dioica* introns,
lie seven and six columns from the second chordate intron, in a region aligned without gaps. The tunicate ancestor
therefore had an intronless PSTK, and the introns of appendicularians were gained in their lineages.

The other appendicularians agree with this (regions of `make pstk-regions`). *O. vanhoeffeni* PSTK is one open reading
frame of 248 codons aligned along the whole *O. dioica* protein, and the PSTK of *O. longicauda*, *M. erythrocephalus*
and *B. stygius* have introns at other positions than the *O. dioica* introns (tblastn with the OdB3 protein, junctions
placed by protein similarity, no transcripts): in *O. longicauda*, an intron of about 70 nucleotides, in phase 0, five
alignment columns before the second *O. dioica* intron; in *M. erythrocephalus* and *B. stygius*, introns near residues
23 and 92 of the OdB3 protein, and in both species an intron in phase 1 near residue 136, which the two giant larvaceans
may share. *O. albicans* PSTK has two introns, of 389 and 2,001 nucleotides, four and seven alignment columns away from
the third and fourth *O. dioica* introns, in another phase for the first: introns gained in its own lineage, as the
intronless ancestral state predicts. The prediction read into the last intron and missed the last exon, which a tblastn
search with the C-terminal region of the *O. vanhoeffeni* protein found 2 kb downstream (E = 10^-5^); the curator added
it to the model (243 residues). Placed by the similarity of the protein to those of *O. dioica* and *O. vanhoeffeni*,
this intron has no placement in the intron classes either, but no transcript tests it. A perfect inverted repeat of 12
nucleotides lies near its two ends (ATGGGATATTAT ... ATAATATCCCAT), and the same six nucleotides (AGAAAT) are found at
both ends: the arrangement of the introns of *F. borealis* that arose from miniature inverted-repeat transposable
elements, with their terminal inverted repeats and target-site duplications [@Henriet2019], and that of a gain by the
insertion of a transposon. Like 20 of the 23 introns of the machinery models of *O. albicans*, it shares sequence with
at least ten other loci of its genome (8,927 loci, over 52% of its length; BLASTN, E at most 10^-5^), so that copies do
not single it out; its two ends, which lie outside an assembly gap of 227 nucleotides, however pair more than those of
shuffled sequences (z = -3.2), about as much as those of the most paired of the 19 other *O. albicans* introns scored
(z = -3.1). <!-- z-scores of this paragraph and of the
PSTK fourth intron: shuffles seeded per intron; O. dioica values from Deigo run 5313acd (OSKA2016 -6.06, OdB3 -6.08,
ranks 7 and 6 of 276, 4 and 3 of the 271 that are not insertions, curation/intron_pairing.py; OKI2018 -0.71), O. albicans values from run 12b9bd9 (2026-11-06): PSTK intron 2 -3.21, against
shuffles of the intron without its 227 N (-3.25 in run 8aac2be, with the N shuffled into the controls); SEPSECS intron
4 -3.11, intron 2 -3.07; 20 O. albicans introns scored of 23 (three longer than 2,500 nt or with a gap at an end).
Earlier runs with a
shared generator: -5.5/-5.8 and -4.1 (13b50c5), -4.55/-6.12 and -3.43 (ccaebc0). -->

PSTK needed more evidence outside *O. dioica*. In *O. albicans* and *O. vanhoeffeni* it was predicted as a gene; in the
fragmented assemblies of *O. longicauda*, *B. stygius* and *M. erythrocephalus*, only as fragments, which Selenoprofiles
discards but which match the PSTK of other appendicularians more closely than its relative KTI12. To test these
fragments, we searched each genome with the PSTK proteins of the other genomes (tblastn) and compared each hit with
KTI12. This supported the loci of *O. longicauda* and *M. erythrocephalus* (38 and 50 bits closer to PSTK than to KTI12)
and the PSTK transcripts of the three *O. dioica* transcriptomes (66 to 84 bits), and left that of *B. stygius*
ambiguous (24 bits), like the divergent PSTK of *O. vanhoeffeni* (13 bits), which the tree places in the PSTK family
(Figure PSTK D); in *F. borealis* it found no PSTK-like sequence at all. PSTK is thus supported in seven of the nine
genomes, ambiguous in *B. stygius*, and possibly missing from the *F. borealis* assembly. Because the gene had no intron
in the tunicate ancestor and has few in appendicularians (above), a complete PSTK can lie on a short contig, and loci of
the fragmented assemblies that align along the protein are acceptable evidence of the gene. We use them for presence
only, and do not count copies. In *O. longicauda*, a second locus, on a contig of 604 nucleotides, has no intron over
residues 12 to 198 but two stop codons in frame, and is 45% identical to the other locus and 37 to 42% to the PSTK of
the other *Oikopleura* species: a pseudogene, or contamination by an unsequenced species. In *O. vanhoeffeni*, the
second locus, predicted as a pseudogene, lies 9.5 kb from the gene on the same contig and is 83% identical to it.

Human proteins do not identify the appendicularian PSTK: the predictions of *O. dioica* and *O. albicans* had no hit in
the human proteome, and that of *O. vanhoeffeni* scored 48 bits, whereas the PSTK predicted with the standard profile in
each of the ten other tunicates was a reciprocal best hit of human PSTK (55 to 108 bits). PSTK itself evolves fast: in
the PSTK tree, which holds both families, the PSTK of two chordates are 2.1 to 2.9 times as distant as their KTI12 (4.2
for *C. intestinalis* and *P. mammillata*; 1.4 for human and mouse, the closest pair), in the kinase domain as in the
C-terminal region, and of 18 selenoprotein and machinery genes with human and *Ciona* reference proteins, PSTK has the
lowest identity between the two (30%; median 49%; `make pstk-rates`). The ratio between *C. intestinalis* and *O.
dioica* (2.4) is that of the other pairs: the appendicularian PSTK combines the fast evolution of the gene with that of
appendicularian proteins [@Berna2012]. The appendicularian PSTK are
PSTK and not KTI12, the only protein homologous to PSTK over its whole length, which diverged from it before the
separation of archaea and eukaryotes [@Sherrer2008]: the two families separated with 100% support in every analysis. The appendicularian proteins (*O.
dioica*, *O. albicans* and *O. vanhoeffeni*, and a second *O. vanhoeffeni* locus predicted as a pseudogene) formed a
clade (100%), apart from the other tunicates, which formed their own clade (97 to 99%). The appendicularian clade
branched as the sister of all other PSTK (86% in the maximum-likelihood tree, 68% under C20, 46% with PMSF, 76 to 92%
after site removal, and 95% without the PSTK of the insect-gut fungi *Smittium* and *Zancudomyces*, Harpellales that use selenocysteine
[@Mariotti2019], which also lie on long branches and branched among the animal proteins, at unresolved positions). With one leaf per appendicularian species (the OdB3 model for
*O. dioica*, and without the *O. vanhoeffeni* pseudogene locus), so that related leaves do not weigh on the topology,
the appendicularian clade kept this position (94%, 92% under C20). Without KTI12, the appendicularian clade grouped
with the PSTK of *Smittium* and *Zancudomyces* (72%), the other long branch of the tree. With 60 more eukaryotic PSTK,
of genera absent from the tree whose PSTK were the closest to the appendicularian proteins in searches of the NCBI
protein database (insects, oomycetes, diatoms, dinoflagellates, green algae and other lineages), the appendicularian
clade grouped with the PSTK of two diatoms, *Nitzschia* and *Fistulifera*, in the maximum-likelihood tree (72%), and
with the Harpellales with PMSF (44%). With seven archaeal PSTK also added, the eukaryotic PSTK, the appendicularian
ones included, formed a clade (99%, 98% with PMSF), within which the appendicularian clade grouped with the PSTK of
the nematode *Gnathostoma spinigerum* and of *Paramicrosporidium saccamoebae*, an endoparasite of amoebae related to
fungi (91%, 88%; Figure PSTK D), whose terminal branches are the longest of the eukaryotic PSTK in these trees. The sister group of the
appendicularian clade thus changed with the sequences sampled. <!-- Run 3662128 (Deigo, 2026-11-04). The first run,
a680eb3, whose eukaryotic set also held a dolphin PSTK (Aethalodelphis, filed under the insects by mistake), placed
the appendicularian clade as the sister of all other PSTK (78%) and, with PMSF, with Gnathostoma and Paramicrosporidium
(78%); with the archaeal PSTK, with these two again (73%, 61%), the eukaryotic PSTK forming a clade (87%, 94%). -->
A long branch next to the
outgroup, or next to another long branch, is the position expected from long-branch attraction [@Bergsten2005], and fast-evolving lineages such as nematodes and
insects branch nearer the base than expected in other trees of the selenocysteine machinery [@Mariotti2019]; a tree grouping the appendicularians with the other
tunicates was not rejected (AU test, p = 0.088). The position of the appendicularian PSTK within the family is
therefore unresolved. Trees of the kinase domain and of the C-terminal region recovered the appendicularian clade (91
to 100% under the selected model and with PMSF); the kinase domain placed it as the sister of all other PSTK (89%, 78%
with PMSF), the C-terminal region next to fly proteins (47%, 61% with PMSF), and neither rejected a grouping with the
other tunicates (AU tests, p = 0.11 to 0.42).

