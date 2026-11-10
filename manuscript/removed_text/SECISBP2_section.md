# SECISBP2/2L section: text removed on 2026-11-10

Charles (2026-11-10) asked to remove the best-hit sentences on the SECIS-binding protein and SEPHS, to blend the
introduction of the machinery trees into the Methods, to remove the SECIS-binding protein tree paragraph except its
conclusion, to rename the section "Reannotation of the selenocysteine machinery of appendicularians", and to add a
SECISBP2/2L section before SECp43 with Figure SBP2. His replacement paragraph: "Only minor changes are needed to the
annotations of the SECIS-binding protein, with the addition of one GC-AG intron in OSKA2016 and OdB3, plus adjustment
of few other junctions. The tree of the SECIS-binding proteins recovers the duplication of SECISBP2 and SECISBP2L in
early vertebrates [@Donovan2009], and confirms the appendicularian protein as their co-ortholog."

## What the removed text holds that the manuscript no longer says

- Best human hits of the SECIS-binding proteins (closer to SECISBP2L by 2 to 38 bits in the nine appendicularian
  genomes; shorter second predictions of B. stygius and M. erythrocephalus closer to SECISBP2; nine of ten other
  tunicates and two amphioxus closer to SECISBP2L, H. roretzi to SECISBP2 by 1 bit): results/tables/stage2_machinery_rbh.tsv
  and stage1_machinery_rbh.tsv.
- Best hits of the selenocysteine SEPHS (SEPHS1 and SEPHS2 within 13 bits) and the statement that best hits do not
  separate these paralogs outside mammals [@Mariotti2015].
- The supports of the SECIS-binding protein tree: kept in the HTML comment of the SECISBP2/2L section.
- The function of vertebrate SECISBP2L in selenoprotein synthesis [@Kiledjian2022; @Dai2022] (now uncited).
- The sampling sentence of the trees (ten other tunicates, three amphioxus for PSTK, the SECIS-binding protein set):
  now in the Methods (Orthology of the machinery genes, Trees), with the reason for the long-branch analyses.

## Superseded statements (do not reuse as they stand)

- "with the addition of one GC-AG intron in OSKA2016 and OdB3" (replacement paragraph): the GC-AG intron is absent
  only from the OSKA2016 ANISEED model; the OdB3 GenBank model is the gold standard unchanged.
- Methods, "the predictions of the appendicularian profiles in the nine appendicularian genomes; and the proteins
  reconstructed from O. dioica transcripts": the machinery trees use the gold-standard proteins of six genomes and the
  final-run predictions of the others, without B. stygius and M. erythrocephalus, and no transcript proteins.

## Removed text

### best hits

The SECIS-binding
proteins of all nine genomes were closer to human SECISBP2L than to SECISBP2, by 2 to 38 bits (a second, shorter
prediction in *B. stygius* and *M. erythrocephalus* was closer to SECISBP2, by 6.5 and 11.5 bits), as were those of nine
of the ten other tunicates and of the two amphioxus species; the *Halocynthia roretzi* protein was closer to SECISBP2 by
1 bit. The selenocysteine SEPHS hit human SEPHS1 and SEPHS2 within 13 bits of each other; best hits do not separate
these paralogs outside mammals [@Mariotti2015]. Best hits thus assign SecS and eEFSec, but not the SECIS-binding protein
or SEPHS2 (nor PSTK, above).

### tree introduction

For these three families, we built maximum-likelihood trees that include the proteins they must be told apart from: the
related kinase KTI12 for PSTK, the SPS1- and SPS2-type proteins of Mariotti et al. [-@Mariotti2015] for SEPHS, and human
SECISBP2 and SECISBP2L for the SECIS-binding proteins. Like all trees of this study, they take the published phylogeny
of the reference sequences as given and ask where the appendicularian sequences fall in it (Methods). The
appendicularian sequences are the gold-standard proteins of
the three *O. dioica* assemblies and of *O. longicauda*, *O. albicans* and *O. vanhoeffeni*, and the predictions of the
final run for *F. borealis* and for the loci without a gold-standard model (among them the selenocysteine SEPHS2 of the
three other *Oikopleura* species); *B. stygius* and *M. erythrocephalus* were left out. To place the appendicularians
among tunicates, the trees include ten other tunicates, seven ascidians and three salps, and two amphioxus species; the
PSTK tree adds a third amphioxus species, and the tree of the SECIS-binding proteins, restricted to deuterostomes and
to full-length proteins, has eleven other tunicates (ten ascidians and one salp), three amphioxus species,
hemichordates, echinoderms and both vertebrate paralogs, with cyclostome proteins (Methods).
Appendicularian proteins evolve fast, those of *O. dioica* about twice as fast as those of *Ciona* [@Berna2012], which
exposes their position in trees to long-branch attraction [@Bergsten2005]. We therefore reanalysed each tree under a
site-heterogeneous model (LG+C60 with posterior mean site frequencies, PMSF) [@Wang2017], after removing the
fastest-evolving sites, and with approximately unbiased (AU) tests of alternative positions (Table LBA; Figure TREES,
PMSF trees).

### SBP2 tree paragraph

In the tree of the SECIS-binding proteins, rooted on mollusc and annelid proteins (Figure TREES), the oikopleurid
proteins form a clade (100% in the maximum-likelihood tree and with PMSF) that groups with the proteins of the other
tunicates (97%, 99%; 96 to 98% after removal of 10 to 40% of the sites). The tunicate proteins are the sister group of
the vertebrate proteins (91%, 96%; 92 to 98%), as tunicates are of vertebrates, and the tree recovers the other
relationships of the deuterostomes sampled: amphioxus with tunicates and vertebrates (76%, 79%; 80 to 88%), and
hemichordates with echinoderms (100% in every full-length analysis). Vertebrate SECISBP2 and SECISBP2L form one clade
(99 to 100%) that also holds cyclostome proteins: a lamprey protein and its hagfish orthologs are the sister group of
the SECISBP2 of jawed vertebrates (98%, 99%), and two other hagfish proteins that of their SECISBP2L (91%, 91%). The
tree thus recovers the duplication of SECISBP2 and SECISBP2L in early vertebrates [@Donovan2009], and
places the appendicularian protein as the ortholog of the single protein of the other tunicates and a co-ortholog of
both vertebrate paralogs. Vertebrate SECISBP2L also supports selenoprotein synthesis: zebrafish
embryos without Sbp2 still make a subset of selenoproteins, and only the loss of both genes abolishes their synthesis
[@Kiledjian2022], and selenoprotein synthesis that depends on SECISBP2L is required for oligodendrocyte
differentiation [@Dai2022].

### notes

Cyclostome proteins of the SECIS-binding protein tree
  (curation/extra_tree_proteins.faa): lamprey LOC116950782 (Petromyzon) with hagfish LOC137431572 (Myxine) and
  LOC147030032 (Eptatretus), sister to jawed SECISBP2 (entry "vertebrate SECISBP2 with lamprey LOC116950782 ...");
  hagfish LOC137424249 and LOC147024609, sister to jawed SECISBP2L (entry "vertebrate SECISBP2L with hagfish ...").
  Upstream region of the L7Ae domain: tunicate clade sister to all other deuterostome proteins (maximum likelihood) or
  to the two LOC137424249/LOC147024609 hagfish proteins (PMSF, 44%) (entry "tunicate proteins", position).

### figure comment

<!-- Figures SBP2 and SEPHS_b: provisional (2026-11-07), placed here for Charles's review of the machinery genes;
they may be removed or aggregated later. Panels drawn by make machinery-reannotation (figures/reannotation_page.py on
figures/inputs/reannotation_<gene>.tsv), from run 0d19c70; the later runs changed only models of other species, which
the panels do not draw. Figure EEFSEC moved to the SEPSECS and EEFSEC section (2026-11-08). -->

### methods appendicularians

the predictions of the appendicularian profiles in the nine appendicularian genomes; and the proteins reconstructed
from *O. dioica* transcripts. The other tunicates were
