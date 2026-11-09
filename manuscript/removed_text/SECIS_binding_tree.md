# SECIS-binding protein tree: text replaced on 2026-11-10

The tree of the SECIS-binding proteins was redesigned on 2026-11-09 (Charles: "Let's fix the tree to simplify the
figure and text", with BLASTP results): deuterostome proteins only, rooted on mollusc and annelid proteins, with more
ascidians, a third amphioxus species, a shark, a frog, a bird, lampreys and hagfishes (commit 600cb0a; Deigo run
f2ea142). The text below described the earlier tree (eukaryote-wide profile sequences, rooted on the non-animal
sequences; commit c7865cd) and was replaced by a shorter text on the new tree.

## What the removed text holds, and how to regenerate it

- The trees of the earlier design: make machinery-phylogeny with the SBP2 options of commit c7865cd (without
  --keep-groups and --outgroup-groups), then the long-branch and domain analyses (results of commit c7865cd in git).
- The position of the *F. borealis* fragment next to insect or nematode proteins in that tree.
- The support of SECISBP2L with the invertebrate proteins in the full-length trees and in the region upstream of the
  L7Ae domain.

## Superseded statements (do not reuse as they stand)

- "The full-length trees thus disagree on the relation of the invertebrate proteins to the two vertebrate paralogs":
  with deuterostomes only and the outgroup of molluscs and annelids, the tunicate proteins are the sister group of
  the vertebrate SECISBP2 and SECISBP2L clade in every full-length analysis (66 to 98%).
- "its relationship to the two vertebrate paralogs is not resolved", "relation to SECISBP2 and SECISBP2L unresolved"
  (Table MACHINERY): resolved as co-orthology.
- The supports of the appendicularian clade with the other tunicates (89%, PMSF 99%), of the L7Ae and upstream
  regions (88%, 94%; 66%, 50%), the AU range p = 0.17 to 0.60 and the alignment sizes (78 sequences, 456 sites; domain
  parts 72 sequences, 276 and 180 columns): values of the earlier set.
- "the data of each part rejected the tree of the other part (p ≤ 0.0002)": now p ≤ 0.003 (SECIS-binding protein).

## Removed text

To place the appendicularians
among tunicates, the trees include ten other tunicates, seven ascidians and three salps, and two amphioxus species
(three for PSTK).

The trees identify the appendicularian SECIS-binding protein as the ortholog of the single protein of the other
tunicates. The proteins of the six oikopleurid genomes formed a clade (100%, 87% after removal of 20% of the sites), which grouped
with the SECIS-binding proteins of the ten other tunicates (89% in the maximum-likelihood tree, 99% with PMSF, 97 to 99%
after removal of 10, 30 or 40% of the sites; after removal of 20%, the *F. borealis* protein joined the oikopleurid
clade, and both grouped with the other tunicates, 48%). The *F. borealis* protein, a predicted fragment of 128 residues, had no stable position (next
to insect or nematode proteins in the other analyses, 18 to 58%). Together with the amphioxus, echinoderm, mollusc and
annelid proteins, the tunicate proteins grouped with vertebrate SECISBP2L in the maximum-likelihood tree (58%), with
PMSF (79%) and after removal of 10 or 20% of the sites (72 and 44%), and vertebrate SECISBP2 fell outside this group.
After removal of 30 or 40% of the sites, vertebrate SECISBP2 and SECISBP2L formed one clade (73 and 67%), with the
invertebrate proteins outside it, the tree expected from a duplication in early vertebrates [@Donovan2009]. The
full-length trees thus disagree on the relation of the invertebrate proteins to the two vertebrate paralogs, and so
do the two parts of the protein (below).
Vertebrate SECISBP2L also supports selenoprotein synthesis: zebrafish embryos without Sbp2 still make a subset of
selenoproteins, and only the loss of both genes abolishes their synthesis [@Kiledjian2022], and selenoprotein synthesis
that depends on SECISBP2L is required for oligodendrocyte differentiation [@Dai2022]. The appendicularian protein is
thus the ortholog of the tunicate protein, closest to SECISBP2L by best hits and in most full-length trees, but its
relationship to the two vertebrate paralogs is not resolved.


of the appendicularian PSTK or SECIS-binding protein with the other tunicates (AU tests, p = 0.17 to 0.60 under the two
models).

For the SECIS-binding protein, the L7Ae
domain placed vertebrate SECISBP2 and SECISBP2L in one clade that excluded the tunicate proteins (88%, 94% with PMSF),
whereas the region upstream of it grouped vertebrate SECISBP2L with the amphioxus and echinoderm proteins (66%; with
PMSF, 50%, with mollusc and annelid proteins in the same clade), as did most full-length trees with all invertebrate
deuterostome proteins.

other were found within groups (appendicularians, other tunicates, vertebrates, protostomes, and the non-animal
SECIS-binding proteins), with one exception:

the trees of commit c7865cd (2026-11-03) for SEPHS and the SECIS-binding proteins, of Deigo run fa286f8
  (2026-11-07, three amphioxus PSTK) for PSTK.

| SECISBP2/SECISBP2L | trees | ortholog of the tunicate protein (89%, C60 99%); relation to SECISBP2 and SECISBP2L unresolved |

The SECIS-binding protein is the ortholog of the single SECIS-binding protein of the other tunicates,
and its best human hit is SECIS-binding protein 2-like (SECISBP2L), as shown below.

The tree of the SECIS-binding proteins included human SECISBP2 and SECISBP2L
and was rooted on the non-animal sequences of the profile.

The alignments had 83 sequences
and 421 sites for PSTK (selected model LG+F+I+R6), 115 and 448 for SEPHS (LG+R6) and 78 and 456 for the
SECIS-binding proteins (LG+I+R5).

columns of either part were left out (one SEPHS and six SECIS-binding proteins), so that both parts were analysed
on the same sequences (PSTK 83 sequences, 298 and 123 columns; SEPHS 114, 226 and 222; SECIS-binding proteins 72, 276
and 180).

each family, the data of each part rejected the tree of the other part (p ≤ 0.0002), but the differences lay in weakly
