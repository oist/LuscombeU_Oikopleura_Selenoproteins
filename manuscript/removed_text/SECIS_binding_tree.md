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


# SECIS-binding protein tree: text replaced on 2026-11-10, set without truncated proteins

Charles (2026-11-10): "There are a few sequences that are obviously truncated, judging from their length. They may be
a distraction for the reader. Can we remove them?" The truncated proteins were left out (curation/tree_exclusions.tsv;
two partial lamprey proteins deleted from curation/extra_tree_proteins.faa), full-length proteins added from BLASTP
searches restricted by taxon (commit db16781; Deigo run a98ee1f, jobs 52954560-2). The text below described the trees of
Deigo run f2ea142 and was replaced by the values of the new run.

## What the removed text holds, and how to regenerate it

- The position of the *F. borealis* fragment next to the oikopleurid clade (49%, 71%): run make sbp2-trees at commit
  07a6b97 (results of commit f2ea142 in git).
- The second lamprey protein (Petromyzon XP_032837153.2, 565 residues) as the sister group of jawed SECISBP2L (72%,
  89%), and the partial *L. planeri* protein (XP_078478270.1) with zebrafish SECISBP2.
- The region upstream of the L7Ae domain grouping the tunicate proteins with vertebrate SECISBP2 and four cyclostome
  proteins (83%, 74%), with 43 sequences.

## Superseded statements (do not reuse as they stand)

- Supports of the oikopleurid clade (94%, 88%), of the tunicate clade (64%, 93%; 88 to 97%) and of tunicates with
  vertebrates (66%, 97%; 89 to 98%): now 100%, 100%; 97%, 99% (96 to 98%); 91%, 96% (92 to 98%).
- "two other hagfish proteins have no stable position": they are now the sister group of jawed SECISBP2L (91%, 91%).
- "neither part rejected the full-length tree (AU tests, p = 0.10 to 0.36)": the region upstream of the L7Ae domain now
  rejects it (p = 0.046; 0.017 and 0.024 under PMSF).
- AU range for the grouping with the other tunicates in the domain trees, p = 0.17 to 0.63: now 0.17 to 0.60.
- Alignment sizes 51 sequences and 881 sites (JTT+R5), domain parts 43 sequences, 563 and 318 columns with eight
  proteins left out: now 989 sites (JTT+F+R5), 49 sequences, 623 and 366 columns, two left out.
- Table MACHINERY row "(64%, C60 93%) ... (66%, C60 97%)": now (97%, C60 99%) and (91%, C60 96%).

## Removed text

PSTK tree adds a third amphioxus species, and the tree of the SECIS-binding proteins, restricted to deuterostomes, five
more ascidians, a third amphioxus species and vertebrate proteins, cyclostomes included (Methods).

In the tree of the SECIS-binding proteins, rooted on mollusc and annelid proteins (Figure TREES), the oikopleurid
proteins form a clade (94% in the maximum-likelihood tree, 88% with PMSF), which the *F. borealis* fragment joins (49%,
71%), and the appendicularian proteins group with those of the other tunicates (64%, 93%; 88 to 97% after removal of
10 to 40% of the sites). The tunicate proteins are the sister group of the vertebrate proteins (66%, 97%; 89 to 98%),
as tunicates are of vertebrates, and vertebrate SECISBP2 and SECISBP2L form one clade (99 to 100% in every
full-length analysis) that also holds cyclostome proteins: a lamprey protein and its hagfish orthologs are the sister
group of the SECISBP2 of jawed vertebrates (98%, 97%), and a second lamprey protein that of their SECISBP2L (72%, 89%);
a partial *Lampetra planeri* protein groups with zebrafish SECISBP2, and two other hagfish proteins have no stable
position. The tree thus recovers

Neither half rejected a grouping
of the appendicularian PSTK or SECIS-binding protein with the other tunicates (AU tests, p = 0.17 to 0.63 under the two
models).

For the SECIS-binding protein, the L7Ae
domain gave the full-length tree (tunicate proteins as the sister group of the vertebrate proteins, 78%, 96% with
PMSF), whereas in the region upstream of it the tunicate proteins grouped with vertebrate SECISBP2 and four cyclostome
proteins (83%, 74%); neither part rejected the full-length tree (AU tests, p = 0.10 to 0.36).

were found within groups (other tunicates, vertebrates, protostomes and amphioxus), with one exception

  PSTK, and of Deigo run f2ea142 (2026-11-10, deuterostome set rooted on molluscs and annelids) for the SECIS-binding
  proteins. Cyclostome proteins of the SECIS-binding protein tree (curation/extra_tree_proteins.faa): lamprey
  LOC116950782 (Petromyzon) with hagfish LOC137431572 (Myxine) and LOC147030032 (Eptatretus), sister to jawed SECISBP2;
  Petromyzon XP_032837153.2, sister to jawed SECISBP2L; L. planeri XP_078478270.1 (LOW QUALITY, 490 residues, 56%
  identical to zebrafish SECISBP2), with zebrafish SECISBP2 (100%); hagfish LOC137424249 and LOC147024609, sister to
  the SECISBP2 group (45%) or to all vertebrate proteins (PMSF, 48%).

and proteins from BLASTP searches at NCBI, one per gene and species (`curation/extra_tree_proteins.faa`): five more
ascidians (*Corella eumyota*, *Ascidiella aspersa*, *Aplidium turbinatum*, *Ascidia mentula* and *Clavelina
lepadiformis*), *B. belcheri*, the SECISBP2L of a shark (*Mustelus asterias*), a frog (*Xenopus laevis*) and a bird
(*Gavia stellata*), three lamprey proteins (two of *Petromyzon marinus*, one of *Lampetra planeri*) and two proteins of
each of two hagfishes (*Myxine glutinosa* and *Eptatretus burgeri*).

realigned. Appendicularian, curated and other tunicate sequences, the subject of the trees, were kept even when
partial.

and 421 sites for PSTK (selected model LG+F+I+R6), 115 and 448 for SEPHS (LG+R6) and 51 and 881 for the
SECIS-binding proteins (JTT+R5).

columns of either part were left out (one SEPHS and eight SECIS-binding proteins), so that both parts were analysed
on the same sequences (PSTK 83 sequences, 298 and 123 columns; SEPHS 114, 226 and 222; SECIS-binding proteins 43, 563
and 318).

| SECISBP2/2L | trees | ortholog of the tunicate protein (64%, C60 93%), which is the sister group of vertebrate SECISBP2 and SECISBP2L (66%, C60 97%) | found | found | found | found | found | found | found | found | found |
