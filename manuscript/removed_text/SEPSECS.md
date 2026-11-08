# SEPSECS subsection: text removed on 2026-11-08

Replaced by a two-sentence version derived from a replacement paragraph by Charles (see the commit that created this
file). The removed text below is verbatim, from the provisional comment of 2026-11-05 to the end of the subsection,
before Figure SEPSECS.

## What the removed paragraphs hold, and how to regenerate it

- Positions and phases of the seven shared introns and of the two introns of the Okinawa species, on a MAFFT L-INS-i
  alignment with the human and *Ciona* proteins, including the assemblies of other individuals (KUM-M3-7f, AOM-5-5f,
  Bar2_p4) and the other oikopleurids: make intron-map (curation/intron_map.py, results/tables/intron_map_SEPSECS.tsv).
- Gain of the two Okinawa introns after the separation of the three *O. dioica* species: present in OKI2018 and
  KUM-M3-7f; one exon in OSKA2016, OdB3, AOM-5-5f, Bar2_p4 and in the OSKA2016 transcripts and OdB3 ESTs that cover
  it; no intron at these columns in *O. longicauda*, *O. albicans* or *O. vanhoeffeni* (make intron-map; transcripts
  in figures/inputs/reannotation_SEPSECS.tsv, make machinery-reannotation-inputs). The Discussion still relies on this
  gain; the evidence is kept there as a comment.
- Transcript support: TRINITY_DN6371_c0_g1_i1 splices both Okinawa introns; TRINITY_DN6371_c0_g1_i4 retains the
  second, with a stop codon in frame (Figure SEPSECS).
- Pairing of the ends of the first Okinawa intron, z = -4.0 against dinucleotide shuffles (results/tables/
  intron_content.tsv, make intron-content; GA-AG introns, unchanged by the AG convention of 2026-11-07).
- Comparison with the other oikopleurids: one shared intron at the column and phase of an *O. albicans* intron (near
  human residue 259; GA-AG in *O. dioica*, GT-AG in *O. albicans*), none at an intron of *O. longicauda* or *O.
  vanhoeffeni*, one intron shared by *O. albicans* and *O. vanhoeffeni* (make intron-map).
- Donor of one shared intron: GC in the Okinawa species, GA in the other two (results/tables/intron_map_SEPSECS.tsv).
- First intron of 423 and 437 nucleotides in OdB3 and Bar2_p4, against 40 to 45 in the four other assemblies, part of
  it repeated at 50 other loci of OdB3 (copies counted by make intron-content, results/tables/intron_content.tsv).

## Superseded statements (do not reuse as they stand)

None.

## Removed text

<!-- Provisional (2026-11-05), for Charles's review: keep or drop. Intron map: make intron-map
(curation/intron_map.py, results/tables/intron_map_SEPSECS.tsv), on the gold-standard models of run 14314dd. -->

The three *O. dioica* species share seven SEPSECS introns, at the same alignment columns and in the same phases in
their six assemblies, and the Okinawa species has two more (Figure SEPSECS; `make intron-map`,
`results/tables/intron_map_SEPSECS.tsv`).
These two introns, of 86 and 75 nucleotides, near residues 63 and 86 of the human protein, are both GA-AG. They are in
OKI2018 and in KUM-M3-7f, another individual of the species, and the OKI2018 transcript TRINITY_DN6371_c0_g1_i1
splices both (another transcript, TRINITY_DN6371_c0_g1_i4, retains the second, which puts a stop codon in frame). In
OSKA2016 and OdB3 the same stretch is one exon, in the assemblies, in the six OSKA2016 transcripts and the OdB3 ESTs
that cover it, and in AOM-5-5f and Bar2_p4. *O. longicauda*, *O. albicans* and *O. vanhoeffeni* have no intron at
these columns (the nearest, an *O. longicauda* intron, lies two columns from the second, in a region aligned without
gaps). The two introns were therefore most likely gained in the Okinawa lineage after it separated from the other two
species. The ends of the first pair more than those of shuffled introns (z = -4.0).

The seven shared introns also differ from those of the other oikopleurids: only one lies at the column and in the
phase of an *O. albicans* intron (near human residue 259; GA-AG in *O. dioica*, GT-AG in *O. albicans*), none at an
intron of *O. longicauda* or *O. vanhoeffeni*, and *O. albicans* and *O. vanhoeffeni* share one other. Within
*O. dioica*, the donor of one shared intron is GC in the Okinawa species and GA in the other two, and the first intron
is 423 and 437 nucleotides long in the North Atlantic assemblies (OdB3, Bar2_p4), against 40 to 45 in the four others,
part of it being repeated at 50 other loci of OdB3.

## Version of 2026-11-08 (first), replaced on 2026-11-08

Merged with EEFSEC into "SEPSECS and EEFSEC" from a second replacement paragraph by Charles. Nothing superseded;
the Okinawa intron counts and sources are kept in the comment of the new subsection. Verbatim:

### SEPSECS

<!-- Section rewritten on 2026-11-08 from a replacement paragraph by Charles; the removed text, with an index of the
extra analyses behind it: manuscript/removed_text/SEPSECS.md. Sources: Table MACHINERY (make machinery-table: reciprocal
best hit in all nine genomes, no other human gene hit); make intron-map (results/tables/intron_map_SEPSECS.tsv): seven
introns at the same alignment columns and in the same phases in the six O. dioica assemblies; the two introns of the
Okinawa species (86 and 75 nt, near human residues 63 and 86) in OKI2018 and KUM-M3-7f, spliced by the OKI2018
transcript TRINITY_DN6371_c0_g1_i1. -->

SEPSECS, which loads selenium on the tRNA, is conserved enough that reciprocal best hits with human proteins identify
it without ambiguity (Table MACHINERY). The three *O. dioica* species share seven introns, and the Okinawa species has
two more, both GA-AG (Figure SEPSECS).

