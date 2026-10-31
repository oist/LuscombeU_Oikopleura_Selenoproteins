# Mass-spectrometry validation of the selenoproteome: lean plan

Draft of 2026-10-31; nothing started. One pool of adults from the Okinawa culture (embryo pools optional), collected
and extracted together with another *Oikopleura* project, is analysed by shotgun proteomics on the OIST instruments by
a trained user; a targeted run follows only if the first one requires it. Transcriptome profiling of the developmental
time course can come later.

## Goal and evidence

The analysis should show which of the 20 predicted selenoproteins are made in *O. dioica*, and where the UGA codon is
read as selenocysteine. No enrichment and no quantification in the first round.

Each gene gets the highest level reached:

| Level | Observation | What it shows |
| --- | --- | --- |
| 1 | any peptide of the protein | the gene is translated |
| 2 | a peptide encoded after the UGA codon | the UGA is read through, not a stop |
| 3 | the peptide holding the selenocysteine: carbamidomethyl-Sec mass, fragments on both sides of the U, selenium isotope envelope in MS1 | selenocysteine at that position |

The cysteine variant of each selenocysteine peptide is searched at the same time: a match to it, rather than to the U
form, would point to a paralog or a wrong model. The criteria of level 3 are those used to accept selenocysteine
peptides by Guo et al. (2018): a high-confidence spectrum after manual assessment, and a precursor with the selenium
isotope signature.

Priorities: KKSU and TXNL1.b (families new to this study), the selenocysteine form of SEPHS2, then the known families.

## What can be detected

In silico digestion of the gold-standard proteins (`make ms-peptides`, `stage2/ms_peptides.py`;
[`results/tables/stage2_ms_peptides.tsv`](../results/tables/stage2_ms_peptides.tsv), five proteases, up to one missed
cleavage, peptides of 7 to 30 residues):

- trypsin reaches the selenocysteine peptide of 17 of the 20 OKI2018 genes; chymotrypsin gives one for the other three
  (GPX_chr1b, TXNL1.b, KKSU), so the two digests together cover all 20;
- 14 genes have at least two tryptic peptides after the selenocysteine (level 2 without the selenocysteine peptide);
- in SELENOK, SELENOO, SELENOS, TXNRD, TXNL1.b and KKSU the selenocysteine lies within two residues of the C-terminus:
  only level 3 can show read-through, and three of these tryptic peptides need one missed cleavage (SELENOS, TXNRD,
  SELENOO);
- 17 of the 20 proteins are below 30 kDa.

| Gene | Mass (kDa) | Tryptic selenocysteine peptide | Tryptic peptides after the Sec | Chymotrypsin, when trypsin has none |
| --- | --- | --- | --- | --- |
| KKSU | 8.2 | none | 0 | KPGKKSU |
| SELENOW_chr1a | 9.4 | IHVIYCGGUGYRPK | 4 | |
| SELENOK | 10.1 | ASAPNMGVGGGUGR | 0 | |
| SELENOW_pair1 | 10.4 | MGVLVHIEFCGIUNHRPK | 6 | |
| SELENOW_pair2 | 10.7 | VNIIYCGVUDHRPK | 4 | |
| SELENOW_PAR | 13.7 | LDFCSSUSYLPK | 5 | |
| SELENOW_chr2 | 13.7 | VHIEYCGSUGYEPR | 3 | |
| SELENOF | 16.4 | SAQFEVCSUK | 2 | |
| TXNL1.b | 18.9 | none | 0 | IGSPIATTNMSELKKVGUG |
| GPX_chr1b | 20.8 | none | 6 | RGNVTIVVNVATNUGLANDHYPQLQNF |
| GPX_chr1a | 21.0 | VTIVTNVASEUGYAK | 6 | |
| GPX_PAR | 21.2 | UGVTDSEYK | 10 | |
| GPX_chr2 | 21.7 | MSEFLGSPTIVVNVASQUGVTDR | 7 | |
| SELENOS | 21.7 | GPSRGGUG | 0 | |
| SELENOT | 22.2 | LGGQNVLAFAICTSUGYR | 6 | |
| DIO_chr1 | 25.0 | YLGEITVLNFGSCTUPPFMNSIK | 3 | |
| DIO_PAR | 27.7 | SFLDENDGLTVLNFGSCTUPPFMASLK | 5 | |
| SEPHS2 | 42.9 | GUGCKVPQEK | 14 | |
| TXNRD | 65.7 | SGAELKNSGCUG | 0 | |
| SELENOO | 67.8 | TEDLCLKVTUSS | 0 | |

OKI2018 gold-standard proteins (the culture is the Okinawa species); OSKA2016 and OdB3 give the same pattern (17 and
16 tryptic selenocysteine peptides). Long selenocysteine peptides are hard to detect: the GPX3 peptide of 44 residues
was not seen by Guo et al. (2018).

## Animals and sample

The culture can provide at most a few hundred individuals per week.

1. Take mature adults of both sexes from the routine culture in batches, as the other project does; rinse in
   filtered seawater, remove the water, snap-freeze in liquid nitrogen, keep at −80 °C; record the number of animals
   per tube.
2. Optional, to agree with the other project: one to two hours in food-free seawater before freezing, so that gut
   algae contribute fewer proteins.
3. Measure protein per animal on a first batch of about 20 animals; this sets the pool size. Target about 50 µg of
   protein per pool, enough for two digests and a repeat.
4. Make two pools from different weeks: pool A is analysed, pool B is kept as a replicate.
5. Optional: embryo pools, shared with the other project, at a stage chosen from the developmental time course (see
   Transcriptome profiling). Each embryo holds little protein, so a pool may need the embryos of many females; the
   number follows from the protein measured on a first pool.

Selenium supplementation, which raised selenoprotein levels in cultured cells (Guo et al. 2018), is not planned.

## Joint design with the other project

Animals are collected once for both projects; one extraction yields their coenzymes and the proteins of both.

- Split extraction: a cold methanol–chloroform–water extraction separates polar metabolites (aqueous phase, for the
  coenzymes), lipids (organic phase) and proteins (precipitate at the interface). The protein pellet dissolves in the
  digestion buffer. The other project's coenzymes set the conditions (light, temperature, redox), as the protein
  pellet tolerates them. IAS lists advice on metabolite extraction and on protein preparation.
- If the other project needs native proteins, or coenzymes bound to them, the frozen pool is split before extraction.
- One proteomics dataset serves both projects: their protein family is searched in the same raw files, with one
  combined protein database.

## Transcriptome profiling (later, not short term)

The unit's developmental RNA-seq time course gives, at no cost in animals, the expression of each selenoprotein and
machinery gene by stage. It would help to choose the embryo stage for the embryo pools, to rank the genes by their
chance of detection, and to test transcript-level claims of the paper with independent libraries:

- expression (transcripts per million) of the 20 gold-standard selenoprotein genes, of SEPHS2, PSTK, SEPSECS, EEFSEC
  and SECISBP2, and of the other project's genes, by stage, in the Okinawa samples mapped to OKI2018 with the
  gold-standard models added to the annotation;
- SELENOI: the share of reads with and without the segment that removes the TGA codon, by stage (an artefact of one
  cDNA library would not recur across stages and libraries);
- coverage of the 3′ UTRs up to the SECIS elements and the poly(A) sites of the operon genes.

The stages and runs of the time course are to confirm (the unit's data catalogue lists six stages of Okinawa and Osaka
animals: eggs, 32 cells, tailbud, hatchling, heartbeat, tailshift).

## Workflow on the OIST instruments

The OIST Instrumental Analysis Section ([Mass Spectrometry](https://groups.oist.jp/ias/mass-spectrometry)) lists an
Orbitrap Fusion Lumos, a Q-Exactive HF, a Q-Exactive Plus and a timsTOF, a training service for LC-MS operation,
advice on sample preparation and data analysis, and peptide synthesis as standards for targeted proteomics.

1. Digestion of the pool A pellet in solution: reduction with DTT, alkylation with iodoacetamide (selenocysteine
   becomes carbamidomethyl-Sec), then a trypsin and a chymotrypsin digest.
2. One data-dependent run per digest. Template from Guo et al. (2018), on a Q-Exactive HF: 120-min gradient, top-20,
   MS1 at 60,000 resolution, MS2 at 30,000. High-pH fractionation only if the proteins are not seen.
3. If proteins are seen but selenocysteine peptides are not: a targeted run (PRM) on the peptides of the table above,
   analysed in Skyline. Synthetic peptides serve as references: Guo et al. (2018) had selenocysteine peptides made
   commercially; their cysteine versions are an alternative, since Sec and Cys peptide pairs have similar fragment
   spectra, retention times and ion mobilities (Si et al. 2025).
4. Only if 1 to 2 mg of protein can be pooled: selenocysteine-selective labelling at low pH and enrichment (Guo et al.
   2018 used 2 mg of protein per sample, labelling at pH 4.0; Bak et al. 2018, pH 5.75).

## Database and search

We search the raw files with Free software.

- Protein database: the OKI2018 RefSeq proteome (GCF_907165135.1); the gold-standard selenoproteins of all three
  assemblies, with U, and the same with C at the selenocysteine position (species-specific peptides also confirm that
  the culture is the Okinawa species); the other project's proteins; common contaminants; decoys made by the search
  tool.
- Search: Comet, which accepts U as a user-defined residue
  ([add_U_user_amino_acid](https://comet-ms.sourceforge.net/parameters/parameters_201201/add_U_user_amino_acid.php),
  no default mass; set to the selenocysteine residue, 150.95364 Da, with ⁸⁰Se). Fixed modification: carbamidomethyl on C
  and U. Variable: Met oxidation, and selenocysteine to dehydroalanine (−81.932170 Da, as in Guo et al. 2018).
  Precursor isotope errors allowed: in the selenium envelope the peak taken as monoisotopic can be the wrong one (Guo et
  al. 2018 allowed one ¹³C error in Mascot for this). Percolator for the false discovery rate (1% at peptide
  level).
- Level 3 by eye for every selenocysteine peptide: fragments on both sides of the U, and an MS1 envelope matching one
  selenium atom (⁷⁴Se 0.89%, ⁷⁶Se 9.37%, ⁷⁷Se 7.63%, ⁷⁸Se 23.77%, ⁸⁰Se 49.61%, ⁸²Se 8.73%; Guo et al. 2018, Figure
  S1), computed for each peptide.
- Raw files and results deposited with the paper (PRIDE).

## Decision points

1. Protein per animal: sets the pool size and whether a 1 to 2 mg pool for enrichment (workflow step 4) is within reach.
2. After the shotgun runs: selenocysteine peptides seen → report per gene; proteins seen without them → targeted run;
   proteins not seen → high-pH fractionation, or pool B with more animals, or stop.

## Questions to settle

At the first consultation with IAS:

- [ ] Who runs the samples: a trained OIST user, or training for one of us?
- [ ] Which instrument (Q-Exactive HF, Orbitrap Fusion Lumos or timsTOF), and how much instrument time.
- [ ] Can IAS synthesise selenocysteine peptides, or only their cysteine versions?
- [ ] Which extraction suits both the coenzymes and the proteins?

With the culture:

- [ ] Protein per animal, measured on about 20 animals, and per embryo.
- [ ] Embryo stage for the embryo pools, after the transcriptome profiling.

With the other project:

- [ ] Their protein family, coenzymes, stage and number of animals.
- [ ] Their extraction constraints (light, temperature, redox; native proteins or not).
- [ ] Data sharing and authorship.

## Sources

- OIST Instrumental Analysis Section, [Mass Spectrometry](https://groups.oist.jp/ias/mass-spectrometry): services,
  sample guideline, instruments (page as saved on 2026-10-31).
- Guo L et al. 2018. Selenocysteine-specific mass spectrometry reveals tissue-distinct selenoproteomes and candidate
  selenoproteins. Cell Chem Biol 25: 1380–1388.
  [doi:10.1016/j.chembiol.2018.08.006](https://doi.org/10.1016/j.chembiol.2018.08.006) (with STAR Methods and
  Supplemental Information).
- Si C et al. 2025. DeepSecMS advances DIA-based selenoproteome profiling through Cys-to-Sec proxy training. Adv Sci
  12: e04109. [doi:10.1002/advs.202504109](https://doi.org/10.1002/advs.202504109)
- Bak DW et al. 2018. A quantitative chemoproteomic platform to monitor selenocysteine reactivity within a complex
  proteome. Cell Chem Biol 25: 1157–1167.
  [doi:10.1016/j.chembiol.2018.05.017](https://doi.org/10.1016/j.chembiol.2018.05.017)
- Comet parameter
  [add_U_user_amino_acid](https://comet-ms.sourceforge.net/parameters/parameters_201201/add_U_user_amino_acid.php).
