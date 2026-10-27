# Overlap check (OpenAlex, 2026-10-16)

Queries: `openalex_queries.sh` (outside the repository). Works citing the appendicularian genome papers (Denoeud et
al. 2010, Naville et al. 2019, Wang et al. 2020, Bliznina et al. 2021, Plessy et al. 2024, through its 2023 preprint) that mention
selenoproteins; works citing Lambert et al. 2002, Selenoprofiles and Secmarker that mention tunicates; topic searches
since 2015 (selenoproteins with *Oikopleura*, appendicularians, larvaceans, tunicates; SECIS in tunicates;
chordate selenoproteomes).

No selenoprotein census of *Oikopleura* or of other appendicularians was found. Published work touching the same
genes or families, to cite and to compare with:

| Work | Relevance |
|---|---|
| Mariotti et al. 2015 (`Mariotti2015`) | selenophosphate synthetase evolution across the tree of life. Reports *O. dioica* (OdB3) as having a single SPS gene, the ancestral Sec-containing SPS2 with a SECIS and the ancestral intron structure, without the Gly isoform found in ascidians: our SEPHS2 is therefore a known Sec gene. Also reports stem–loops (SRE) next to the UGA of UGA-containing SPS genes |
| Esposito et al. 2021 (`Esposito2021`) | thyroid hormone toolkit (TPO, DIO, THR) in chordates, from UniProt and ANISEED proteins. Reports a single DIO homolog per non-vertebrate chordate (none in one species), and the absence of the catalytic globular domain and its RPLYXXFGS motif from all tunicate DIOs, which they take as a possible sign of different catalytic sites while calling for better structural annotation. Species names are missing from the full-text extraction used here: whether *O. dioica* is among the species, and with how many DIOs, is to be checked in the article's tables and figures. If it is, our two Sec-containing *O. dioica* DIO genes, and the annotation errors we document, bear directly on both statements |
| Darras 2021 (`Darras2021`) | review of deiodinase evolution: deiodinases in all chordates; in tunicates, cites the Sec-containing ascidian enzyme of *Halocynthia roretzi* and partial ascidian sequences suggesting at least two deiodinase homologs. No statement about *Oikopleura* |
| Shepherdley et al. 2004 (`Shepherdley2004`) | first tunicate deiodinase characterised (*H. roretzi*): Sec-dependent outer-ring activity, SECIS in the 3′ UTR |
| Jiang et al. 2012 (`Jiang2012`) | selenoproteomes of six invertebrates (sponge, placozoan, sea anemone, limpet, annelid, amphioxus); no tunicate |
| Trenz et al. 2021 (`Trenz2021`) | glutathione peroxidase evolution; Sec and Cys forms |
| Baclaocos et al. 2019 (`Baclaocos2019`) | SELENOP evolution in metazoans |
| Mariotti et al. 2012 (`Mariotti2012`) | vertebrate selenoproteomes, for the chordate comparison |
| Frese et al. 2024 (`Frese2024`) | *Ciona* proteome across embryogenesis (mass spectrometry; possible check of Sec peptides) |

| Ticó et al. 2026 (`Tico2026`) | annotation classes used here; for vertebrates, RefSeq (with a dedicated selenoprotein pipeline) is correct for most selenoprotein genes, Ensembl and GenBank for few: a contrast with the *O. dioica* RefSeq annotation |

Not yet identified: Sreelatha et al. 2018 (SELENOO/SelO as an AMPylase; the title used in the query was wrong).
Still to check: the species list of Esposito et al. 2021 (see above).

## Second round (OpenAlex, 2026-10-19)

735 distinct works: works citing Mariotti et al. 2015 (tunicate topics), Jiang et al. 2010, Mahdi et al. 2015,
Donovan and Copeland 2009, Zarsky and Tachezy 2015 (tunicate topics), Berna et al. 2012 and 2015, Lobanov et al.
2007 (marine topics) and Kocot et al. 2018, and 18 topic searches since 2015 (peroxisomes and redox in Oikopleura,
machinery genes, selenium in plankton, misannotation, fast evolution, operons, long branches). No title or abstract
mentions appendicularians together with selenium, selenocysteine or a machinery gene: no prior report of
appendicularian selenoproteins, of SECp43 loss in Oikopleuridae or of the second appendicularian SEPHS. Works found
and useful: Manta et al. 2022 (review of the selenophosphate synthetase family), Dai et al. 2022 (SECISBP2L drives
selenoprotein synthesis in oligodendrocytes), Berna et al. 2014 (fast-evolving tunicates), Ferrandez-Roldan et al.
2019 (O. dioica and gene loss), Delsuc et al. 2018 (tunicate phylogenomic framework), Szantho et al. 2023
(compositionally constrained sites drive long-branch attraction). No published example of a selenoprotein missed in
an omics analysis because of misannotation was found.

## Third round (OpenAlex, 2026-10-27)

1,316 distinct works: the works behind the remaining [REF] placeholders, prior work for steps 3 to 5 (SECIS-based
discovery of selenoproteins, SECIS structure and its variants, TXNL1 and its PITH domain), works citing Walczak et
al. 1996, Grundner-Culemann et al. 1999, Rivas et al. 2017, Lambert et al. 2002, Castellano et al. 2004 and Andersen
et al. 2009 (on topics), 19 topic searches since 2015, and works on the core topics since 2026-10-01. No title or
abstract mentions appendicularians together with selenium, selenocysteine, SECIS or a machinery gene: still no prior
report of appendicularian selenoproteins, of a selenocysteine in TXNL1 or of the KKSU family. Metadata of the cited
works were taken from PubMed.

- Placeholders filled: C20/C60 profiles, Le et al. 2008 (`Le2008`, Bioinformatics; the title search had returned
  another paper of the same authors); slow-fast, Brinkmann and Philippe 1999; AU test, Shimodaira 2002. Still open:
  an example of operon genes merged into one gene model (no such example found).
- SECIS-based discovery and comparative prediction, cited in step 4: Kryukov et al. 1999, Lescure et al. 1999,
  Castellano et al. 2001; SGP2 (Parra et al. 2003) and its use for TGA-containing genes (Castellano et al. 2004).
- SECIS structure, cited in the model description and the Discussion: Krol 2002 (review of SECIS forms), Latrèche et
  al. 2009 (kink-turn; internal loop and helix 2 modulate recoding activity), Chapple et al. 2009 (SECISaln), Gao et
  al. 2023 (covariation artefacts, including structural alignments of non-homologous sequences, their issue 11),
  Korotkov et al. 2002 (SELENOM SECIS with cytidines in place of the unpaired adenosines) and Novoselov et al. 2007
  (GGGA in place of AUGA in apicomplexans).
- C-terminal selenoproteins: Haft and Tolstoy 2025 (prokaryotic families with the selenocysteine near one end were
  missed by earlier searches), cited with KKSU.
- Reviews and surveys by Ticó and Mariotti: the metazoan selenoproteome (Annu Rev Anim Biosci 14, 2026) and the
  vertebrate selenoproteome (BMC Genomics 27, 2026), cited in the Introduction.
- SECIS-binding proteins: Kiledjian et al. 2022 (zebrafish Secisbp2l supports selenoprotein synthesis without Sbp2),
  cited with Dai et al. 2022 in the orthology section.
- TXNL1 and the proteasome: Andersen et al. 2009 (Txnl1 binds Rpn11 of the 26S proteasome); Gao et al. 2025 (cryo-EM
  of TXNL1 on the 19S particle, via PSMD1, PSMD4 and PSMD14; proteasome binding needed for the
  ubiquitin-independent degradation of TXNL1 after exposure to metal- or metalloid-containing oxidants); Arkinson et
  al. 2025 (in the degrading proteasome, the C-terminal tail of TXNL1 covers the catalytic groove of Rpn11 and
  coordinates its active-site zinc); Lee et al. 2026 (PITH domain on RPN11, thioredoxin domain towards RPN2 and RPN13;
  His289, the last residue, coordinates the RPN11 zinc; a 19S state without visible TXNL1 C-terminal tail; PDB entries of X. Chen et al.,
  "Human 19S proteasome bound to TXNL1 PITH domain [...]", probably from this study). Lead for the Discussion: in *O. dioica*, the
  selenocysteine lies at the position of human Gly283 and is followed by one residue, where the human tail that
  reaches the RPN11 zinc has six; the *O. dioica* protein also lacks the thioredoxin domain. Also found: a 2026
  review of TXNL1 (Wu et al., Front Cell Dev Biol, doi:10.3389/fcell.2026.1876585), not cited.
- Not found by title: the "ancient origin" SECIS paper of Chapple et al. and the PITH domain paper of Scott et al.
  2005.
