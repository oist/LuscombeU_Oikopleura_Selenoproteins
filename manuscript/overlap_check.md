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
