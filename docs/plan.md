# Plan of the analysis

## Stage 1 — published tools and parameters (justifies stage 2)

For the nine appendicularian genomes:
- Selenoprofiles 4.5.8 with the default metazoan profiles; Secmarker for tRNA-Sec.
- Generic protein survey: miniprot alignments of profile sequences, Sec codon = in-frame TGA matching the reference
  Sec context.
- SECIS: Rfam RF00031 (Infernal) on 3′ regions of the Sec genes.
- Public annotations of the three *O. dioica* assemblies compared with the Sec codons (Ticó et al. 2026 classes).
- Transcript evidence: OKI2018 transcriptome, OSKA2016 transcripts, OdB3 ESTs (splice-aware alignment; Sec codon in
  an exon; screens for Sec-extended and C-terminal Sec reading frames with a SECIS).
Outputs that justify stage 2: genes missed by every tool; RF00031 sensitivity on appendicularian 3′ regions;
annotation error classes.

## Curation — baseline from the three *O. dioica* assemblies

Tables in `curation/`: one row per Sec gene and assembly, with the Sec codon, the evidence (tool, transcript, cross-
assembly comparison) and notes on assembly errors (recorded, not fixed).

## Stage 2 — final parameters (delivers the results)

- *Oikopleura* SECIS model built from baseline 3′ regions (Oikopleuridae seeds, OSKA2016 held out as control);
  held-out, shuffled, non-target (other tunicates) and covariation tests; threshold choice.
- Lineage reference proteins (curated *O. dioica* sequences) and profile HMMs for genes missed in stage 1
  (SELENOK, PSTK), each tested with the target genome left out.
- Re-calling on the nine genomes with the final parameters; census per genome; machinery; Cys homologs.
- Operons in the three *O. dioica* assemblies (collinear genes at most 500 nt apart; Plessy et al. 2024).

## Figures

Built from `results/tables/` only.

## Notes for the manuscript

- State that `OKI2018_I68` (NCBI GCF_907165135.1) is the `OKI2018_I69` assembly of Bliznina et al. 2021 (I69: the
  sequenced individual; I68 is a typo in the NCBI record), and that the NCBI name is used for consistency.

## Decisions and open items (2026-10-16)

- Order of the results: stage 1; Sec machinery; complete pool of Sec genes (annotated from proteins and transcripts);
  SECIS model trained on that pool (with same-locus and 90% identity redundancy filters, family balance, and
  leave-one-family-out evaluation).
- tRNA-Sec: standard tools are enough (Secmarker and RF01852 find at least one supported locus in every genome); no
  special-casing.
- Multiple copies of a gene (tRNA-Sec loci, eEFSec, SEPHS2) in some appendicularian assemblies may come from
  incomplete haplotype purging or from contamination: the B. stygius assembly may be contaminated by M.
  erythrocephalus (Pichon et al. 2019, F1000Research, doi:10.12688/f1000research.21551.2; cite where copy numbers or
  B. stygius results are discussed).
- Fragmented assemblies: the pipeline may miss genes whose locus is broken by the assembly (current long-read
  assemblies do not have this problem); the analysis shows such genes as fragments (for example PSTK in
  O. longicauda, B. stygius and M. erythrocephalus) to validate the baseline. At the end, review the assemblies and
  drop those of too low quality that add little phylogenetically.
- Open: which genes form the machinery set (minimum: tRNA-Sec, SEPHS2, PSTK, SecS, eEFSec, SBP2; possible: SECp43,
  SCLY).
- SECIS model, when rebuilt from the complete pool (restoring what the earlier analysis behind the slides had):
  1. 90% identity filter on the seeds (and one seed per genomic locus);
  2. leave-one-species-out evaluation;
  3. non-selenoprotein 3' regions and shuffled genomes as backgrounds;
  4. whole-genome and non-target (Ciona, Styela, Thalia) scans, to set the thresholds;
  5. covariation of the base pairs.
- Machinery (in progress): appendicularian Selenoprofiles profiles (<family>_app) built from the standard profiles
  plus appendicularian proteins (stage 1 predictions, curated transcripts), OdB3 held out; tested by Selenoprofiles
  on the nine genomes against stage 1; then deposited (Zenodo) and passed to the pipeline as extra profiles.

## Decisions (2026-10-24)

- **Genome sets.** Stage 2 is developed and trained on the six *Oikopleura* assemblies only: OKI2018, OSKA2016 and
  OdB3 (*O. dioica*), *O. longicauda*, *O. albicans* and *O. vanhoeffeni*. *B. stygius*, *M. erythrocephalus* and
  *F. borealis* are held out and used at the end of the manuscript, after stage 2, for the census over as many
  appendicularians as possible (and as tests of how far the stage 2 tools generalise beyond *Oikopleura*).
  Reasons: their contribution to stage 2 is small (SECIS model: 1 of 73 seeds, from *B. stygius*; machinery profiles:
  8 sequences, 5 from *F. borealis*); their assemblies are fragmented (*B. stygius*, *M. erythrocephalus*: short N50,
  fragment predictions) or noisy (*F. borealis*: 188 selenoprotein predictions, 31 % without a human hit, many
  cysteine homologues grouping with bacteria or fungi); holding out *F. borealis* removes the circularity of the
  SECp43 profile, which contains the *F. borealis* protein it then recovers; and most splice junctions of
  *F. borealis* are non-canonical (Charles Plessy), which the gene prediction tools handle poorly.
  *F. borealis* stays essential for its phylogenetic position (only non-oikopleurid appendicularian): the SECp43
  loss in the oikopleurid ancestor, SELENOL and selenocysteine SELENOU, PSTK absent or missed, and the timing of the
  GPX and SELENOW expansions belong to the final census discussion.
- **No recomputation now.** Existing results stay as they are. The change is applied the next time each step is
  rerun for another reason: SECIS model retrained without *B. stygius* and *M. erythrocephalus*; machinery profiles
  rebuilt and retested without the three genomes; selenoprotein profiles built without them from the start; the
  manuscript statements that rely on them (SECp43 loss, PSTK fragments, their columns in Table MACHINERY, the
  *F. borealis* SECIS evaluation) moved to the census section. Trees: decide at the next rerun whether the three
  genomes stay as descriptive taxon sampling (they are not used to train any tool) or are removed.
- **Bacterial-like sequences: possible horizontal gene transfer, strict standard.** Similarity to bacterial sequences
  is "interesting if true": tunicates acquired cellulose synthase by a lateral gene transfer from a prokaryote at the
  base of the lineage (Sagane et al. 2010, Development, doi:10.1242/dev.044503, PMID 20335363; two CesA genes in
  *O. dioica*), so other transfers cannot be excluded. But no false positive is allowed: a candidate is worth
  reporting only after every attempt to disprove it has failed. Contamination is the default explanation. Checks
  to pass, at least: the gene lies on a contig with animal genes and appendicularian composition (GC, coverage,
  codon usage), not on a short or bacterial-like contig; it has introns, or spliced transcripts from the species;
  it is found in independent assemblies or individuals (for *O. dioica*: OKI2018, OSKA2016, OdB3 and the
  transcriptomes); a foreign-contamination screen of the assembly (for example NCBI FCS-GX, or Tiara) does not flag
  the contig; and its placement among bacteria holds under site-heterogeneous models, with vertical descent rejected
  by AU tests. One HGT screen, at the end of the manuscript, for all genomes used in the census.

## Plan of the next steps (agreed 2026-10-24)

Done: stage 1; machinery census (with one selenoprotein, SEPHS2); hit validation, round 1 (reciprocal best hits and
trees of the selenoprotein families). Analysis order below; the manuscript may present the steps in another order.
Proteomic data that may become available later serve only for the final validation, at the end of the manuscript;
no earlier step waits for them.

1. **SECIS search in the validated hits**, with the existing models only (Rfam RF00031 and the current model);
   descriptive, trains nothing. A hit is left out only on evidence against it: best human hit in an unrelated
   family, placement among bacteria or fungi, or a held-out genome (*B. stygius*, *M. erythrocephalus*,
   *F. borealis*). Lack of evidence is not a reason: family-level validation counts the related human genes (MIEN1
   for SELENOW, SEPHS1 for SEPHS2) and the tree placement in the family, and the curated *O. dioica* genes are always
   kept (SELENOS has no human hit). Checked on the 73 seeds of the current model: the first wording (best hit
   outside the expected genes, or no hit) would have removed 14, of which 13 trusted selenoproteins (8 chr2-type
   SELENOW with best hit MIEN1, including SELENOW_chr2 of OKI2018 and OdB3; 2 *O. vanhoeffeni* SELENOW of the same
   clade with no hit; 3 selenocysteine SEPHS2 with best hit SEPHS1) and 1 *B. stygius* GPX; the rule above removes
   only the latter.
2. **Annotation fix** of OKI2018, OSKA2016 and OdB3 (gold standard): selenocysteine codons, stop codons and 3' UTR
   ends from transcripts; GFF3 with the selenocysteine codons marked and protein sets; later the RefSeq/GenBank
   corrections. OKI2018 SELENOS is annotated as a disrupted allele (the OKI2018 transcriptome carries both alleles,
   so the 13-nt deletion is real; whether this individual carries both is unknown); OdB3 lacks the SELENOO
   C-terminus; OSKA2016 is the only assembly with all 18 genes intact and is the reference for SELENOS. Done before
   step 3 because SECIS windows anchored on the selenocysteine codon may miss the SECIS of SEPHS2 and SELENOT.
3. **New SECIS model** from the six *Oikopleura* assemblies (*O. longicauda* included): 3' regions from the stop
   codons of step 2; one SECIS per locus; 90 % identity filter on the SECIS hits (vsearch or CD-HIT). On the 73
   current seeds, 90 % merges 73 into 61 (OdB3 haplotype copies, the GPX_chr1b SECIS of OKI2018 and OdB3, copies
   within *O. longicauda*, *O. albicans* and *O. vanhoeffeni*, and the OKI2018 SELENOW pair, whose overlapping
   windows may hold one SECIS twice). Evaluation leave-one-family-out and leave-one-species-out (the three
   *O. dioica* assemblies held out together: they share near-identical SECIS). Final comparison, on the same held-out
   sets and backgrounds (random windows, shuffled sequences, whole-genome scans, non-target tunicates), with Rfam
   RF00031, Rfam merged with appendicularian sequences, and the current model, to check that updating an existing
   model is not the better solution.
4. **UGA/SECIS screens** in every genome where possible: transcript-based (OKI2018 RefSeq and Trinity
   transcriptome, OSKA2016 annotation and transcripts, OdB3 annotation and ESTs); genome-based (SECIS hit, then an
   open frame through a TGA upstream) for *O. longicauda*, *O. albicans* and *O. vanhoeffeni*; human and a few
   tunicates with known selenoproteomes as positive controls of the screen.
5. **TXNL1** (thioredoxin-like protein 1, candidate; curation/candidates.tsv), with a rule fixed in advance: counted
   only if the new model scores its SECIS above background in the three *O. dioica* assemblies and the TGA is
   conserved; otherwise reported as a candidate. Also check the TGA position in the other tunicates. Public
   proteomics: none found for appendicularians (PRIDE API, 2026-10-24: no project for Oikopleura, larvacean or
   appendicularian, against 9 for Ciona and 1 for Botryllus; published mass spectrometry on *O. dioica* covers house
   proteins, histone modifications and metallothioneins, with no deposited data found). Proteomic data that may
   become accessible later: final validation only.
6. **Appendicularian selenoprotein profiles** from the gold standard, each tested with an assembly held out; the
   machinery profiles rebuilt without *F. borealis* (decision of 2026-10-24).
7. **Stage 2**, with handling of non-canonical introns (the SELENOT selenocysteine exon is lost at a GA-AG intron):
   the six *Oikopleura* assemblies, then the three held-out genomes. Stage 1 and stage 2 measured against the gold
   standard on held-out assemblies only (the gold standard feeds the stage 2 profiles).
8. **Hit validation, round 2**: AU tests for SELENOW_PAR (SELENOW against SELENOV), SELENOF, SELENOO, SELENOT; human
   MIEN1 in the SELENOW references; Ensembl identifiers assigned to taxonomic groups; held-out genomes in the trees
   or not.
9. **Census of all appendicularians** and comparison with the other tunicates, after a contamination screen
   (FCS-GX or Tiara) and the single strict HGT screen.
10. **SECIS and operons.**
11. **Discussion**; Zenodo deposits (models, profiles, gold-standard annotations); pipeline release with the stage 2
    models; methods written as each step finishes.

## Manuscript style

- Write "selenocysteine" when the amino acid is meant alone; keep "Sec" only inside names (tRNA-Sec, eEFSec, SecS)
  and in the first definition, "selenocysteine (Sec)".
- Give the gene symbol (human HGNC symbol) at the first mention of each gene, and the full name at the first mention
  of a symbol (for example "SECIS-binding protein 2 (SECISBP2)", "selenoprotein T (SELENOT)"); expand every acronym at
  its first use (for example "selenocysteine insertion sequence (SECIS)").
- Orthology of the machinery genes (to establish before the paper; outside the pipeline, in stage2/orthology/):
  1. first pass for every predicted gene: reciprocal best hit against a reference proteome that contains the related
     families (human RefSeq proteome, from the control genomes), recording the best and second-best hits;
  2. for the ambiguous cases of step 1 (decided 2026-10-17: reciprocal best hits first, phylogeny only where they
     are ambiguous), a maximum-likelihood tree (IQ-TREE) of the appendicularian proteins with orthologs from
     other tunicates and vertebrates and with the closest related family as outgroup (PSTK with KTI12; eEFSec with
     eEF1A and other translation GTPases; SEPHS2 with the SPS1-type paralogs of Mariotti et al. 2015; SECISBP2 with
     SECISBP2L; SECp43 with other RRM proteins such as TIA1); fast-evolving appendicularian sequences need models that
     resist long-branch attraction (for example LG+C20 or site removal);
  3. diagnostic features where they exist: selenocysteine or cysteine at the catalytic site of SEPHS2 (SPS2-type),
     the domain IV extension of eEFSec, the domain architecture (Pfam);
  synteny is of little use in O. dioica (scrambled genomes), except between its cryptic species.
- Discussion: a selenoprotein gene can be missed for three reasons that a survey must distinguish: the prediction
  method fails, the annotation pipeline fails, or the assembly is wrong at that locus (moved out of the introduction).
- Introduction TODO: confirm by a full-text screen of the large surveys that SEPHS2 (Mariotti et al. 2015) is the only
  appendicularian selenoprotein gene reported.
- Introduction, to revisit once the results are complete: the first paragraph (machinery; consider a table of the
  machinery genes with their roles and symbols); the importance of the O. dioica selenoproteome (placeholder in the
  text); the last paragraph (concrete contributions).
- Gene symbols: TRU-TCA1-1 (human tRNA-Sec gene) appears once, in the Methods, for text searches.
- Background noted: community metabolomics of marine plankton (Llewellyn et al. 2015 Prog Oceanogr, keywords
  "meta-metabolomics"; Heal et al. 2021 mSystems) and meta-omics models of plankton ecosystems (Stec et al. 2017);
  zooplankton-derived dissolved organic matter linked to zooplankton meta-transcriptomes (bioRxiv
  10.1101/2024.10.27.620474). Our genome annotations of zooplankton feed such approaches.
- Definition (agreed 2026-10-17): "the machinery" is the six components dedicated to selenocysteine and required for
  it: tRNA-Sec, SEPHS2, PSTK, SEPSECS (synthesis), EEFSEC, SECISBP2 (insertion). SECp43 (TRNAU1AP) is reported as an
  associated factor: not required in mouse liver and neurons (Mahdi et al. 2015), homologs in lineages without
  selenocysteine (yeasts, land plants in the Selenoprofiles profile), but reported by Selenoprofiles by default and
  absent from all oikopleurids (discussion). Other borderline genes (SEPHS1, SECISBP2L, SCLY, SARS1, ALKBH8, TRIT1,
  RPL30, nucleolin, EIF4A3) are excluded by the specificity criterion; SECISBP2 versus SECISBP2L is an orthology
  question for the reciprocal best hits.
- Reciprocal best hits (decided 2026-10-17): reference proteome = human RefSeq proteins (GCF_000001405.40, fetched
  with the genome; current, with gene symbols and the related families such as KTI12 and SECISBP2L); the RefSeq
  proteomes of OKI2018 and of the other annotated genomes serve for the reverse searches. Optional cross-check: the
  OIST Swiss-Prot BLAST/DIAMOND databases (2020_06, BioinfoUgrp modules), which are older and site-specific.
- Reciprocal best hits, first results (2026-10-18, human RefSeq proteome):
  - SecS (SEPSECS), eEFSec (EEFSEC) and SECp43 in F. borealis: reciprocal best hits; extra eEFSec copies only forward.
  - SECISBP2: in every tunicate, appendicularian or not, the single SBP2-family gene is closest to human SECISBP2L
    (margins -2 to -38 bits; H. roretzi +1), whereas vertebrates have both paralogs, cleanly separated. The tunicate
    gene is therefore a co-ortholog of the vertebrate SECISBP2 and SECISBP2L (duplication in vertebrates, to check
    in the literature, for example Donovan and Copeland 2009); the manuscript must name it accordingly.
  - SEPHS2: the selenocysteine-containing SEPHS of tunicates (and of zebrafish, weakly) are closer to human SEPHS1
    or only marginally to SEPHS2: best hits cannot separate SPS1 and SPS2 outside mammals (Mariotti et al. 2015).
    Use the residue at the catalytic site (selenocysteine or cysteine: SPS2 type) as the criterion.
  - PSTK: other tunicates are reciprocal best hits (55-108 bits); the O. dioica, O. albicans and B. stygius PSTK have
    no blastp hit in the human proteome at E <= 1e-3, O. vanhoeffeni 48 bits: too divergent for best hits against
    human; the absence test (tunicate references, KTI12 outgroup) and a phylogeny are the evidence.
  - Phylogenies needed: PSTK (with KTI12), SEPHS (SPS1/SPS2, with the Mariotti et al. 2015 sequences), SBP2 family
    (placement before the vertebrate duplication).
- Phylogenies, provisional (2026-10-18; sandbox, one CPU, LG+G4/R4 only; to redo on Deigo with the full model
  search and a site-heterogeneous model):
  - PSTK (154 sequences, rooted on KTI12): appendicularian PSTK form a clade (100) inside a tunicate group (74), not
    with KTI12; the PSTK/KTI12 split has only 54; the appendicularian clade is sister to insect-gut fungi (83), a
    long-branch artefact.
  - SEPHS (131, midpoint): vertebrate SEPHS1/SEPHS2 separate (100 each); tunicate selenocysteine copies (89) and
    glycine forms (97) separate; appendicularians have a second, non-selenocysteine SEPHS (Ala/Ser; Olon, Oalb, Ovan,
    OKI2018; 100) sister to the ascidian glycine form (95); the appendicularian SEPHS2 with selenocysteine (all nine)
    sit on a long branch outside the animal group. To check: the non-selenocysteine copy in OSKA2016 and OdB3; the
    Mariotti et al. 2015 SPS1/SPS2 sequence set.
  - SECISBP2 family (118, rooted on non-animal sequences; ClipKIT gappy 0.5, 304 columns; length filter 0.2 so that
    the domain-only tunicate predictions enter): oikopleurid proteins (99) nest inside the tunicate clade (87), sister
    to Thalia; tunicates group with protostome and echinoderm SECIS-binding proteins (39), and the vertebrate
    SECISBP2 and SECISBP2L clades (100 each) are not resolved against them (34-41): the tunicate gene is a single
    SECIS-binding protein that cannot be assigned to either vertebrate paralog. F. borealis is misplaced (13).
- Phylogenies, final (2026-10-19, Deigo, 8 CPUs; full ModelFinder: PSTK JTT+F+R6, SEPHS LG+R6, SECISBP2 family
  LG+I+R5; plus LG+C20+F+G for PSTK and SEPHS in results/trees/C20):
  - PSTK: appendicularian PSTK monophyletic (100) and inside PSTK, not KTI12 (PSTK/KTI12 split 94 with JTT+F+R6, 52
    with C20); their position inside PSTK is unstable, a long-branch effect (next to birds with JTT+F+R6, next to
    insect-gut fungi with C20, 81). Conclusion: PSTK, position within the family unresolved.
    CORRECTION (2026-10-19): the PSTK/KTI12 split has 100 with both models. The values 94 and 52 above, and the
    earlier provisional 54, came from a display error: Bio.Phylo moves supports by one branch along the rerooting
    path. Supports are now read from the unrooted bipartitions, and the figure scripts put them back on their
    bipartitions after rerooting (reroot_keeping_supports). Other values checked from the unrooted trees: PSTK
    appendicularian clade 100; SEPHS appendicularian selenocysteine clade 100 (C20 60), non-selenocysteine clade 100,
    with the tunicate glycine forms 94 (C20 95), tunicate selenocysteine copies 91 (C20 56); SECIS-binding proteins:
    oikopleurid clade 100, with Nematostella 70, other tunicates 90.
  - SEPHS: the non-selenocysteine appendicularian copies (Olon, Oalb, Ovan, OKI2018) group with the tunicate glycine
    forms (94; C20 95). The appendicularian selenocysteine SEPHS2 form their own clade (100, with F. borealis) at the
    base of the eukaryotic SEPHS (90 is the support of the whole eukaryotic group, not of a grouping with the glycine
    forms: corrected after drawing the collapsed tree); their position is unresolved (long branch). The
    selenocysteine copies of Styela, Halocynthia and Thalia group together (91; C20 56).
  - SECISBP2 family: the oikopleurid proteins form a clade (100, with the O. dioica profile sequence), sister to
    Nematostella (70), not with the other tunicates (90, among protostome and echinoderm proteins, 18): unlike the
    provisional tree, no support for grouping appendicularians with other tunicates. Conclusion: a member of the
    SECIS-binding protein family, orthology to the tunicate gene and to the vertebrate paralogs unresolved.
- Display trees (figures/phylogeny_<family>_collapsed.*, figures/phylogeny_collapsed.py): nodes with ultrafast
  bootstrap < 50 collapsed into polytomies; clades of a single taxonomic group (config/taxon_groups.tsv, assigned
  by hand from genus names) drawn as one tip; appendicularian and other tunicate sequences always shown.
- Long-branch analyses (stage2/machinery_lba.py, make machinery-lba; to run on Deigo): LG+C60+F+G with PMSF,
  slow-fast removal of the fastest sites (10 to 40 %), removal of non-target long branches (insect-gut fungi in
  PSTK), and an AU test of the grouping of appendicularians with the other tunicates. Background: O. dioica proteins
  evolve about twice as fast as those of Ciona, with homogeneous divergence across amino acids and a fast-evolving
  cysteine (Berna et al. 2012, doi:10.1016/j.ympev.2011.11.013; Berna and Alvarez-Valin 2015,
  doi:10.1016/j.margen.2015.07.007); the cysteine result also matters for the selenocysteine/cysteine homologs.
- Taxon sampling of the machinery trees (2026-10-19): added four tunicates (Botryllus schlosseri GCF_051294905.1,
  Phallusia mammillata GCF_965637545.1, the salps Pegea sp. GCA_034620365.1 and Salpa thompsoni GCA_047495815.1)
  to break the appendicularian long branch, and two cephalochordates (Branchiostoma floridae GCF_000003815.2,
  B. lanceolatum GCF_035083965.1) as chordate anchors. No echinoderm or hemichordate: the deuterostome clade itself
  may be an artefact of unequal branch lengths and compositional heterogeneity (Serra Silva et al. 2025, Current
  Biology, doi:10.1016/j.cub.2025.06.045), whereas chordates are well supported. SEPHS tree: the SPS1/SPS2-labelled
  sequences of Mariotti et al. 2015 (Supplemental Material S8; curation/mariotti2015_sps.faa). PSTK and
  SECIS-binding proteins: RefSeq proteins replace the partial Selenoprofiles predictions where a genome has them.
  Profile and reference sequences capped at three per taxonomic group (and SPS type).
- Machinery trees with the extended taxon sampling (2026-10-19, Deigo, commit 149cf6d; collapsed figures redrawn):
  - PSTK (86 sequences, 691 sites, LG+F+R6): PSTK/KTI12 split 100. Appendicularian clade 100 (OKI2018, OSKA2016,
    OdB3, O. albicans, O. vanhoeffeni; the O. longicauda, B. stygius and M. erythrocephalus fragments are below the
    length filter), sister to all other PSTK (91; C20 91). Other tunicates (ascidians and salps) form a clade (99)
    elsewhere in PSTK. Long-branch analyses: appendicularian clade 100 and PSTK/KTI12 100 in every analysis; never
    with the other tunicates; sister to the other PSTK with PMSF LG+C60 (90), 10-30 % of fast sites removed (86-92)
    and without the insect-gut fungi (89); with 40 % removed, next to the fungi (35). AU test: the tree with
    appendicularians and other tunicates together is not rejected (p-AU 0.275, delta logL 5.7). Conclusion: PSTK;
    position within the family unresolved; grouping with the other tunicates not rejected.
  - SEPHS (129 sequences, 976 sites, LG+R6): appendicularian selenocysteine SEPHS2 clade of 11 (all nine genomes
    except B. stygius, with F. borealis, M. erythrocephalus and the O. dioica SPS2 of Mariotti et al. 2015) 83 (C20
    89), position unresolved (32; C20 56). Non-selenocysteine copies (OKI2018, O. longicauda x2, O. albicans,
    O. vanhoeffeni) 100, inside the tunicate glycine forms (95; C20 100), next to Phallusia (78; C20 73). Ascidian
    selenocysteine copies 100; with the salps and Thalia 93 with C20 only. All tunicate SEPHS except the
    appendicularian selenocysteine clade: 78 (C20 92). Reading: the tunicate selenocysteine/glycine duplication
    precedes the divergence of appendicularians, so the appendicularian selenocysteine SEPHS2 is expected next to the
    ascidian and salp selenocysteine copies; it sits outside them on a long branch. Long-branch analyses pending.
  - SECIS-binding protein family (90 sequences, 391 sites, LG+R5): Oikopleuridae clade (with B. stygius and
    M. erythrocephalus and the O. dioica profile sequence CBY19450.1) 100, next to Nematostella (48), outside a group
    of the other tunicates, vertebrates and other bilaterians (71); F. borealis unresolved (34). Other tunicates 100,
    branching before the vertebrate SECISBP2/SECISBP2L duplication (vertebrate clade 86). Long-branch analyses
    pending.
- stage2/machinery_lba.py corrected (2026-10-19): (1) the p-AU value was read from a line of the rate table (0.1646
  for PSTK) instead of the AU table (0.275); (2) leaves were assigned from the groups table, which left out the
  reference and profile sequences of appendicularians and tunicates (O. dioica SPS2 and CBY19450.1, Botryllus SPS1
  and SPS2, Ciona XP_002129526) and put both appendicularian SEPHS paralogs in one set; the AU constraint for SEPHS
  and the SECIS-binding proteins therefore forced these sequences out of their clades. Now: leaves assigned by taxon
  (config/taxon_groups.tsv), the appendicularian clade is the largest appendicularian clade of the ML tree, the
  other tunicates are those of the orthologous form (SEPHS glycine forms left out), and a constrained tree that does
  not satisfy the current constraint is recomputed with its AU test. PSTK summary regenerated (trees unchanged).
- Long-branch analyses of SEPHS and the SECIS-binding proteins (2026-10-20, Deigo, commit 72da760;
  results/trees/lba/<family>.summary.tsv):
  - SEPHS: the appendicularian selenocysteine SEPHS2 clade holds with PMSF LG+C60 (79) and with 30 % of fast sites
    removed (84); with 10, 20 and 40 % removed, F. borealis leaves it (the other ten stay together, 100). Its position
    is unstable: at the base of the animal SEPHS in the ML, C20 and PMSF trees (32-56), and outside a group that
    includes kinetoplastid, stramenopile and amoebozoan sequences in the slow-fast trees, a position not possible
    for an animal gene. The tunicate SEPHS clade without it is supported in every analysis (ML 78, C20 92, PMSF 97),
    with the appendicularian non-selenocysteine copy inside the glycine forms (95-100; next to Phallusia 73-98)
    and the ascidian and salp selenocysteine copies together with C20 (93) and 40 % removal (83). AU test: the tree
    in which the appendicularian selenocysteine SEPHS2 groups with the ascidian and salp selenocysteine copies is
    rejected (p-AU 0.014, delta logL 52; LG+R6); in that constrained tree the glycine forms move next to the
    vertebrates. Reading: the tunicate duplication into selenocysteine and glycine forms precedes the divergence of
    appendicularians (their non-selenocysteine copy is a glycine-form ortholog); the appendicularian
    selenocysteine SEPHS2 is too divergent for its position to be inferred, and the AU rejection is made under a
    site-homogeneous model, exposed to the same long-branch effect. Its assignment rests on the selenocysteine at
    the catalytic site (SPS2 type; Mariotti et al. 2015) and on the presence of the separate glycine-form ortholog.
  - SECIS-binding protein family: the Oikopleuridae clade holds in every analysis (99-100). With PMSF LG+C60 it
    groups with the other tunicates (81), and with 40 % of fast sites removed too (71, F. borealis joining it,
    71); in the ML tree and with 10-30 % removed it sits next to Nematostella (below 50). AU test: grouping with the
    other tunicates not rejected (p-AU 0.503, delta logL 0.1). The vertebrate SECISBP2 and SECISBP2L form one clade
    that excludes the tunicates in every tree (86-100). Reading: the appendicularian protein is the ortholog of the
    single tunicate SECIS-binding protein gene, which branches before the vertebrate duplication (a pro-ortholog of
    both SECISBP2 and SECISBP2L); support comes from the site-heterogeneous model, not from the ML tree.
  - Summary of the three families: PSTK and SECIS-binding protein: appendicularian placement with the other
    tunicates not rejected (p-AU 0.275 and 0.503), supported for the SECIS-binding protein by PMSF C60 (81).
    SEPHS2: placement unresolved; the non-selenocysteine copy settles the timing of the tunicate duplication.
- SEPHS orthology AU tests (2026-10-20, Deigo, commit 559d9c2; results/trees/lba/SEPHS2.{ml,pmsf}.*): both positions
  of the appendicularian selenocysteine SEPHS2 clade are rejected under both models. Inside the clade of all tunicate
  SEPHS: p-AU 0.0067 (LG+R6, delta logL 43) and 0.019 (LG+C60+F+G PMSF, delta logL 34); next to the ascidian and salp
  selenocysteine copies, with the glycine forms and the second appendicularian copy as their sister: p-AU 0.012
  (delta logL 51) and 0.0082 (delta logL 39). Checks: no appendicularian selenocysteine SEPHS2 fails the IQ-TREE
  composition test (median p 74 %); its path from the common node of the tunicate and vertebrate SEPHS (2.2
  substitutions per site) is shorter than that of the second, non-selenocysteine appendicularian copy (3.0), which is
  placed inside the tunicate glycine forms in every analysis. Reading: under both models the data reject the
  appendicularian selenocysteine SEPHS2 as a member of the tunicate SEPHS clade, and neither composition nor the
  longest branch explains it. Two readings remain: a divergent ortholog of the tunicate selenocysteine copy whose
  signal these models do not recover, or a separate SEPHS2 lineage (with the ortholog of the tunicate selenocysteine
  copy lost in appendicularians). The manuscript states the rejection and leaves the origin open; the protein is a
  selenocysteine SEPHS2 by its catalytic residue in any case.
- Curation of the machinery tree sets (2026-10-21, after review of the PSTK alignment by eye: a truncated lancelet
  PSTK, rare insertions in mammals): stage2/machinery_phylogeny.py now (1) leaves out the human stage 1 predictions
  (the human RefSeq proteins represent human; the PSTK prediction carried an N-terminal extension and the SEPHS
  predictions included pseudogenes); (2) keeps identical sequences once (genome or RefSeq names preferred); (3) after a
  first alignment, removes sequences covering less than half of the core columns, except appendicularian, curated and
  tunicate sequences, then realigns; (4) trims with ClipKIT gappy 0.9 instead of smart-gap (smart-gap had kept 224 of
  691 PSTK columns and 426 of 976 SEPHS columns occupied by at most 6 sequences, for example the 75-residue insertion
  of mammalian KTI12, which is in the canonical human and mouse proteins). Removals are listed in
  results/trees/<family>.removed.tsv. Preview on the previous sets (local MAFFT 7.505): PSTK 83 sequences, 413
  columns (removed: human prediction, duplicate human KTI12, Branchiostoma floridae PSTK fragment of 106 residues);
  SEPHS 117, 448 (6 human predictions, 6 duplicates); SECIS-binding proteins 82, 427 (2 human predictions; fragments
  of Nematostella (78 residues, the ML sister of the oikopleurid clade), Helobdella, Trichoplax, Daphnia, Bodo and
  Sphaeroforma). Trees, long-branch analyses, figures, tables and manuscript numbers to recompute.
- Machinery trees after curation (2026-10-22, Deigo, commit 806e317; alignments PSTK 83 x 419, LG+F+I+R5; SEPHS
  117 x 454, LG+R6; SECIS-binding proteins 81 x 444, LG+I+R6). Removed: duplicate human KTI12 and the B. floridae
  PSTK fragment; 6 SEPHS duplicates; 7 partial SECIS-binding proteins (KOO29690.1, Daphnia, Bodo, Sphaeroforma,
  Trichoplax, Helobdella, Nematostella); human predictions (not listed in this run's removed.tsv; recorded from now on).
  - PSTK: PSTK/KTI12 100 and appendicularian clade 100 everywhere; other tunicates 97-100; the appendicularian clade
    now groups with the PSTK of the insect-gut fungi (56 ML, 74 C20, 84 PMSF, 60-77 site removal), and without the
    fungi it is sister to all other PSTK (94): long-branch attraction between two long branches. AU p 0.268.
  - SEPHS: oikopleurid selenocysteine SEPHS2 (10 sequences, 7 genomes) 100 everywhere; F. borealis joins under C20
    and PMSF only (83, 76). Placement: base of animal SEPHS under C20 and PMSF (39, 56), outside a group with
    kinetoplastid, stramenopile and amoebozoan proteins in ML and every site-removal tree. Tunicate SEPHS without it
    82-95; glycine forms with the second appendicularian copy 75-99; Phallusia sister 79 (ML), 95 (PMSF). AU: all
    rejected, more strongly than before (inside tunicate SEPHS 0.011 ML, 0.0045 PMSF; orthologous position 0.0047,
    0.0009; generic 0.0057). No SEPHS fails the composition test; paths from the tunicate-vertebrate common node:
    appendicularian selenocysteine copy 1.3, second copy 2.5 substitutions per site.
  - SECIS-binding proteins: Oikopleuridae clade 99-100; F. borealis joins it in all but the 40 % site-removal tree
    (47-77). Appendicularians with the other tunicates: 70 ML, 93 PMSF, 53-68 site removal (oikopleurids with
    tunicates 76 at 40 %); AU p 0.586. Vertebrate SECISBP2+SECISBP2L clade excluding tunicates and amphioxus 79-98.
    The earlier ML placement next to Nematostella was due to a 78-residue fragment.
  - Provisional Results, Table MACHINERY and Methods updated with these numbers; collapsed figures redrawn (title on
    two lines, legend below the axis label). The C20 run now uses --max-per-group 3 for its sequence collection too
    (its trees always used the copied trimmed alignments; only its unused faa and groups files had all sequences).
- Repository brought in line with the curated machinery analysis (2026-10-22): Table MACHINERY is now generated by
  stage2/machinery_table.py (make machinery-table, also run at the end of make machinery-lba) from the census, the
  reciprocal best hits, the long-branch table and the SEPHS trees, with the orthology column and the second SEPHS row
  computed from them; the manuscript table is pasted from results/tables/stage2_machinery_table.md. Manuscript: the
  machinery paragraph calls the appendicularian protein a SECIS-binding protein co-orthologous to SECISBP2 and
  SECISBP2L. READMEs: top-level steps (actual make targets), data sources (all 24 genomes, Mariotti et al. 2015
  references) and layout; stage2/README.md (curation, trees, long-branch analyses, table, run times); figures and
  results/tables indexes. Still open: results/tables/stage1_machinery_rbh.tsv (make machinery-rbh on Deigo); the
  second SEPHS copy in OSKA2016 and OdB3.
- Domain analysis of the machinery trees (2026-10-22, local run, 2 CPUs; stage2/machinery_domains.py, make
  machinery-domains): each trimmed alignment split at a domain boundary of the human protein (InterPro entries of
  Q8IV42, Q99611, Q96T21; curation/machinery_domains.tsv), the parts analysed on the same sequences. PSTK kinase
  1-263 (294 columns) / C-terminal (125), 82 sequences; SEPHS N-terminal 1-238 (234) / C-terminal (220), 116; SECIS-
  binding protein upstream of L7Ae 1-631 (253) / L7Ae 632-854 (191), 76.
  - In every family, each part rejects the tree of the other (AU p <= 0.0023), but no clade with >= 95 % support in
    one domain tree is incompatible with one so supported in the other, except arrangements within taxonomic groups
    and among prokaryotic and alveolate SEPHS.
  - PSTK: appendicularian clade 100 (kinase) and 97 (C-terminal), without the fungi; sister to all other PSTK (57) or
    next to dipteran PSTK (54); appendicularians with tunicates not rejected (0.44, 0.32). C-terminal data reject the
    full-length tree (0.0001), kinase data do not (0.19).
  - SECIS-binding protein: L7Ae tree as the full-length tree (oikopleurids with tunicates 81; vertebrate SECISBP2 +
    SECISBP2L 95); upstream region: oikopleurids next to vertebrate SECISBP2, with tunicates (67, 67). With tunicates
    not rejected (0.38, 0.49). F. borealis (315 aa) elsewhere in both.
  - SEPHS: N-terminal data reject no tested position (with tunicate Sec copies 0.45, inside all tunicate SEPHS 0.33,
    orthologous position 0.19); C-terminal data reject the orthologous position (0.0034; others 0.08). The N-terminal
    tree groups both appendicularian copies with dipteran and louse SPS2 (84; second copy with them 94), against the
    C-terminal tree (second copy in tunicate glycine forms 77-83; insect SPS2 with Rhodnius SPS2 100): the only
    conflict >= 90 % involving an appendicularian protein. The rejection of SEPHS2 orthology rests on the C-terminal
    half.
  - Provisional Results and Methods updated; display trees figures/phylogeny_<family>_collapsed_domain_<domain>.pdf;
    figures/phylogeny_collapsed.py got --formats.
- Schematic of the selenocysteine machinery (2026-10-22): figures/machinery_scheme.py (make machinery-scheme), drawn
  from scratch with SBGN-style glyphs (metabolite, RNA, protein, complex; conversion, catalysis, binding); PDF and
  editable SVG; draft legend in figures/README.md. Biochemistry checked against the abstracts of Xu2006 and
  SmallHoward2006. Note: SmallHoward2006 reports SECp43 complexes with SLA/LP (SEPSECS) and SPS1, not SPS2; the
  manuscript sentence "in complexes with SecS and SEPHS2" needs checking.
- Schematic placed in the introduction as Figure SCHEME (image from figures/machinery_scheme.svg, legend after it);
  SBGN reference added (LeNovere2009); the SVG now names sans-serif as fallback font.
- Deigo runs of 2026-10-23 (commit 872680d): make machinery-rbh and make machinery-domains-pmsf.
  - results/tables/stage1_machinery_rbh.tsv now in the repository: PSTK of the ten other tunicates reciprocal best
    hits of human PSTK (55-108 bits); SECIS-binding proteins closer to SECISBP2L in nine of ten tunicates and both
    amphioxus (H. roretzi SECISBP2 +1 bit); tunicate SEPHS mostly SEPHS1, Sec copies of B. schlosseri and S. clava
    SEPHS2 within 5 bits. Cited in the provisional Results. stage2_machinery_rbh.tsv: same values, prediction names
    without the _app suffix; Table MACHINERY unchanged.
  - Domains under LG+C60+F+G (PMSF): PSTK as under ML (with tunicates not rejected, 0.64 and 0.30; appendicularians
    next to the fungi 42 or to flies, Trichinella and Aphanomyces 70). SECIS-binding protein: L7Ae oikopleurids with
    tunicates 91, vertebrate SECISBP2+SECISBP2L 96; with tunicates not rejected (0.38, 0.51). SEPHS: the ML N-terminal
    grouping of the appendicularian copies with insect SPS2 disappears (second copy next to the tunicate glycine
    forms in both halves, 39-65); both halves reject the ortholog position (0.016, 0.011), neither rejects the
    position inside all tunicate SEPHS (0.22, 0.071). One strong conflict (95 vs 95): tunicate SEPHS with vertebrates
    (C-terminal) against vertebrates with amphioxus and protostomes (N-terminal).
  - Makefile fix: the PMSF domain figures were written as phylogeny_<family>_collapsed_domain_.pdf ($$d_pmsf read as
    one shell variable); now ${d}; misnamed files removed, figures redrawn.
- Orthology of the selenoproteins (2026-10-23, agreed: "apply it to the other proteins, then revisit"): the machinery
  approach generalised. machinery_rbh.py takes --families (sets machinery, selenoproteins; several human genes per
  family, margin to the best other gene); machinery_phylogeny.py adds the human RefSeq proteins of each selenoprotein
  family (outgroups GSR for TXNRD, SELENOM for SELENOF), takes the appendicularian predictions only from appendicularian
  genomes, skips families without appendicularian sequence or profile, and writes PDF only with --formats;
  machinery_lba.py skips families without tree; new stage2/selenoprotein_orthology.py (appendicularian clades, sister
  group, nearest vertebrate genes, best hits). Families: DIO, GPX, MsrA, SELENOF, SELENOK, SELENOL, SELENOO, SELENOS,
  SELENOT, SELENOU, SELENOW, TXNRD (those with a selenocysteine prediction in an appendicularian genome, or in the
  baseline). Tested on synthetic data (human sequences from the InterPro API, mutated copies as predictions). To run
  on Deigo: make selenoprotein-rbh, selenoprotein-phylogeny, then selenoprotein-pmsf; selenoprotein-orthology after.
  Known gap: no appendicularian SELENOK in stage 1.
- Selenoprotein orthology, first pass (2026-10-24; Deigo commit 07488b7; summary by make selenoprotein-orthology,
  results/tables/stage2_selenoprotein_orthology.{tsv,md}, now with labels and curated O. dioica gene names; trees read
  in figures/selenoproteins/*_collapsed_pmsf.pdf):
  - TXNRD: appendicularian clade (9 genomes, 100) with the tunicate TXNRD1/3-type proteins (89-92), co-ortholog of
    TXNRD1 and TXNRD3 (90-95); no appendicularian protein in the TXNRD2 clade, where ascidians and salps have one.
  - GPX: the four O. dioica Sec GPX in two appendicularian clades (chr1a+chr1b; PAR+chr2), both in the group of
    vertebrate GPX4 and tunicate GPX4-like proteins (70-75 PMSF); best hits GPX4 for most Sec predictions. No
    appendicularian protein in the GPX1/2 or GPX3/5/6 groups (tunicates have GPX3-like). Cys GPX7/8-like proteins in
    O. longicauda and O. vanhoeffeni (97). F. borealis Sec GPX clade outside the GPX4 group.
  - SELENOW: three lineages. chr1a, pair1 and pair2 with Halocynthia SELENOW in the SELENOW clade (expanded, tandem
    pair). chr2 (all nine genomes, 95-99) with tunicate, amphioxus, insect and fish SelW-like proteins whose best human
    hit is MIEN1. PAR (OKI, OdB3, Olon, Mery; 97-100) sister to mammalian SELENOV (79-82), best hit SELENOV (ambiguous);
    SELENOV arose in placental mammals (Mariotti et al. 2012), so the grouping is suspect: to test.
  - DIO: DIO_PAR and DIO_chr1 among the tunicate DIOs (with protist deiodinase-like and sea urchin sequences nested,
    weak supports); invertebrate chordate DIOs nearer DIO1 than DIO2/DIO3; best hits mixed DIO1/2/3. Unresolved.
  - SELENOF, SELENOO, SELENOT: reciprocal best hits of the human gene; trees place the appendicularian clade outside
    the animal or tunicate clade (SELENOF between plants and animals; SELENOO next to Phytophthora SelO; SELENOT next
    to nematodes), the pattern of long branches: to test (AU) before calling them orthologs.
  - SELENOK: the curated O. dioica SELENOK is not in stage 1; the appendicularian predictions in the tree are
    glycine/tyrosine/proline homologues, not it. SELENOS: only the OSKA2016 protein, no human hit.
  - Families absent from the Oikopleuridae but in other tunicates: SELENOH, SELENON, MSRB, FrnE (stage 1); SELENOL and a
    Sec SELENOU only in F. borealis; appendicularian MsrA and SELENOU (PRXL2-like) are cysteine forms.
  - F. borealis, and to a lesser extent B. stygius and M. erythrocephalus, carry many cysteine homologues that group
    with bacteria or fungi (GPX, DIO-like thioredoxins, SelO, SelW-like, GSR-like): likely contaminants; exclude from
    the census or flag.
  - Profile sequences named by Ensembl identifiers are unclassified and kept in full (SELENOO): add these codes to
    config/taxon_groups.tsv before final trees.
  Next (to agree): AU tests for SELENOW_PAR (SELENOW vs SELENOV), SELENOF, SELENOO, SELENOT; MIEN1 added to the SELENOW
  references; a curated SELENOK; contamination filter.
- Step 1 done (2026-10-24, local; stage2/secis_validated_hits.py, make secis-hits; results/tables/stage2_secis_hits*.tsv):
  111 selenocysteine predictions in the six Oikopleura assemblies (OKI 12, OSKA 14, OdB3 21, Olon 28, Oalb 20,
  Ovan 16), none left out by the rule; 14 within 500 nt of a contig end (11 in O. longicauda). RF00031: 41 with a
  hit, 15 at 20 bits or more, 4 at 25; best random window per genome 13.2-22.6; no cysteine homolog at 20 bits.
  Current model: OSKA2016 (held out) 12 of 14 at 25 bits, best random window 18.1, no cysteine homolog at 20.
  Families with the least SECIS evidence: SELENOO, SELENOS, SELENOT (none at 20 bits with RF00031), SELENOF (1 of 7).
  Provisional paragraph and Methods added to the manuscript.
- Step 2 started (2026-10-24): curation/gold_standard.py (make gold-standard), tested on synthetic data (two
  strands, a GA-AG intron, an unstranded transcript, a transcript with an insertion, a prediction-only gene, an
  annotation starting after the selenocysteine codon). minimap2 -x splice moved a GA-AG junction by two bases and
  created a false insertion; -un (no GT-AG preference) aligns it exactly, so the script uses -x splice:hq -un. To run
  on Deigo; then review each gene against the pilot (SELENOT GA-AG intron and selenocysteine exon, SEPHS2 first exon,
  OKI2018 SELENOS allele, OdB3 SELENOO C-terminus, SELENOK, the SELENOW genes without RefSeq model).
- Step 2, first Deigo run reviewed (2026-10-24, commit 98a5e95; 54 curated genes, 52 models). Four faults, checked on
  the ENA sequence of OKI2018 (OU015566-OU015570, = RefSeq NC_137291-NC_137295) and of the OdB3 scaffolds:
  (1) with -un, minimap2 places a junction anywhere among the placements that give the same spliced sequence: 53 of
  the 110 introns of the OKI2018 and OdB3 models become GT-AG, GC-AG or AT-AC by such a shift (mRNA unchanged); of
  the 18 others, GA-AG is one of the equivalent placements of 12; (2) the protein was read in the transcript and the
  GFF3 coding sequence on the genome, which disagree when they differ by indels (OKI DIO_chr1: 13-nt insertion in
  the transcript, 110 aa; OKI SELENOO: 1 nt, 27 aa); (3) a complete frame ranked before a longer one; (4) the 3' UTR
  and after-stop FASTA files were excluded by the *.fa rule of .gitignore. Other observations: TXNRD is 600 aa in
  OKI2018 and OSKA2016, with an N-terminal segment carrying a glutaredoxin-like CPFC motif (domain to be verified);
  OKI SELENOW_pair1 and SELENOW_pair2 were both modelled from one transcript (TRINITY_DN2389_c0_g1_i2), which holds
  the two frames (relevant to step 10); OSKA2016 SEPHS2 and OdB3 SELENOO had neither transcript nor prediction.
  Revision of curation/gold_standard.py: junction placement by the script (fewest mismatches in +/- 8 nt, then
  canonical sites), normalisation of every exon chain to canonical sites among equivalent placements; candidates
  from transcripts, annotated mRNAs and the stage 1 prediction, all read on the genome sequence (internal check that
  the GFF3 coding sequence gives the protein); ranking by curated context (either side of the U, since the curated
  context was read on contiguous genome and runs into introns for GPX_chr1a, SEPHS2 and SELENOT), length within 90%
  of the longest (from the first methionine when the 5' end is open), completeness, transcript support, annotated
  chain; an open 5' end starts at the first ATG, or at the first codon for a fragment without ATG (3' ESTs);
  assembly or transcript protein decided by identity to the same gene in the other assemblies (margin 0.02),
  transcript-based models marked with the GFF3 attribute exception; annotation_fix column; transl_except comma
  percent-encoded; the FASTA outputs committed. Tested on the ENA sequence with simulated transcripts (a 13-nt
  insertion in transcripts, a 13-nt deletion in an assembly copy, a 1-nt insertion, 5'-truncated ESTs) and
  RefSeq-style annotations of each status. A PubMed search found no report of GA-AG introns in O. dioica (two hits, both on
  F. borealis: Henriet et al. 2019, doi:10.1016/j.cub.2019.07.092; Soo et al. 2026, doi:10.1093/nar/gkag659); this
  was wrong, see the next entry. To rerun on Deigo (the
  alignments are cached in work/stage2/gold_standard/).
- Intron classes of O. dioica (2026-10-24, references provided by C. Plessy): Denoeud et al. 2010, supplementary
  section 10.1 (doi:10.1126/science.1194167): in the reference annotation (73,254 introns, 72.6% validated by
  cDNAs), GT-AG 87.6%, GA-AG 9%, GC-AG 2.7%, GG-AG 0.7%; no AT-AC or U12 intron and no minor spliceosome; GA-AG
  introns mostly 55-100 nt (GT-AG 35-55 nt), with a C...AAAG acceptor consensus; 24.2% of UTR introns and 11.2% of
  CDS introns are not GT-AG; the annotation allowed all G*-AG donors. Frey and Pucker 2020 (Cells 9:458,
  doi:10.3390/cells9020458): GA-AG more frequent than GC-AG in O. dioica and the copepod Eurytemora affinis. The
  PubMed record of Frey and Pucker 2020 lost the italic species names of its abstract ("observed inand"), which
  explains why the search missed it. curation/gold_standard.py now uses these classes: splice-site cost
  log10(f(GT-AG) / f(site)) / 2 (GA-AG 0.49, GC-AG 0.76, GG-AG 1.05, other sites 1.47, in mismatches; AT-AC no
  longer preferred), equivalent placements moved to the most frequent class, and the table reports atypical
  introns (GA-AG, GC-AG, GG-AG) apart from other splice sites, which point to an alignment, assembly or transcript
  problem. On the test set (ENA sequence, chains of the first run): the SELENOT intron placed as GA-AG in OKI2018
  and OdB3 (as in the pilot), and the DIO_PAR intron too; other splice sites remain in OKI GPX_chr1a (1),
  SELENOO (4) and SEPHS2 (1), to review on the rerun. Relevant to step 7 (stage 2 with non-canonical introns):
  G*-AG donors for O. dioica.
- Short introns (2026-10-24, noted by C. Plessy): O. dioica introns peak at 40-44 nt, with a few percent below 35 nt
  (Denoeud et al. 2010, figure S22; the cDNA mapping used exonerate --minintron 25, so shorter introns were not
  measurable). minimap2 2.26 -x splice:hq -un aligns introns of 20 nt or more as introns and those of 16 and 18 nt as
  deletions (synthetic test, GT-AG and GA-AG alike); the container has minimap2 2.30 (to check on Deigo; a release
  binary could not be run here). curation/gold_standard.py now reports deletions of 10 nt or more with an equivalent
  placement on a GT-AG, GA-AG, GC-AG or GG-AG site as possible introns (notes column); the junction window (+/- 8 nt)
  and the end tolerance (10 nt) stay below the shortest intron sizes. The pipeline notes of the pilot
  (docs/oikopleura_sec_annotation_agent.md in the pipeline repository) already ask to allow GA-AG and short introns
  and record that lowering the minimum intron length of exonerate did not recover the missed Sec exons. For step 7:
  check the minimum intron lengths of exonerate and GeneWise as called by Selenoprofiles 4.5.8 (source not
  reachable here).
- Step 2, second Deigo run reviewed (2026-10-24, commit 150b560): 53 of 54 genes modelled (OdB3 SELENOO: no
  transcript, annotation or prediction holds the codon, as in the pilot). All 32 assembly-based models of OKI2018
  and OdB3 translate to their protein on the ENA sequence. Agreements with the pilot: OKI2018 SELENOS, the assembly
  lacks 13 nt (transcript protein, 196 aa); SELENOT and SEPHS2 need only a new CDS start and the Sec declaration on
  the RefSeq exon chains; DIO_PAR, GPX_PAR and SELENOF are the RefSeq models with Sec declared and need no change;
  the SELENOT intron is GA-AG in OKI2018, OSKA2016 and OdB3. Transcript proteins chosen where the assembly lacks or
  adds bases: OSKA2016 SELENOW_PAR (2 nt; identity 0.91 against 0.72) and OdB3 GPX_PAR (1 nt; 0.99 against 0.91).
  OKI2018 SELENOW_PAR and SELENOW_pair1/pair2 lie in the 3' part of the RefSeq mRNAs of other genes
  (XP_081636423.1; XP_081628722.1, 5' UTRs of 1,976, 2,818 and 2,336 nt), as does OdB3 SELENOW_PAR in CBY23864.1:
  transcripts that hold another coding sequence upstream (step 10, operons). Faults: (1) many O. dioica mRNAs have
  short 5' UTRs without an in-frame stop (trans-splicing), so ranking on a complete frame tied transcript models
  with 3'-truncated annotated ones and the annotation preference then chose the truncated model (OdB3 SELENOF,
  141 aa without stop instead of 144 aa); (2) the annotation preference came before length (OKI2018 TXNRD: RefSeq
  545 aa chosen over the 600-aa transcript model); (3) identity as a fraction favoured a partial prediction with a
  GT-CG junction over the EST protein (OdB3 DIO_PAR); (4) chosen predictions had no UTRs (OKI2018 SELENOW_chr1a);
  (5) the OdB3 annotation splits exons at the TGA. Revision: ranking by context, 90% length band, stop codon, ATG
  start, transcript support, chain with UTRs, length, then annotation (upstream stop reported, not ranked); identity
  counted over the part of each reference protein that both proteins cover; UTR exons grafted from compatible
  transcripts onto chosen CDS-only models; touching blocks merged; candidates table per assembly
  (curation/gold_standard/<genome>.candidate_models.tsv) for review. To review after the rerun: OKI2018 DIO_chr1 (model
  from a transcript with 56 mismatches and a TT-AA junction; RefSeq starts after the Sec), OKI2018 SELENOO (RefSeq
  587 aa; transcripts incompatible), and the OSKA2016 genes, all reported as new exon chains against the OSKA2016
  annotation.
- Step 2, third Deigo run reviewed (2026-10-24, commit 330d19f), with the candidate tables and the ENA sequence:
  31 assembly-based models of OKI2018 and OdB3 translate to their protein. Checked against the orthologues of the
  other assemblies: OKI2018 has assembly errors in three genes, where the transcript protein is used: SELENOS
  (13 nt), DIO_chr1 (frameshift in the 3' part: the assembly frame ends ...WCENWRKF, the transcript protein
  ...VEYVGGPGPFGYDLDEMVRKLEEILKDN as in OSKA2016 and OdB3; RefSeq bridges it with a 2-nt gap and gives another
  C-terminus; the 7-nt terminal exon of the transcript alignment is an artefact of the broken region) and SELENOO
  (5 nt; the transcript protein has 536 residues identical to OSKA2016 against 519 for the RefSeq chain, which
  lacks a 15-aa segment, DTTVIRDKFYDGNAK, present in OSKA2016). OKI2018 TXNRD: transcript model of 600 aa chosen
  (RefSeq 545 aa). OKI2018 DIO_PAR: RefSeq places an intron 11 nt from the transcript-supported GA-AG intron
  (DPYAAKKPIR against DPYGQTKPIR; orthologues D[ST]DAQLKPIR, not decisive), and the transcripts lack the RefSeq
  first exon. OSKA2016: the annotation has CDS rows only and its models splice around the TGA (for example GPX_PAR,
  542884 joined to 543323), so new exon chains are needed for all genes. OdB3: DIO_PAR, eight ESTs give the same
  complete 243-aa protein and differ from the assembly by one base, whose frame breaks (only a 5' fragment
  prediction with a GT-CG junction was a candidate); DIO_chr1, the EST lacks the first 12 nt and the prediction
  reads through a 37-nt intron; SEPHS2, ESTs end with frame errors and the assembly protein is kept. Changes:
  transcript chains extended with annotated or predicted exons beyond their ends (tested on the ENA sequence: OKI
  DIO_PAR 243 aa with the GA-AG intron; OdB3 DIO_chr1 220 aa, MWAQIF...SPIQCY as in OSKA2016); transcript protein
  for comparison = the one shared by the most transcripts; a complete transcript protein replaces a fragmentary
  assembly protein more than 10% shorter unless its identity over the shared part is lower by more than 0.02;
  annotation_fix reports intron moves.
- Step 2, fourth Deigo run reviewed (2026-10-24, commit 0450c82): OKI2018 DIO_PAR now has the transcript-supported
  GA-AG intron and the RefSeq first exon (243 aa; annotation_fix: move 1 intron of XP_081636384.1); OdB3 DIO_chr1
  220 aa; OdB3 DIO_PAR uses the complete EST protein (244 aa); 30 assembly-based models translate on the ENA
  sequence. New faults: (1) the transcript protein "shared by the most transcripts" was a truncated one for OdB3
  GPX_PAR (167 aa instead of 183: several ESTs start at the same internal ATG); (2) extensions beyond the 5' end of
  a complete transcript model added N-terminal residues without support (OSKA2016 SELENOW_pair1 103 aa against 93
  in OKI2018 and OdB3; OKI2018 TXNRD 611 against 600), which matters for trans-spliced mRNAs whose 5' end is the
  true end; (3) a predicted exon brought a TC-AT junction without transcript support (OdB3 TXNRD). Changes: such an
  N-terminal extension is kept only when its residues align to the same gene in another assembly (half of them at
  40% identity; tested on the proteins of the fourth run: kept for OdB3 DIO_chr1, OKI2018 and OSKA2016 DIO_PAR,
  rejected for OSKA2016 SELENOW_pair1 and the TXNRD extensions, whose residues differ between assemblies), and kept
  but marked as not tested when no other assembly has the gene; a column start_in_transcripts; models with a splice
  site outside the intron classes that no compatible transcript shares rank below the others; the reference
  transcript protein excludes transcripts with more than 20 indels and prefers, within 95% of the longest, the one
  that contains the most other transcript proteins.
- Step 2, fifth Deigo run reviewed (2026-10-24, commit 16a1fcd): as expected for OdB3 GPX_PAR (183 aa), OSKA2016
  SELENOW_pair1 (93 aa), OKI2018 TXNRD (600 aa), OSKA2016 SELENOS (9-aa extension rejected) and OdB3 TXNRD (no
  unsupported junction; EST completed by the prediction, 609 aa); 30 assembly-based models translate on the ENA
  sequence. Remaining fault: OdB3 SEPHS2 had no stop codon (323 aa; OKI2018 393, OSKA2016 391): the OdB3
  annotation splits the gene at the TGA (CBY08599.1, N-terminal part spliced around the TGA by a 46-nt intron;
  CBY08598.1, C-terminal part with the stop codon), and extensions used only models holding the codon. Change:
  transcript chains are extended with any overlapping annotated mRNA or prediction, when the extension changes the
  coding sequence (tested on the ENA sequence: OdB3 SEPHS2 392 aa, C-terminus as in OSKA2016); annotation_fix
  reports merges of compatible split models.
- Step 2 done (2026-10-24, sixth Deigo run, commit 151eea3): OdB3 SEPHS2 complete (392 aa), merges of split OdB3
  models reported; other proteins unchanged; 30 assembly-based models translate on the ENA sequence. Summary: 53
  models for 54 curated codons (OdB3 SELENOO: none), 52 complete (OdB3 GPX_chr1a: 5' fragment), 50 with a 3' UTR
  (median 300 nt); transcript proteins where the assembly lacks or adds bases: OKI2018 SELENOS, SELENOO, DIO_chr1;
  OSKA2016 SELENOW_PAR; OdB3 DIO_PAR, GPX_PAR; 195 introns, 19 GA-AG, 6 GC-AG, 1 GG-AG (13%), 10 other sites in 6
  models. Provisional Results and Methods sections added to the manuscript; Biopython cited (Cock2009). Next: step
  3, SECIS model from the gold-standard sequences after the stop codons (curation/gold_standard/*.after_stop.fa) and
  the stage 1 predictions of O. longicauda, O. albicans and O. vanhoeffeni.
- Step 3 started (2026-10-24): stage2/secis_model_training.py (make secis-v2). Windows of 2,000 nt after the stop
  codons (gold standard for the three O. dioica assemblies; stage 1 predictions kept by step 1, extended to the stop
  codon, for O. longicauda, O. albicans and O. vanhoeffeni), shortened at the next selenoprotein coding sequence
  (OKI2018 SELENOW_pair2: 197 nt before SELENOW_pair1; minimum window 100 nt); negatives: random windows of the six
  assemblies and of Cint, Csav, Scla, Tdem (training and evaluation seeds apart), dinucleotide-shuffled positives,
  windows after cysteine homologs. Bootstrap from RF00031, one SECIS per locus, 90% identity filter (Biopython global
  alignment with free end gaps; no vsearch or CD-HIT in the image); merged model per fold from the RF00031 seed
  sequences and the fold's seeds aligned together to RF00031 (no --mapali: it requires the Rfam alignment to give
  the consensus columns of the Rfam model exactly). Folds per family and per species; RF00031 and the current model
  scored on the same windows and thresholds (current model optimistic: trained on most held-out sets). Whole-genome
  scans are left to step 4 (genome-based UGA/SECIS screen). Tested locally (Infernal 1.1.5 from Debian) on the ENA
  sequence of OKI2018 and OdB3, a pseudo-genome with stage 1-style predictions (four of six truncated predictions
  extended to the gold-standard stop; the two others are OdB3 genes whose assembly has an extra base) and a synthetic
  non-target genome, with a stand-in for RF00031 (rfam.org is not reachable here).
- Reannotation figures (2026-10-24): figures/reannotation_inputs.py gathers, for one gene, the public models that
  overlap the gold-standard models (candidate tables of step 2; RefSeq CDS located by finding the RefSeq protein in
  the translated exon chain; GenBank models from the ENA EMBL records), reference proteins and their Pfam domains
  (InterPro API), with the source of each value; figures/reannotation_figure.py draws them (SEPHS2, DIO_PAR; human
  DIO2 and Halocynthia roretzi Q6U6H1, a complete tunicate deiodinase with U133, as references for DIO_PAR; most
  tunicate deiodinases in UniProt are fragments of 80-130 aa that stop at the TGA). While drawing DIO_PAR, the
  annotation status of curation/gold_standard.py was found to report a codon inside an intron of the annotated
  coding sequence as "coding sequence starts after the codon" (OSKA2016 DIO_PAR and SEPHS2); it now reports
  "codon in an intron of the coding sequence (skipped)", class of Tico et al. 2026; the table changes at the next run
  of make gold-standard.
- Reannotation figures, layout (2026-10-24, comments of Charles): panel A aligned on the selenocysteine TGA in every
  assembly (one TGA line); splice sites marked by the colour of the intron chevron (GA-AG, GC-AG, GG-AG; other
  sites black), so that drawn lengths stay unchanged; the sites of public and gold-standard introns are read on the
  genome by figures/reannotation_inputs.py, or taken from the candidate table; public models must overlap the
  coding sequence of the gold-standard model (not only its UTR); Pfam short names. SEPHS2 and DIO_PAR to be redrawn.
- Third example, TXNRD (2026-10-24): RefSeq XP_081643508.1 (OKI2018) has other 5' exons, without the glutaredoxin
  domain (its first 49 residues do not align to the gold-standard protein, identical from residue 104 on), and its
  coding sequence ends at the TGA; the ANISEED model (OSKA2016) has an extra exon inside the GA-AG intron and ends in
  an intron before the TGA, and is joined to 6 or more exons upstream; GenBank CBY19735.1 (OdB3) has the GA-AG
  intron but one intron away from the ESTs, and ends at the TGA. Gold-standard defect found: the OdB3 model took its
  5' part (beyond the ESTs, which cover the 3' part) from the Selenoprofiles prediction rather than from GenBank
  (609 against 600 aa, ranked by length): an N-terminal extension of 12 aa not found in OKI2018 or OSKA2016 and a
  GT-AG intron in place of the GA-AG intron of GenBank and of the OSKA2016 transcript (VNAKATE, residues 105-111
  of the OSKA2016 protein, becomes GLIQ). The N-terminal test did not apply because the EST's own model is a fragment. New rule in
  curation/gold_standard.py: among the candidates equal to the chosen one on every criterion before protein length,
  the one whose protein aligns best to the same gene in the other assemblies (mean global alignment score) replaces
  it when its score is higher by 20 (about 10 bits). From the candidate tables (OKI2018 and OdB3, local genomes),
  the rule changes only OdB3 TXNRD (mean score 3031 against 2986; OKI2018 TXNRD: 2945 against 2940, unchanged).
  Also: annotation_status now separates "coding sequence ends at the codon (TGA read as stop)" from "codon inside
  the coding sequence, not declared" (a CDS that ends with the TGA as stop codon was counted as inside: GenBank
  CBY19735.1, RefSeq XP_081643508.1); the candidate table gives the splice sites of the annotated models that are
  not candidates and the annotated CDS rows (new column annotated_cds). Tested on a reconstruction of the OdB3 and
  OKI2018 TXNRD loci (transcripts spliced from the genome, annotation and prediction chains of the tables). The
  results paragraph on the public annotations needs its counts again after the rerun. The TXNRD figure will be
  drawn from the rerun (a preview was made with the expected OdB3 model).
- Multi-gene page (2026-10-24, request of Charles): figures/reannotation_page.py draws an A4 portrait page, three
  genes per row (columns), each with its gene models and proteins; assembly names and gold-standard protein labels
  once per row on the left; gene structures half as thick (CDS 1.4 mm, UTR 0.55 mm); public models named, with the
  curator's note, in a caption wrapped above them (letters a, b, c when several share a row); reference proteins
  labelled with their accession on the same line above their bar; less space between proteins; key at the bottom of
  the page. One x scale per gene. A test with the three genes repeated shows that three rows (nine genes) fit on the
  page. Previewed with SEPHS2, DIO_PAR and TXNRD (inputs regenerated with the new inputs script; expected OdB3
  TXNRD model).
