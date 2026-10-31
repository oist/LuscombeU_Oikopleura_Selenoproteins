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
- Second row of the page (2026-10-24): SELENOT (RefSeq and GenBank: same exon chain, coding sequence from an
  internal ATG 62 codons after the TGA; ANISEED: TGA in an intron of the coding sequence, extra 5' exon, ends in the
  GA-AG intron), DIO_chr1 (RefSeq "low quality protein", corrected model starting 4 codons after the TGA; ANISEED:
  TGA exon skipped, joined to 5 exons upstream; GenBank: two models split at the TGA, CBY06720.1 ending with the TGA
  as stop and CBY06719.1 starting 4 codons after it; OKI2018 gold model mapped from the transcript), SELENOW_PAR (no
  annotated coding sequence: in the 3' UTR of the RefSeq and GenBank mRNAs of another protein, 592 and 589 aa; no
  OSKA2016 model; OSKA2016 gold model mapped from the transcript). References: human SELENOT, DIO3, SELENOW
  (UniProt); Halocynthia Q6U6H1; Ciona SELENOT and SELENOW from the stage 1 predictions on the Ciona intestinalis
  genome (results/trees/selenoproteins; UniProt has only fragments or cysteine homologs). Commands with the
  curator's notes: figures/reannotation_inputs.sh.
- Gold-standard defect found with SELENOW_PAR (2026-10-24): in OKI2018 and OdB3 the model was the annotated mRNA of
  the upstream protein (XP_081636423.1, CBY23864.1), which holds SELENOW_PAR in its 3' part, because annotated chains
  were preferred among candidates with the same protein; the transcripts start 11 nt before the start codon. Same
  for SELENOW_pair1 and pair2 (OKI2018, mRNA of XP_081628722.1) and SELENOW_pair2 (OdB3, CBY21828.1). New rule:
  an annotated mRNA whose coding sequence does not overlap the candidate's is ranked after the other candidates of
  the same length and is not preferred as annotated chain; the notes report that it holds the gene. Tested on a
  reconstruction of the OKI2018 locus (transcript chosen, note written). The operon question (step 10) is unchanged:
  the transcripts of the pair genes still hold both open reading frames.
- Page, other changes: public exons cut at the window instead of dropped (the RefSeq and GenBank mRNAs of
  SELENOW_PAR were not drawn); stretches longer than 1,500 nt without any drawn exon drawn 200 nt long with break
  marks (the 5.2 kb intron of OKI2018 DIO_chr1, 176 nt in OSKA2016); proteins placed by local alignment, the
  OSKA2016 and OdB3 gold-standard proteins through the OKI2018 one, unaligned residues drawn as a line (global
  alignments with free end gaps placed the weakly conserved N-termini of the three DIO_chr1 proteins differently);
  "other frame" only for public coding sequences that overlap the gold-standard one; RefSeq proteins corrected for
  the assembly located by their first and last 15 residues.
- Gold-standard rerun (2026-10-24, pushed by Charles, commit f9ae243) with the three new rules: OdB3 TXNRD from
  the ESTs and the GenBank model (600 aa, GA-AG intron), OKI2018 and OdB3 SELENOW_PAR, OKI2018 SELENOW_pair1 and
  pair2 and OdB3 SELENOW_pair2 from the transcripts (notes: the annotated mRNA of another protein holds the gene),
  as expected. One more change from the similarity rule: OSKA2016 SELENOW_pair2 is now comp18928_c0_seq3 (94 aa, one
  exon) instead of comp18928_c0_seq2 (101 aa, an intron with a GC-CG site; mean score 468 against 332). Annotation
  status: 12 coding sequences end at the TGA read as stop (see results/tables/gold_standard.tsv),
  12 skip the TGA in an intron, 9 start after it, 4 hold it without declaration, 3 declare it, 3 end before it, 10
  genes have no annotated coding sequence. The results paragraph on the public annotations must be rewritten with
  these counts.
- Page layout, second version (2026-10-24, comments of Charles): the notes moved from captions above the models to
  one paragraph per gene between the gene models and the proteins (public models of each assembly with letters as
  on the models, notes on the gold-standard models, accessions of the reference proteins); SECIS drawn as a flat
  green box on the model line at its position, bit score after the end of the model; public and gold-standard tracks
  closer (3.4 mm); one protein row per species, labelled on the left (human; C. intestinalis), the reference
  accessions in the paragraphs; all tunicate references are now the C. intestinalis stage 1 predictions
  (results/trees/), so that every row of genes has the same species rows (Halocynthia Q6U6H1 and UniProt H2Y2A4
  no longer used); assembly names horizontal on the left. Third row: GPX_PAR (RefSeq declares the selenocysteine,
  same model as the gold standard; a RefSeq non-coding RNA ends in the TGA exon; ANISEED skips the TGA exon;
  GenBank starts 14 nt before the TGA in another frame for 41 codons, the OdB3 assembly having one extra base),
  SELENOK (no RefSeq or ANISEED model; GenBank reads the TGA as stop), SELENOS (RefSeq loses the frame where the
  OKI2018 assembly lacks 13 bases; ANISEED reads the TGA as stop; GenBank splits the gene in two). Nine genes fill
  the page; figures/reannotation_page.{pdf,svg,png} committed (PNG at 200 dpi, under 1 MB).
- Page 2 of the reannotation figure (2026-10-24): GPX_chr1a, GPX_chr1b, GPX_chr2, SELENOF, SELENOO, SELENOW_chr1a,
  SELENOW_chr2, SELENOW_pair1, SELENOW_pair2; public and gold-standard tracks of an assembly closer (2.4 mm against
  4.2 mm between assemblies). Human references: GPX4 for the three GPX genes, SELENOF, SELENOO, SELENOW; for
  SELENOW_chr2 the best human hit is MIEN1 (Q9BRT3), a cysteine homolog: the page marks the residue of the first
  reference aligned to the gold-standard U when the reference has no selenocysteine. Gene with no gold-standard model
  in one assembly (SELENOO, OdB3: exon of the TGA not assembled) drawn with a note.
- Gold-standard defects found while preparing page 2, and rules added to curation/gold_standard.py (rerun needed):
  1. OdB3 GPX_chr1a was the Selenoprofiles prediction (146 aa, 5' fragment without ATG, reading 22 nt into an
     intron), while GenBank splits the gene at the TGA into CBY23346.1 (ending with the TGA as stop) and CBY23347.1
     (starting 13 codons after it), with no EST to join them. New candidates: an annotated model holding the codon
     and ending at most 30 nt after it, joined to the next annotated model downstream (at most 300 nt further), the
     stretch between them read as exon. For GPX_chr1a: 183 aa, ATG, stop, GT-AG introns, mean alignment score to
     OKI2018 and OSKA2016 913 and 976 against 683 and 741 for the prediction.
  2. OdB3 SELENOW_chr1a was the prediction (100 aa, reading through an intron) instead of the 15 spliced ESTs (86 aa,
     GC-AG intron, like OKI2018 and OSKA2016): the prediction was longest and the ESTs fell under 90% of its length.
     New criterion, before the length criterion: not contradicted by the transcripts (more supporting transcripts
     incompatible with the chain than compatible with it); contradicted chains are also left out of the longest
     length. From the candidate tables, only OdB3 SELENOW_chr1a changes.
  3. OSKA2016 SELENOW_chr2: the only transcript reads through into the next ANISEED gene (S40.g16479), its 3' UTR
     holding five introns and that gene's exons. New rule: 3' UTR exons after the exon of the stop codon are left out
     from the first one that overlaps the coding sequence of another annotated model. Only this model has such
     exons.
  4. The annotated models listed near a gene were those within 5 kb of the codon; now also those overlapping any
     candidate (OSKA2016 SELENOO, 8.3 kb, has two more ANISEED models, S11.g10191 and g10192, not yet in the table).
  Tested on reconstructions of the OdB3 loci (GPX_chr1a, SELENOW_chr1a) and on the other local tests; the read-through
  rule by a unit test (no OSKA2016 genome here). Page 2 was previewed with the expected models; OKI2018 SELENOO is
  mapped from the transcript onto an assembly that lacks 5 bases: its introns with other sites are where RefSeq has
  GT-AG or GA-AG introns (alignment placement around the missing bases).
- Gold-standard rerun (2026-10-24, commit 54fa4b4): OdB3 GPX_chr1a (183 aa, the two GenBank models joined), OdB3
  SELENOW_chr1a (spliced ESTs, 86 aa, GC-AG intron) and OSKA2016 SELENOW_chr2 (3' UTR cut at the next gene) as
  expected; OdB3 DIO_chr1 now from the two GenBank models joined (same protein; UTRs change); OSKA2016 SELENOO now
  lists the three ANISEED models that split the gene. Regression: OKI2018 SELENOO became the 70-aa C-terminal frame
  of the assembly. The RefSeq chain (587 aa), chosen before, is contradicted by all four transcripts, so the
  transcript chain read on the assembly was chosen; its frame is broken by the 5 missing bases and starts at an
  internal ATG, so it counted as complete and the transcript protein (600 aa, identity 0.79 against 0.80 on the
  shared 70 residues) was not used. New rule: at a similar identity, the complete transcript protein is also used
  when the assembly protein is less than half as long. From the tables, only OKI2018 SELENOO is in this situation.
  Also: annotation_status separates "codon inside the coding sequence, in another frame" (OdB3 GPX_PAR and SEPHS2)
  from codons in frame without declaration (OKI2018 GPX_chr1a and GPX_chr2, translated as X by RefSeq). Rerun needed;
  page 2 is committed after it (preview made with the expected OKI2018 SELENOO, the model of commit ce64671).
- Results section on the gold standard rewritten with the counts of the rerun (and the two expected changes above):
  53 models, all complete; 170 introns, 27 atypical (16%); GA-AG in all three assemblies for DIO_PAR, SELENOT,
  SELENOW_chr2 and TXNRD; public treatment of the codon per model: 3 declared, 2 in frame as X, 12 read as stop, 12
  skipped in an intron, 2 in another frame, 9 starting after, 3 ending before, 10 without coding sequence; genes split
  or joined by the annotations. Figure legend REANNOTATION added. There is no third page: the 18 genes fill two.
- Gold-standard rerun (2026-10-24, commit 0f885f8): OKI2018 SELENOO again mapped from the transcript (600 aa); OdB3
  GPX_PAR and SEPHS2 now "codon inside the coding sequence, in another frame"; no other change. The counts of the
  rewritten Results section are those of this table. Page 2 committed.
- Page layout (comments of Charles): the paragraph moved below the proteins, so that gene models and proteins form a
  grid; thin grey rulers between rows of genes, above the key and between columns.
- Step 3, first run of make secis-v2 (2026-10-24, commit 3dbd9f3, gold standard of that day): 111 positive windows
  (O. dioica 53 from the gold standard; O. longicauda 22, O. albicans 20, O. vanhoeffeni 16 from the stage 1
  predictions). 19 windows of the three other species have no hit with any model (10 of them shorter than 2 kb:
  contig ends); all O. dioica windows have hits. Pooled, species folds (each species held out from training):
  new model 85/111 above the training threshold, 0.43 evaluation negatives above it per Mb; merged 67, 0.43/Mb;
  current 82, 0.75/Mb (optimistic: trained on most of these sequences); RF00031 27, 0.11/Mb. Family folds: new 82
  (0.95/Mb), current 82 (0.75), merged 64 (0.50), RF00031 26 (0.02). In the species folds, the new model misses 2
  windows that the current one detects (OSKA2016 GPX_PAR, O. vanhoeffeni GPX.9) and detects 5 that it misses
  (OKI2018 SEPHS2, OSKA2016 and OdB3 SELENOT, OdB3 SELENOF and SEPHS2). Final model: 75 seeds (OKI2018 18, OSKA2016
  12, OdB3 8, O. longicauda 20, O. albicans 9, O. vanhoeffeni 8), 67 consensus positions. Problem: the best
  evaluation negative scores 44-48 bits with the new and current models in every family fold and in the O. dioica
  fold, more than most positives: probably an O. dioica random window that holds a SECIS element. The random windows were
  not kept away from the positive windows (a 2 kb window on the gene strand overlapping one of the 18 positive
  windows of an assembly is expected about 0.4 times per evaluation set of this size; the training windows, which
  set the thresholds, have the same problem). Changes: random windows that overlap a positive window are drawn
  again; new table stage2_secis_v2_negatives.tsv with the evaluation negatives above the training threshold and
  their coordinates (high-scoring ones that remain are candidate SECIS elements of selenoproteins not yet known).
  The rerun also takes the positive windows from the current gold standard (OdB3 GPX_chr1a: window after the stop
  codon of the joined model, 119 nt further; OSKA2016 SELENOW_pair2). Tested locally on the reduced genomes (no
  random window overlaps a positive window; negatives table written). Provisional conclusion: the new model is the
  candidate to replace the current one (more sensitive than the merged model at the same false-positive rate in
  the species folds, much more sensitive than RF00031); to confirm with the rerun, then choose its threshold.
- Reannotation pages reordered by family (2026-10-25, request of Charles), families kept together in reading order:
  page 1 GPX_PAR, GPX_chr1a, GPX_chr1b, GPX_chr2, DIO_PAR, DIO_chr1, TXNRD, SELENOO, SEPHS2; page 2 SELENOW_PAR,
  SELENOW_chr1a, SELENOW_chr2, SELENOW_pair1, SELENOW_pair2, SELENOT (Rdx family), SELENOF, SELENOK, SELENOS (ER
  membrane or lumen). New target make reannotation-pages with this order.
- Step 3, rerun of make secis-v2 (2026-10-25, commit cbc204c; random windows kept away from the positive windows):
  the 44-48 bit negatives are gone; the best evaluation negatives are now 19.7-24.3 bits (24.1 bits with the new
  model, a S. clava random window). Pooled, species folds: new 84/111 above the training threshold, all above the
  best evaluation negative of their fold, 3 negatives above the threshold (0.32/Mb; random windows of S. clava and
  C. savignyi, 21.6-23.1 bits); current 82 (0.64/Mb; optimistic); merged 66 (0.32/Mb); RF00031 27 (0.11/Mb).
  Family folds: new 81 (0.72/Mb), current 82 (0.50), merged 65 (0.25), RF00031 26 (0.02). O. dioica fold: new 49/53,
  no negative above 19.7 bits (current 45/53, 1 negative); misses OKI2018 and OSKA2016 SELENOO, OSKA2016 GPX_PAR and
  SELENOF. Held-out scores of the new model, species folds: >= 20 bits 84, >= 25 bits 64, >= 30 bits 43. Of the 27
  windows not detected, 19 have no hit with any model (9 of them shortened by a contig end). Final model: 76 seeds
  (OKI2018 18, OSKA2016 12, OdB3 9, O. longicauda 20, O. albicans 9, O. vanhoeffeni 8), CLEN 66, 7 rounds, best
  training background 17.8 bits in the last round; seeds start a median of 127 nt after the stop codon (quartiles
  67 and 238); the 53 O. dioica windows score 25.6-59.1 bits with it (training data).
  Decision: the new model (results/models/oikopleura_secis.cm) replaces appendicularia_secis.cm for the SECIS of
  genes; threshold 20 bits in the 3' windows of genes (3 of 4,691 evaluation negative windows of 2 kb passed it in
  the species folds). The threshold of genome-wide scans is left to step 4. Local check, not in the pipeline (6 Mb
  of OKI2018 random windows away from the positive windows, Infernal 1.1.5): 3 hits >= 20 bits (0.5/Mb); the best,
  28.1 bits at NC_137291.1:11878408 (-) and 27.5 bits at NC_137294.1:5395714 (-) (29.9 bits with the first model),
  have no stage 1 prediction within 20 kb: candidates for step 4. Downstream analyses that still use the first
  model (stage 1 SECIS figures, make secis-hits) to revisit in step 4.
  Figures and manuscript switched: reannotation pages and single-gene figures use the new model (open box below 20
  bits; none in the 53 models); new Results and Methods sections; the first-model section condensed to one
  paragraph (Figure SECIS-2 to redraw). While switching, the SELENOW pair showed two defects: (1) the SECIS of
  SELENOW_pair2 was that of SELENOW_pair1 (best hit of the 2 kb window, which reaches the next gene; own SECIS 23-92
  nt after the stop codon, 35.8-42.7 bits): the figure scripts now shorten the window at the next gold-standard
  coding sequence, as the training did; (2) in OKI2018 both models come from one transcript (TRINITY_DN2389) that
  holds the two genes, so that the 3' UTR of SELENOW_pair2 (892 nt, the maximum of the table) holds the coding
  sequence of SELENOW_pair1, and the 5' UTR of SELENOW_pair1 (492 nt) that of SELENOW_pair2. New rule in
  curation/gold_standard.py: UTRs are cut at the coding sequence of another gene of the table (3' UTR before the
  next one, 5' UTR after the previous one), with a note. Local test (reconstructed OKI2018 and OdB3 loci): only these
  two models change (OKI2018 SELENOW_pair2 3' UTR 197 nt, SELENOW_pair1 5' UTR 197 nt, one exon less each); no other
  gold-standard model reaches another gene's coding sequence. Expected after make gold-standard: 168 introns (the
  two dropped introns were counted twice), 3' UTR median 303 nt, range 4-695; the Results numbers are to update then.
- Reannotation pages, similarity to OdB3 (2026-10-25, request of Charles): right of each protein panel, a column with
  the percentage of residues similar to the OdB3 gold-standard protein (BLOSUM62 > 0, selenocysteine read as
  cysteine; global alignment with free end gaps, BLOSUM62, gap open 10, extension 0.5; denominator: the columns from
  the first to the last aligned pair, internal gaps counted as differences). OdB3 as reference: first published
  assembly, and OKI2018 is the outgroup of the three, so similarity to OKI2018 would miss the shortest distance.
  Values: human 44-70%, C. intestinalis 44-74%, OKI2018 88-99%, OSKA2016 94-100% (SELENOK: OKI2018 97%, OSKA2016
  94%). No value for SELENOO (no OdB3 model).
- make gold-standard rerun (2026-10-25, commit 7cd13a4): only the OKI2018 SELENOW pair changed, as in the local test
  (SELENOW_pair2: 3' UTR 197 nt, 2 exons; SELENOW_pair1: 5' UTR 197 nt, 2 exons; notes added). Table: 168 introns,
  27 atypical (16%); 3' UTRs median 303 nt, range 4-695. Manuscript numbers updated; figure inputs of the two genes
  regenerated (only the splice-site rows change); pages redrawn.
- Step 4 started (2026-10-25): stage2/uga_secis_screen.py, make secis-scan (heavy) and make uga-screen (light).
  Design sent to Charles before implementation. Genome-wide scans with the Oikopleura model of the six Oikopleura
  assemblies and Cint, Csav, Scla, Tdem, each with a copy shuffled in 1 kb windows (dinucleotides preserved); the
  genome-scan threshold is the lowest integer score with at most one expected shuffled hit in the largest
  Oikopleura assembly (pooled rate). Transcript screen: frames through a TGA (ATG, or 5' open, to a stop codon that
  ends before a SECIS hit of 20 bits or more), at least 40 codons, internal (>= 10 codons after the TGA) or
  C-terminal; frames that overlap a longer frame, or a longer ordinary ORF, in another reading frame are left out.
  Positive controls: the gold standard in each O. dioica assembly; human and C. intestinalis RefSeq mRNAs (Oikopleura
  model and RF00031; known = transl_except selenocysteines). Genome-based screen: the first version (a TGA in an open
  stretch of the 5 kb before a SECIS hit, longest stretch kept) picked random open stretches rather than the
  selenocysteine exon in the local test (SELENOW_chr1a, SELENOF, GPX_chr1b missed), so it was replaced by a
  comparative screen (as SGP2, Castellano et al. 2004): tblastx between species of the 5 kb upstream of each
  genome-scan hit; a TGA aligned to a TGA (E <= 1e-5, >= 10 aligned residues on one side, no stronger alignment of the
  same regions in another frame, open stretch >= 20 codons). Candidates annotated by blastp against the human RefSeq
  proteome (residue aligned to U) and the gold-standard proteins; O. dioica candidates grouped across assemblies.
  Local test (reduced OKI2018 regions with 7 gold-standard genes, RefSeq-like annotation and Trinity-like
  transcripts from the gold standard; OdB3 contigs with 5 genes standing in for another species): transcript
  screen 7/7 in both sets (rank 1 each; 3 other frames from the 3' extension of SELENOW_pair2 into SELENOW_pair1);
  comparative screen 5/5 in each genome, with 9 other candidates each, expected from identical sequences of the same
  species (intergenic TGAs align); the step-3 locus NC_137291.1:11878408 scores 28.1 bits in the genome scan.
  Next: run make secis-scan, then make uga-screen, on Deigo.
- Step 4, held-out genomes (2026-10-25, request of Charles): B. stygius, M. erythrocephalus and F. borealis are
  scanned too (and their shuffled copies; the genome-scan threshold still comes from the Oikopleura copies only) and
  screened by the comparative screen against all species; they are not partners of the Oikopleura assemblies, so the
  Oikopleura candidates do not depend on them (decision of 2026-10-24). Reference genes of a held-out genome: its
  stage 1 selenocysteine predictions (all left out as held out). Tested with the OdB3 test contigs as a held-out
  stand-in: OKI2018 no comparative candidate (its only partner was held out), the stand-in 14 with the 5 genes.
- Step 4, first make uga-screen on Deigo (2026-10-25): formatdb failed on the regions of the comparative screen
  (probably empty regions: SECIS hits at contig starts in fragmented assemblies), and the single genome threshold
  came out at 33 bits, because it was set for the largest assembly (one expected false hit in it). Changes: regions
  with fewer than 60 bases (A, C, G, T) left out; failed commands now print their error output (and formatdb.log);
  thresholds per genome from the pooled shuffled rate: genome threshold (at most one expected false hit in that
  genome; genome-wide SECIS predictions) and comparative threshold (at most 100 expected false hits, not below
  15 bits: the regions of the comparative screen, where the second species gives the specificity); the summary
  table gives both, with the hits and the expected false hits at each; the hits table lists the hits at the
  comparative threshold, flagged when at the genome threshold.
- Step 4, first complete make uga-screen (2026-10-25, commit 50c6990). Thresholds per genome (pooled shuffled rate):
  genome threshold 28-29 bits in the O. dioica assemblies, 31-34 bits elsewhere; comparative threshold 18-25 bits.
  Assembly sizes (both strands): O. dioica 113-141 Mb; O. longicauda 618, O. albicans 732, O. vanhoeffeni 1,287,
  B. stygius 793, M. erythrocephalus 1,748, F. borealis 286, Cint 281, Csav 401, Scla 754, Tdem 1,604. Genome hits at
  the genome threshold: OKI 23, OSKA 22, OdB3 29, Olon 30, Oalb 15, Ovan 19 (expected false hits 0.4-0.9); the
  non-target tunicates are at the level of their shuffled copies (Cint 267 hits at 15 bits against 344 shuffled):
  the Oikopleura model does not see ascidian SECIS elements genome-wide.
  - Transcript screen: all 18 gold-standard genes found in the OKI2018 RefSeq and Trinity sets, 15 in the OSKA2016
    transcripts, 7 in the ANISEED models (models missing or broken), OdB3 GenBank and ESTs most of them. Controls:
    human RefSeq 18 of 25 selenoprotein genes with RF00031, 3 with the Oikopleura model; C. intestinalis 15 of 17 with
    RF00031, 9 with the Oikopleura model.
  - Comparative screen: 14-22 known genes per O. dioica assembly, 20 of 28 O. longicauda, 12 of 20 O. albicans, 9 of 16
    O. vanhoeffeni stage 1 selenocysteine predictions; few other candidates (OKI 2, OSKA 3, OdB3 1, Olon 1, Oalb 5,
    Ovan 6). Among them SELENOK of O. albicans, O. vanhoeffeni and M. erythrocephalus (C-terminal ...GGUGR, partners
    the O. dioica SELENOK and each other), the gene that stage 1 missed in every appendicularian; TXNRD of
    M. erythrocephalus; SELENOW_PAR, SELENOW_chr1a and GPX of O. albicans or O. vanhoeffeni not in stage 1.
  - New candidate family: a 68-70 aa protein with a hydrophobic segment and the C-terminus ...RRKF(T/S)NWKPGKKSU,
    TGA followed by a stop codon, in all three O. dioica assemblies (RefSeq XM_081776455.1 on NC_137292.1 at
    12254043, ANISEED S22.g13919, GenBank GSOID_T00008777001, Trinity, OSKA2016 transcripts, 7 ESTs; SECIS 28-37 bits)
    and in O. albicans (genome, SCLG01000111.1, ...RRRFCNWTPGKKSU, SECIS 31.9 bits); no human hit.
  - Held-out genomes: a family of short open stretches (...UIA(N/D/S)HT(K)LL(P/S)ETTSSGKLSLKN...) on five
    M. erythrocephalus contigs, partnered with a B. stygius contig start: a repeated element, not a gene.
  - Defects found and fixed in stage2/uga_secis_screen.py: the spliced leader of O. dioica mRNAs holds a TGA in one
    frame (...HF-U-VRFRLSN), which gave 5' fragment frames with the U at codons 2-6 (now at least --min-before 10
    codons before the TGA; the earliest known selenocysteine is codon 13); other frames of a transcript that holds a
    known gene, and other TGAs of gold-standard genes (genome overlap, or protein identity >= 90 % to a gold protein of
    any O. dioica assembly), are no longer counted as candidates (column status); the summary counted haplotype
    copies and stage 1 names as recovered genes (now gold-standard genes, stage 1 selenocysteine predictions or
    annotated selenocysteines, as fits the genome); human hits now carry their RefSeq description; new column
    other_species (candidates of other species whose U aligns to the U). To rerun: make uga-screen (minutes).
- Step 4, second make uga-screen (2026-10-25, commit 0b4b4ef). Gold-standard genes found by the transcript screen:
  OKI2018 Trinity 18/18, RefSeq 15/18 (no RefSeq model for SELENOK and SELENOW_chr1a; the RefSeq SELENOS frame is
  lost where the assembly lacks 13 bases), OSKA2016 transcripts 15 (no DIO_PAR, DIO_chr1, GPX_chr1a), ANISEED 9, OdB3
  ESTs 15 (no GPX_chr1b, SELENOW_pair2, SEPHS2), GenBank 12; comparative screen 14-15 per assembly (SELENOO, SELENOS,
  SELENOT, DIO_chr1 lack a partner in another species at the comparative threshold). The first summary's 18/18 for
  RefSeq counted stage 1 names. The cross-species column missed the KKSU family: the U is its last residue and was
  read as X, so local alignments stopped before it; U is now read as C in these comparisons (checked on the table:
  the O. dioica KKSU candidates and the O. albicans one now link). New status "homolog of a known selenoprotein"
  (partner at a known selenocysteine codon, or U aligned to the U of a gold-standard or human selenoprotein), so
  that "candidate" holds only proteins of no known selenoprotein family; --known-identity lowered to 80 % (an OSKA2016
  SELENOW_pair1 transcript at 89 % was a candidate).
- Step 4, third make uga-screen (2026-10-25, commit d2c5bfb): candidates supported in another species (U aligned to
  U) are only the KKSU family (all O. dioica sets and O. albicans) and the M. erythrocephalus/B. stygius repeat. The
  transcript screen also finds TXNL1 (thioredoxin-like protein 1, best human hit NP_004777.1, 40-41 % identity) in
  all three O. dioica assemblies at the codon of curation/candidates.tsv (OKI2018 NC_137292.1:13041527, OSKA2016
  S27:260385, OdB3 FN653019.1:362479): a 168-codon frame ending ...KKVGUG and a stop codon, SECIS 29.0-31.8 bits in
  the annotation transcripts (23.6-31.8 bits in the transcript assemblies), 67-95 nt after the stop codon. By the
  rule fixed in advance for step 5 (SECIS above background in the three assemblies, TGA conserved), it qualifies;
  its other species are to check in step 5. A weaker case: CPLANE2-like (best human hit NP_112169.2, ciliogenesis and
  planar polarity effector 2), C-terminal TGA, SECIS 37.6-42.8 bits, in OSKA2016 and OdB3 only.
  Classification changes: homolog status propagated through the partner of the comparative screen (SELENOK of
  O. vanhoeffeni through that of O. albicans); comparative candidates upstream of the SECIS of a known selenoprotein
  gene (its coding sequence ends at most 3 kb before the hit; for example a TGA of the gene upstream of OSKA2016
  GPX_chr1a, whose SECIS scores 59.1 bits) get the status "upstream of the SECIS of a known gene".
- Step 4 done (2026-10-26, fourth make uga-screen, commit 6243f71): classification as intended (comparative screen:
  candidates are the KKSU family in OKI2018, OSKA2016, OdB3 and O. albicans, and the M. erythrocephalus/B. stygius
  repeat; ten copies of known selenoproteins absent from stage 1; four TGAs upstream of the SECIS of known genes).
  Genome-scan hits at the genome threshold after no known gene in OKI2018: the KKSU SECIS (NC_137292.1:12253777,
  36.6 bits), NC_137292.1:13040099 (31.6 bits, 1.4 kb upstream of the TXNL1 codon; similar hits upstream of TXNL1 in
  OSKA2016 and OdB3; the TXNL1 SECIS itself, 70 nt after its stop codon, lies after a stage 1 Ostsp2 prediction
  labelled unaligned), NC_137295.1:2369276 (30.1 bits) and the step-3 locus NC_137291.1:11878408 (28.1 bits).
  Provisional Results and Methods written (section "UGA/SECIS screens outside the known families"). Next: step 5
  (TXNL1; the KKSU family in the other species and other databases).
- Step 5 started (2026-10-26): stage2/candidate_genes.py (make candidate-genes). Candidates in curation/candidates.tsv
  (TXNL1; KKSU, the family of step 4, provisional name after its C-terminus, rows for OKI2018, OSKA2016, OdB3 and
  O. albicans) and curation/candidates.faa (the OKI2018 frames of the screen, and the O. albicans KKSU open stretch).
  (1) The rule fixed in advance for TXNL1, applied to both: TGA at the curated position and a SECIS of at least
  20 bits (Oikopleura model, 2 kb after the next stop codon) in the three O. dioica assemblies -> counted, otherwise
  candidate; the genome thresholds of step 4 given for comparison. (2) tblastn against the six Oikopleura assemblies,
  the three held-out appendicularians, the ten other tunicates and the two amphioxus: best locus, codon aligned to the
  U (read beyond the alignment end when it stops at most 5 codons before the U, flagged), codons to the next stop,
  SECIS after it (Oikopleura model and RF00031). (3) blastp against the RefSeq proteomes (vertebrates, tunicates,
  amphioxus): residue aligned to the U and length after it. Local test (OKI2018 chromosome 2 and OdB3 FN653019.1):
  contexts of candidates.tsv reproduced (MSELKKVGUG*APS; KKSU NWKPGKKSU*MLY); SECIS OKI2018 TXNL1 31.0 bits, OdB3
  29.0, OKI2018 KKSU 36.6.
- Step 5, make candidate-genes on Deigo (2026-10-26, commits 09bdc95 and 4b13c80).
  - Rule (fixed in advance for TXNL1, applied to KKSU): both counted. TXNL1: TGA, then one codon and a stop codon,
    in OKI2018, OSKA2016 and OdB3 (context MSELKKVGUG*APS / MNELKKVGUG*APS as curated); SECIS 31.0, 30.7, 29.0 bits
    after the stop codon (genome thresholds of step 4: 28, 28, 29). KKSU: TGA followed directly by a stop codon in the
    three assemblies; SECIS 36.6, 31.3, 30.8 bits.
  - KKSU in other genomes (tblastn): found in all six Oikopleura assemblies and in B. stygius and M. erythrocephalus
    (49-63 % identity to the OKI2018 protein), always with the TGA followed directly by a stop codon (TAA; TAG in
    M. erythrocephalus); SECIS after it in the six Oikopleura assemblies (O. longicauda 36.5, O. albicans 31.9,
    O. vanhoeffeni 21.6 bits; RF00031 17.4 at most), none after the B. stygius and M. erythrocephalus hits; no hit in
    F. borealis, the ten other tunicates, amphioxus, or any RefSeq proteome (human included). The stop codon right
    after the TGA is conserved across species whose proteins are only 50-60 % identical, as expected if the TGA is read
    and the next codon ends the protein.
  - TXNL1 in other genomes: no TGA at the position of the O. dioica selenocysteine. The codon just after the end of the
    alignment is a glycine codon in O. albicans, O. vanhoeffeni, M. erythrocephalus, F. borealis and the three salps
    (read in frame beyond the alignment, which stops just before that position; the 6 codons that follow to a stop
    may be intron in species whose protein continues); the position is not aligned in O. longicauda, B. stygius, the
    ascidians and B. lanceolatum; no SECIS after any of these loci. RefSeq TXNL1 of human, mouse, zebrafish, ascidians,
    Thalia and B. lanceolatum: 285-293 aa, with the C-terminal PITH domain; the O. dioica C-terminus does not align
    with them in blastp. Correction (2026-10-26, from the reannotation plot): the O. dioica protein (168 aa) aligns with
    the PITH domain of human TXNL1 (residues 120-283), not with its thioredoxin domain, which it lacks; its
    selenocysteine stands where human and C. intestinalis TXNL1 have the glycine of their C-terminus ...V(G)KKGESH
    (human G283). The 6 codons to a stop codon read after the glycine codon in O. albicans, O. vanhoeffeni,
    M. erythrocephalus, F. borealis and the salps are this KKGESH end, not intron: in those species the TGA position
    holds a glycine and the protein ends 6 residues later, as in human.
  - Defect: tblastn found no hit at all in S. clava, B. schlosseri, P. mammillata and B. floridae, although blastp finds
    TXNL1 in the S. clava and B. schlosseri proteomes (E 1e-38, 8e-31): to diagnose (their cached rows.json must be
    removed before a rerun).
- Baseline decision on TXNL1 and KKSU deferred (2026-10-26, Charles: see their reannotation plots first). Prepared:
  curation/baseline_candidates.tsv (the two genes in the format of the baseline, status "candidate") and make
  gold-candidates (curation/gold_standard.py on them, outputs in curation/gold_candidates and
  results/tables/gold_candidates.tsv, the baseline and its gold standard untouched); figures/reannotation_inputs.py
  takes no UniProt reference if none exists (KKSU), and figures/reannotation_page.py then uses the first
  gold-standard protein as the axis of the protein panel.
- Reannotation plots of the candidates (2026-10-26, gold-candidates run of Charles, commit 4b3445d): TXNL1 168 aa in the
  three assemblies (4 coding exons; public models with the same exon chain end at the TGA read as stop); KKSU 68-70 aa
  (3 exons; ANISEED and GenBank models end at the TGA read as stop; in OKI2018 the RefSeq mRNA of XP_081632581.1, a
  173 aa protein, holds KKSU in its 3' UTR, as in an operon). figures/reannotation_inputs_candidates.sh builds their
  inputs (--table option added to figures/reannotation_inputs.py); figures/reannotation_page.py places a C-terminal
  selenocysteine a few residues beyond the end of the local alignment, so that it is drawn (TXNL1: at human G283).
  curation/gold_candidates/*.fa were ignored by git (needed for the SECIS boxes): .gitignore exception added.
- Decision (2026-10-26, Charles, after the plots): TXNL1 and KKSU join the baseline (curation/baseline_sec_genes.tsv,
  status sec, note on how they entered), so that all computations include them for now; they can be removed later if
  doubts remain. In figures, same style but their own figure (make reannotation-pages: reannotation_candidates), to
  signal that they need the reader's attention. The separate candidate run (baseline_candidates.tsv, make
  gold-candidates, curation/gold_candidates) is removed; their figure inputs join figures/reannotation_inputs.sh. To do:
  make gold-standard (20 genes per assembly), then the figure, the manuscript counts, and make uga-screen (they become
  known genes). Not rerun now: the SECIS model of step 3 (trained on the 18 genes); a later make secis-v2 would include
  their windows.
- Gold standard with 20 genes (2026-10-26, commit 3fb327e, run by Charles): only the six TXNL1/KKSU rows added, the
  54 earlier rows unchanged. 60 codons, 59 models, all complete; 57 with a 3' UTR (median 303, 4-695); 184 introns,
  27 atypical (14.7%); TGA read as stop in 17 public models (OKI 3, OSKA 5, OdB3 9); 11 without annotated CDS (OKI
  KKSU added: in the 3' part of RefSeq XM_081776455.1, mRNA of a RING finger protein 11 homolog, XP_081632581.1).
  Figure CANDIDATES drawn (make reannotation-pages: figures/reannotation_candidates.*; new options of
  reannotation_page.py: --fit, key right below the genes and page cut below it, 125 mm; --sources, last line of the
  key). Manuscript: gold-standard counts updated (the SECIS model section now says it was trained on the 53 models
  of the 18 known genes); step 4 text corrected (KKSU is not annotated in RefSeq, it lies in the 3' UTR of another
  mRNA); new Results section on the two genes (rule, TXNL1 PITH domain and the Gly of other species, KKSU in the
  Oikopleuridae) with the legend of Figure CANDIDATES and Methods "Candidate genes". To do: make uga-screen
  (step 4 numbers), then update the step 4 section; diagnose the tblastn genomes without hits.
- Step 4 rerun with the 20-gene gold standard (2026-10-26, commit 837c4cd): TXNL1 and KKSU found in all six O. dioica
  transcript sets; KKSU by the comparative screen in the three assemblies (TXNL1 has no partner TGA); the 18 known
  genes as before. Genome-scan hits after gold-standard genes: 20, 18, 18; OKI2018 hits after none: three (1.4 kb
  upstream of the TXNL1 codon, NC_137295.1, NC_137291.1). The only candidates supported in another species are the
  B. stygius / M. erythrocephalus repeat. Manuscript step 4 section updated. make gold-standard rerun: unchanged.
- tblastn without hits in Scla, Bsch, Pmam, Bflo (2026-10-27): the databases are complete (formatdb 380, 220, 311 and
  432 sequences; Statistics_db-len, in residues for tblastn, x3 = 377, 496, 240 and 513 Mb; the same ratio to the
  gzipped file size as Blan and Tdem). No technical failure: no locus at E <= 1e-5 (a 168-aa query against 377 Mb
  needs about 51 bits per HSP). Written into the step 5 text; the proteomes are the evidence there.
  candidate_genes.py: the proteome residue aligned to the U is also read at the same offset beyond the end of the
  local alignment when it stops fewer than --extend (5) residues before the U (column residue_note "extended"), as in
  the figure (human TXNL1 G283). To do: make candidate-genes (the genome searches are cached).
- SECIS model description (2026-10-27; stage2/secis_model_description.py, figures/secis_model_figures.py, make
  secis-description; repository files only). Seeds: 76 from 10 families and 6 species; median identity 0.45 over the
  consensus positions; 25 of 2,850 pairs >= 0.7 (23 O. dioica orthologs, 0.72-0.91; two GPX paralog pairs in Oalb
  and Ovan); 54 groups at 0.7, largest 3 (training already dropped seeds >= 90% identical). Structure: helix 1 (5 bp),
  internal loop with UGA (5') and GA (3') in all seeds, helix 2 (10 bp), AAA in a 5-nt loop below helix 3 (3 bp) in
  all seeds: the form 2 arrangement (Grundner-Culemann 1999). Pairing in the seed alignment (83/83/94% canonical,
  MI q < 0.05 at 18/18 pairs) is matched by decoys aligned to the same model (83/78/88%, 12/18 pairs): the alignment
  makes these statistics, so R-scape on the seeds would not test conservation (non-homologous sequences aligned to
  the model); R-scape was not run (building it here was refused; not needed). Evidence instead from the orthologous
  SECIS of the three O. dioica assemblies aligned by sequence: single changes from a canonical pair kept it 18/27
  (expected 0.29 with the transition rate of unpaired positions; p = 5.5e-5), double changes 10/11 (expected 0.43;
  p = 0.0015); core motif positions never changed (152 comparisons), paired 18% and other unpaired 28% of positions
  differ; same with OSKA or OdB3 as reference. Genome scans at the genome threshold, by class: O. dioica 10 of 74 hits
  without a gene upstream (28.1-36.9 bits; 2.4 expected), 4 of them 0.8-1.4 kb upstream of the TXNL1 codon in the
  three assemblies; other Oikopleura 21 of 64 (33.0-46.7 bits; 1.8 expected), some in pairs as from repeats; held-out
  7 of 12; non-target tunicates 3 of 12 (the other 9 after stage 1 selenoproteins: the step 4 sentence "the model does
  not recognise" corrected). Controls: human selenoprotein mRNAs up to 31.5 bits, others up to 22.8; C. intestinalis
  37.3 and 28.3. Manuscript: two paragraphs on the model, Figures SECIS-MODEL and SECIS-SCORES, Methods "Description of
  the SECIS model"; references Walczak1996, GrundnerCulemann1999, Rivas2017.
- TXNL1 and KKSU SECIS in training? Assessment (2026-10-27): they have the motifs and helix 3 of the seeds, score
  29.0-36.6 bits (lower half of the training windows, 25.6-59.1), and are as distant from the seeds (nearest 0.61-0.68)
  as seeds of different genes (median 0.62): they would add two families (6 of 117 windows). Against: both were
  found by this model (self-reinforcement; their family folds would not be independent), and retraining changes
  every model-dependent number (steps 3 to 5). Proposal: not now; include them in one final rebuild if stage 2 adds
  genes. To examine: the 31 Oikopleura genome-threshold hits without a gene upstream (results/tables/
  stage2_secis_scan_classes.tsv): repeats, missed selenoprotein genes, other frames; and the conserved hits upstream
  of TXNL1.
- Figures SECIS-MODEL and SECIS-SCORES revised (2026-10-27, Charles's comments): model figure 2x2 (A structure with
  a sequence logo of the seeds and the relative entropy of the model emissions per position, mean 0.740 bits as
  cmstat; D on its right; B and C below); B labels give 76 seeds and the 38 O. dioica seeds of 17 genes with
  orthologs; scores panel B with Okabe-Ito colours, one marker per class, transparency for the two frequent classes,
  less jitter, the genome-threshold tick in the legend. The human RF00031 controls classed as other mRNAs above 38
  bits are GPX1 (NM_001329502.2) and DIO1 (NM_001039715.3, NM_001324316.2) transcripts whose frame does not end at an
  annotated selenocysteine: their hit is the SECIS of the gene (text added). Decision (Charles): the TXNL1 and KKSU
  SECIS stay out of the model. New: stage2/secis_orphan_hits.py (make secis-orphans) for the 48 hits of the classes
  "nothing found upstream" and "stage 1 prediction without selenocysteine". Literature round 3: openalex_round3.sh
  (outside the repository).
- Hits without a selenoprotein gene (2026-10-27, make secis-orphans, two runs; the second reports the best human and
  the best gold-standard blastx hit). O. dioica (10): 2 copies of the DIO_chr1 SECIS on extra OdB3 contigs, 2 repeats
  (25 and 53 copies), 4 = one element conserved in the three assemblies between a MARS1-like gene and TXNL1 (in the
  5' UTR of the OSKA TXNL1 transcript), OSKA S54 conserved (26.5 and 27.0 bits in OKI and OdB3), OKI NC_137295.1
  without counterpart. Other Oikopleura (26): 8 copies of known families by blastx (DIO, GPX, SELENOK, SELENOT,
  SELENOW; all five stage 1 predictions without selenocysteine), among them an O. longicauda SELENOK (TGA aligned to
  the selenocysteine, 69 nt before the hit, after a CPLANE2-like region as in O. dioica); 1 repeat; 4 O. albicans hits
  after an OGG1-like region on two haplotype contigs; 13 without protein hit. Held-out: 5 of 9 copies of GPX, SELENOT,
  SELENOW. Figure SECIS-SCORES B: final classes (copy of a known selenoprotein 19, repeated sequence 3, nothing found
  30) in stage2_secis_scan_classes.tsv (final_class, evidence); manuscript paragraph rewritten. To consider: the
  MARS1-TXNL1 and S54 elements; O. longicauda SELENOK for the census (step 9).
- Figure SECIS-SCORES B (2026-10-27, Charles): genome thresholds kept strict (one expected false hit per genome);
  hits below them, down to the comparative threshold, drawn small (orange after a selenoprotein gene) to show the
  background; legend below the panel. Figures embedded in manuscript/manuscript.md (PNG, SVG for the scheme and the
  tube map): PIPELINE (figures/pipeline_metro.svg copied from the pipeline repository), TREES, SECIS-1,
  REANNOTATION, SECIS-MODEL, SECIS-SCORES, CANDIDATES; legends added for PIPELINE, TREES and SECIS-1.
- Discussion, provisional paragraph "What the *Oikopleura* SECIS model has learned" (2026-10-27, Charles approved).
  Kept for later at his request: cite earlier SECIS work across lineages in that paragraph; test the model on
  *F. borealis* (held out, outside *Oikopleura*).
- Literature, third round (2026-10-27; OpenAlex script of round 3, run by Charles; 1,316 works; abstracts read
  outside the repository; metadata from PubMed): see manuscript/overlap_check.md. 22 entries added to
  manuscript/references.bib; [REF] placeholders filled for C20/C60 (Le2008), slow-fast (Brinkmann1999) and the AU
  test (Shimodaira2002); the operon example stays open. New citations in the Introduction, the orthology of the
  SECIS-binding proteins, the SECIS model description, step 4, step 5 and the Discussion. Lead for the Discussion:
  the human TXNL1 C-terminal tail (six residues after Gly283, ending with His289) coordinates the RPN11 zinc in the
  degrading proteasome (Arkinson2025, Lee2026); the *O. dioica* selenocysteine sits at the position of Gly283,
  followed by one residue, and the protein lacks the thioredoxin domain.
- TXNL1 homologs (2026-10-27, make txnl1-homologs, run by Charles; results/tables/stage2_txnl1_homologs*.tsv). The
  selenocysteine gene is not the TXNL1 ortholog. Each assembly also has a full-length TXNL1 (thioredoxin domain + PITH,
  285 residues: OKI XP_081633344.1, NC_137292.1 at 10.70 Mb, minus strand, 2.34 Mb from the selenocysteine gene;
  OdB3 CBY22771.1, FN653019.1, 1.02 Mb from it; OSKA S14.g11310, a 756-residue model whose best human hit is GPR107,
  probably a GPR107-like gene merged with TXNL1). Human TXNL1 and the full-length gene are reciprocal best hits
  (OKI and OdB3: 205-206 bits, 40% identity over residues 9-287; back to TXNL1 at 205 bits, TXN 91). The
  selenocysteine gene hits human TXNL1 first (132-134 bits, margin 79-82 over PITHD1), but its best hit in its own
  assembly is the full-length gene (64-68% identity over residues 7-166, against about 40% with human TXNL1); next, a
  PITHD1 homolog (33-34%). Thioredoxin part of human TXNL1 (1-117): best protein hit the full-length gene in OKI
  (69 bits), a DNAJB4-like protein in OSKA and OdB3 (72 bits; full-length TXNL1 second, 68-71 bits); the domain alone
  does not separate TXNL1 from TXN (blastx of the thioredoxin exons: TXN 61, TXNDC2 57). No thioredoxin-domain hit
  near the selenocysteine gene. Reading: a PITH-only paralog of TXNL1, from a duplication in the tunicate or
  appendicularian lineage. To do: name; whether other appendicularians have the PITH-only copy (step 5 compared the
  best locus per genome and RefSeq TXNL1 proteins, i.e. full-length orthologs); residue of the full-length O. dioica
  TXNL1 at the selenocysteine position; PITH tree; the OSKA merged model as a possible example of merged genes.
  Genome-locus grouping of txnl1_homologs.py splits or merges some loci (overlapping HSPs across short introns);
  the protein searches are the reliable part.
- make txnl1-paralogs written (2026-10-27; stage2/txnl1_paralogs.py, tested on synthetic data): the full-length TXNL1
  and the selenocysteine copy in the other genomes and proteomes (structure, codon at the selenocysteine position,
  SECIS), a PITH tree rooted on PITHD1, and the neighbourhood of the full-length copy (GPR107-like gene, merged
  OSKA2016 model, transcripts spanning both). Gene symbol of the selenocysteine copy: Charles proposes PITHD1L1;
  to settle after the tree (a symbol based on PITHD1 fits if the copy groups with PITHD1, one based on TXNL1 if it
  groups with TXNL1, as its best hits suggest).
- TXNL1 paralogs, results (2026-10-27, make txnl1-paralogs, two runs by Charles). Residue aligned to the
  selenocysteine: glycine followed by six residues in the full-length O. dioica TXNL1 (OSKA2016 and OdB3 ...GEAGEAH,
  OKI2018 ...GEVGERE) and in human TXNL1 (...GKKGESH, His289). Other genomes: a full-length TXNL1 with a glycine codon
  and six codons in O. albicans, O. vanhoeffeni, M. erythrocephalus, F. borealis and the three salps (tryptophan in
  H. roretzi); no PITH-only copy outside O. dioica, except a 123-residue PITH-only TXNL1 protein in C. intestinalis
  (XP_078490523.1, LOC144746873; glycine and six residues; groups with Ciona TXNL1, 96), absent from C. robusta.
  The O. longicauda PITH-only fragments (contig ends) group with Oikopleuridae TXNL1, not with the selenocysteine copy.
  PITH tree (rooted on PITHD1, 15 sequences, including the O. dioica PITHD1 homolog XP_081636562.1): the
  selenocysteine copy is in the TXNL1 clade, sister to the O. dioica full-length TXNL1 (92), the pair sister to the
  TXNL1 of the other Oikopleuridae (95), F. borealis outside (appendicularian clade 100): a duplication in the
  O. dioica lineage. Neighbourhood: a GPR107-like gene upstream of the full-length TXNL1 on the same strand (579, 610
  and 2,152 nt between the hits in OKI2018, OSKA2016 and OdB3), two models in OKI2018 and OdB3, one merged model in
  OSKA2016, no transcript aligned to both. The supports drawn in the first figure were read from nodes after
  rerooting, which moves them along the path to the new root; corrected (supports read by bipartition; the values
  above are off that path). Symbol: PITHD1L1 would place the gene in the PITHD1 family, which the tree contradicts;
  to discuss with Charles.
- Names (2026-10-27, Charles): TXNL1.a for the full-length O. dioica TXNL1 and TXNL1.b for the selenocysteine copy,
  after the tunicate nomenclature guidelines (Stolfi et al. 2015: paralogs sharing one human ortholog take .a, .b;
  "-like" discouraged). Capitalisation of O. dioica gene symbols deferred. Applied to the manuscript (step 5 rewritten
  with the results of make txnl1-homologs and make txnl1-paralogs; Figure TXNL1-TREE and Methods "TXNL1 homologs and
  paralogs" added) and to Figure CANDIDATES (figures/reannotation_page.py --labels TXNL1=TXNL1.b). The internal
  identifier stays "TXNL1" in curation/ (baseline, candidates, gold standard), the scripts and the result tables until
  the capitalisation is settled, so that both changes are made in one pass. Corrections made on the way: the
  P. mammillata and B. floridae proteomes hold no TXNL1 (only a PITHD1 homolog; the earlier text said they held TXNL1);
  TXNL1.b has no CxxC motif.
- Renaming done (2026-10-27, Charles: "go ahead"): the internal identifier of the selenocysteine gene is now TXNL1.b in
  curation/ (baseline, candidates, gold standard: identifiers and gene column only; the family stays TXNL1), in
  results/tables/gold_standard.tsv (gene column), in the figure inputs (reannotation_TXNL1.b.*) and in the script
  defaults. A rerun of make gold-standard must leave curation/gold_standard/ unchanged (check). The other tables
  change at the reruns of uga-screen, secis-description, candidate-genes, txnl1-homologs and txnl1-paralogs; the cache
  of candidate-genes is redone when the candidate names change. Figure CANDIDATES: TXNL1.a (public model aligned on
  its codon aligned to the selenocysteine, and protein) drawn with TXNL1.b, from
  figures/inputs/reannotation_TXNL1.b.paralogs.{tsv,fa} written by make txnl1-paralogs. Tree supports checked on the
  rerun: 92, 95, 100 and 97 (Ciona TXNL1 with the C. intestinalis PITH-only protein).
- Reruns after the renaming (2026-10-27, Charles): curation/gold_standard reproduced unchanged; the other tables changed
  in the names only; Figure CANDIDATES drawn with TXNL1.a. Introns (make txnl1-introns, new; repository files only): the
  three coding introns of TXNL1.b (phase 0) are at the positions of three of the five coding introns of TXNL1.a, in the
  three assemblies; TXNL1.a has one more intron in the shared region (phase 1, codon 157) and one in the thioredoxin
  part (phase 2, codon 67): a duplication of genomic DNA, not a retrocopy. Coincidence question (Charles): selenoproteins
  are mostly oxidoreductases, many with a thioredoxin fold, but the TXNL1.b selenocysteine replaces a glycine at the
  start of the C-terminal tail, not a redox cysteine, and the thioredoxin domain is lost; possible Discussion item.
- Position of the selenocysteine codon (2026-10-27, question of Charles: are TGA codons near exon boundaries more than
  by chance? make sec-positions, new; repository files only). OSKA2016, 20 genes (18 with coding introns): no excess
  near junctions (mean mid-rank 0.58 against 0.5, one-sided p = 0.89; 5 codons within 30 nt of a junction against 6.2
  expected; no codon split by an intron, 0.2 expected); the two closest are SEPHS2 (1 nt) and SELENOT (2 nt). Strong
  excess at the C-terminus: 6 of 20 genes have the codon among the last 3 before the stop (KKSU, TXNL1.b, TXNRD,
  SELENOS 1-2; SELENOK, SELENOO 3), against 0.4 expected (P = 2e-6; one gene per family, 6 of 12, P = 2e-8); test
  chosen after looking, and four of the six families have a C-terminal selenocysteine in other animals too. Both new
  genes are C-terminal: KKSU (TGA followed directly by a stop codon) fits Charles's scenario (an ancestral TGA stop
  becomes a Sec codon when a SECIS is gained); TXNL1.b's TGA replaces the GGA glycine codon of TXNL1.a (one
  substitution). Possible test on Deigo: tandem stop codons (TGA followed by a stop) at the end of O. dioica genes,
  and SECIS-like hits after them. Not in the manuscript yet.
- Tandem stop codons (2026-10-28, make tandem-stops, run by Charles). In the three O. dioica assemblies, 32-36% of
  genes end with TGA; 535 to 700 genes per assembly (3.5-4%) end with a TGA followed by another stop codon at +1 or +2,
  as often as the background frequency of stop codons predicts (enrichment at +1 after TGA 0.94-1.02; after TAA and TAG
  0.45-0.83). SECIS-like hits of at least 20 bits in the 2 kb after these genes: 0, 1 and 1 (0-0.17%), as after TAA or
  TAG genes (0.07-0.09%; Fisher p 0.39-1); known selenoprotein genes: 67-88%. The only genes of the class with a hit
  at the genome threshold are known C-terminal selenoproteins whose selenocysteine was annotated as a stop (TXNRD,
  SELENOS, KKSU, TXNL1.b; a SELENOW pair copy on an extra OdB3 contig; S. clava TXNRD2): no new candidate. The path is
  available in hundreds of genes, but leaves no excess of SECIS-like elements after them (limited power: such hits are
  rare after any gene). Side observation: relative to TAA and TAG, TGA is more often followed by a stop at +1 in
  tunicates and Branchiostoma, not in vertebrates (not tested further).
- Stage 2 prepared (2026-10-28, steps 6 and 7; request of Charles: profiles from the gold standard, the SECIS model and
  non-default parameters, to recommend for new appendicularian assemblies). Test assemblies added on Charles's offer:
  KUM-M3-7f, AOM-5-5f and Bar2_p4 of Plessy et al. 2024 (Zenodo 10.5281/zenodo.10241527, CC0; config/genomes.tsv role
  `test`; fetch.sh now extracts members of .tar.gz archives), new assemblies of the three cryptic species, the main use
  case. Reference for them: gold models lifted with minimap2 (stage2/truth_sets.py; checked locally by lifting
  OKI2018 onto itself: 20 of 20 identical, 3 read as disrupted on the assembly, as in the gold standard; lifts across
  cryptic species find nothing, as intended). Profiles (stage2/build_profiles.py, which replaces
  build_machinery_profiles.py; make machinery-profiles keeps its old options): builds E0 (gold standard), LOO_<g>
  (without one assembly) and final (plus validated stage 1 predictions of O. longicauda, O. albicans and
  O. vanhoeffeni); new profiles KKSU and TXNL1b (no dot: Selenoprofiles reads the family up to the first dot) with
  p2g_refiltering x.label != 'pseudo' instead of the AWSI filter; B. stygius, M. erythrocephalus and F. borealis held
  out of all builds (decision of 2026-10-24, so the machinery is rebuilt without them here). Checked locally without
  Selenoprofiles (--no-build): all 59 gold proteins align their selenocysteine to the selenocysteine column of the
  standard profiles; mafft --keeplength drops 0 to 11 residues per protein. SECIS: the Oikopleura model with GA 20 added
  (stage2/cm_thresholds.py), so the pipeline uses --cut_ga; the pipeline's -T 20 alone would leave inclusion at
  E <= 0.01. Selenoprofiles 3' sequences cut to 2,000 nt (default 6,000), the window length of the 20-bit threshold.
  Runs (stage2/run/): default (stage 1 settings on the new assemblies), eval, eval_miniprot (-C 0.5, for exons at
  GA-AG introns), loo_OKI/OSKA/OdB3, final; benchmark stage2/benchmark.py (tested locally on the stage 1 sites table:
  13, 14 and 13 of 20, 20 and 19 genes in OKI2018, OSKA2016 and OdB3). Pipeline release 0.2.0 has every option
  needed; the samplesheet profiles column is quoted ("eukarya,machinery,KKSU,TXNL1b"; nf-schema reads CSV quotes).
  Not tested locally: selenoprofiles build (installing it here was refused), minimap2 on the real new assemblies.
- Stage 2 first eval run (2026-10-28, Charles): Selenoprofiles failed on B. stygius with the E0 GPX profile
  ("'>' not supported between NoneType and int" in p2g_refiltering): "selenoprofiles build -load" writes back only the
  keys of its parameter list, so go_terms was dropped from GPX.fa.config, and the GPX filter x.go_score() > 0 then
  compares None (no GO terms). build_profiles.py now appends the keys of the standard configuration that the build
  did not write (only go_terms, in GPX, among the families built here; checked on the 15 standard configurations).
  The machinery profiles are not affected (none has go_terms). Rebuild with make stage2-inputs, then resume the runs.
- Stage 2 reference sets, first version (2026-10-28, make stage2-truth on Deigo): lifts of the gold models complete
  for 15 to 16 of 20 genes in KUM, AOM and Bar2 (in-frame indels in 1 to 2); 4, 4 and 2 partial (a coding segment end
  or the codon not aligned: DIO_chr1, GPX_PAR, GPX_chr1a, SELENOO in KUM; DIO_chr1, GPX_chr1b, SELENOO, SELENOW_pair2
  in AOM; DIO_chr1, SELENOS in Bar2) and 1 or 2 disrupted (AOM GPX_PAR, protein identity 0.02 to the gold protein;
  Bar2 GPX_chr1a, identity 1.00, and GPX_chr1b). The benchmark counted only complete lifts. Changed: segment ends in an
  indel are placed from the nearest aligned base (30 nt), the codon at the lifted position is recorded, and the
  benchmark counts every reference whose codon is a TGA (disrupted models are expected in Nanopore assemblies: indels);
  structure compared for complete models only; regions and PAF kept in work/stage2/truth/ for inspection.
- Stage 2 benchmark, first results (2026-10-28, runs default, eval, eval_miniprot, loo_* by Charles). Reference
  codons verified TGA: 20, 20, 19 in OKI2018, OSKA2016, OdB3; 19, 18, 19 in KUM, AOM, Bar2 (lifts). Sec codons
  found (Selenoprofiles label selenocysteine or miniprot): stage 1 13, 14, 13; stage 1 settings on the new assemblies
  12, 12, 11; eval 20, 19, 18 and 19, 17, 18; leave-one-assembly-out (fair test) OKI 20, KUM 19, OSKA 19, AOM 17,
  OdB3 18, Bar2 18. Every reference found by Selenoprofiles has a SECIS hit at the 20-bit threshold. The only gene
  missed in the O. dioica assemblies is SELENOS (OSKA2016, OdB3, AOM, Bar2), with no Selenoprofiles prediction at the
  locus; miniprot with -C 0.5 recovers it in OSKA2016, OdB3 and Bar2 (not AOM) and changes nothing else in the found
  counts. Selenoprofiles found SELENOS in OSKA2016 in stage 1 (standard profile) and in OKI2018 and KUM with LOO_OKI,
  not with E0: the standard SELENOS refiltering requires 5 of 9 residues conserved in vertebrates (the gold proteins
  have 3 or 4) and an AWSI z-score above -2; which clause removes it is tested by make stage2-diagnose (new).
  SELENOT and SEPHS2 are found by Selenoprofiles only, never by miniprot (their selenocysteine codons are 2 and 1 nt
  from a donor site, where miniprot can move the junction to avoid the in-frame stop codon). Structures (eval,
  O. dioica): median coding-base sensitivity and precision 1.00, identical introns in 6 to 10 of 15 to 18 genes.
  Sec sites outside the references in the new assemblies: genes whose lift missed the codon (AOM GPX_chr1b and
  SELENOO, probably KUM DIO_chr1), SELENOO in Bar2 (no OdB3 reference), copies (a second SELENOK 10 kb from the first
  in KUM; DIO copies in Bar2; the OSKA2016 GPX_PAR copy on S12 known from step 4; OdB3 haplotype contigs), and miniprot
  calls 3 nt from the codon (SELENOF, GPX_chr1b); no site outside a selenoprotein gene.
- Filter diagnosis (2026-10-28, make stage2-diagnose, run by Charles): with E0, Selenoprofiles predicts SELENOS in
  OSKA2016, OdB3, KUM, AOM and Bar2 with the label selenocysteine (coverage 0.98-0.99, identity with the profile
  0.29-0.30, AWSI z-score -0.84 to -1.01), and only the motif clause of the standard SELENOS filter (5 of 9 residues
  conserved in vertebrates; the O. dioica proteins have 3 or 4) removes it. The predictions that passed earlier (OSKA2016
  with the standard profile, OKI2018 and KUM with LOO_OKI) were shorter (coverage 0.68-0.72) and aligned so that 5
  residues matched. No other SELENOS candidate of the six assemblies has the label selenocysteine (unaligned, pseudo,
  proline or arginine). Change (stage2/profile_config_overrides.tsv, applied by build_profiles.py): the motif is
  required only for predictions not labelled selenocysteine. The two O. albicans references lost with E0 (DIO.9,
  GPX.24) are fragments of 55 and 49 residues removed by the default length/coverage filter (len > 60 or coverage >
  0.4): with E0 their coverage is 0.20 and 0.25 (0.23 and 0.52 with the standard profile, plus a second exon for DIO.9);
  kept as they are. Next: rebuild the profiles, rerun eval, eval_miniprot and the loo runs, then the benchmark.
- Stage 2 benchmark with the SELENOS change (2026-10-29, rerun by Charles). Every O. dioica reference codon is found
  in eval and in the leave-one-species-out runs: 20/20 in OKI2018 and OSKA2016, 19/19 in OdB3, 19/19, 18/18, 19/19 in
  KUM, AOM, Bar2 (stage 1 settings: 11 to 14). SELENOS is now found by Selenoprofiles in OSKA2016, OdB3, KUM, AOM and
  Bar2 (SECIS 39.3 to 43.8 bits) and by miniprot only in OKI2018. Every codon found by Selenoprofiles has a SECIS hit of
  at least 20 bits. miniprot -C 0.5 changes no count of found genes any more (it adds miniprot support to SELENOS in
  three assemblies, removes it from OKI2018 GPX_chr1a, and moves single extra sites in O. albicans and
  O. vanhoeffeni): not adopted; the pipeline's miniprot options are kept. New with the change: a SELENOS
  selenocysteine prediction in O. longicauda (SCLD01158551.1) and in O. albicans (SCLG01000103.1), without a SECIS hit
  of 20 bits in 2 kb: to validate in step 8 (census). Recommended settings written (stage2/recommended/README.md;
  make stage2-bundle builds the deposit archive). Next: run final (census input), make stage2-bundle.
- Stage 2 final run (2026-10-29, Charles; final profiles, pipeline miniprot options) and deposit archive built
  (work/stage2/appendicularia_selenoprotein_annotation.tar.gz). Same codons found as eval in the six O. dioica
  assemblies (20, 20, 19; 19, 18, 19). In O. longicauda, O. albicans and O. vanhoeffeni, miniprot now finds 22, 16 and 14
  of the stage 1 references (eval: 9, 10, 3), as expected since their own proteins are queries (not a test); the two
  O. albicans fragments stay missed. All sites new in final are miniprot calls 3 to 9 nt from a codon found at the same
  gene or copy (mostly DIO, with the more divergent DIO queries of the other species), plus one B. stygius SELENOW
  fragment and one M. erythrocephalus GPX miniprot site; six eval sites in the other species are no longer called.
  Pipeline improvement to consider (release 0.3.0): merge miniprot sites within a few codons of another site of the
  same strand into one, and ship the stage 2 profiles and SECIS model as an appendicularian option. Stage 2 complete;
  next: manuscript section on stage 2, then step 8 (hit validation, round 2) and step 9 (census, from the final run).
- Manuscript (2026-10-29): stage 2 Results section ("Settings for appendicularian genomes recover the selenoproteome
  of O. dioica in new assemblies", Table STAGE2) and Methods (Stage 2 profiles; Test assemblies and reference genes;
  Stage 2 runs and benchmark), provisional; test assemblies added to the genome Methods, recommended settings to the
  data availability. Inferences flagged in the text: KUM DIO_chr1 as the extra DIO site ("probably"), miniprot moving
  junctions at SELENOT and SEPHS2 ("not tested"). Next: step 8 (hit validation, round 2).
- Step 8 started (2026-10-29): hit validation, round 2 (code; trees to run on Deigo).
  - Taxonomy: new stage2/taxonomy.py, shared by machinery_phylogeny.py, machinery_lba.py, selenoprotein_orthology.py,
    secis_validated_hits.py, machinery_domains.py, machinery_table.py and figures/phylogeny_collapsed.py. Groups from
    the genome code or genus at the start of a name, the Ensembl species prefix of a protein identifier and the
    UniProt species mnemonic (new config/taxon_codes.tsv), then any token that is a known genus (about 210 genera
    added to config/taxon_groups.tsv; new groups Cartilaginous fishes, Hemichordates, Haptophytes, Rhizarians,
    Metamonads, Land plants, Viruses), then descriptive bacterial names. In the first-pass selenoprotein trees, 8
    leaves stay unclassified (no species in the name). Side effects of a rerun at this commit (results of earlier
    commits not recomputed): (1) machinery and selenoprotein trees: sequences formerly unclassified now count in the
    limit of three per group, so the profile subsampling changes; (2) Table MACHINERY: unchanged (checked); (3) step 1
    (make secis-hits): 17 more hits would be "left out" as "nearest relatives non-animal": the selenocysteine DIO of
    OdB3 (5, current-model SECIS 50 to 55 bits), O. longicauda (5), O. albicans (4) and O. vanhoeffeni (2), and
    O. albicans SELENOO.4. In the first-pass DIO tree the selenocysteine DIO clades, curated DIO_PAR and DIO_chr1
    included, have alveolate and amoebozoan DIO as nearest relatives (bootstrap 59 to 66), sequences that the
    earlier taxonomy left unclassified. Step 1 tables kept as committed; DIO added to the AU tests.
  - Trees, round 2 (make selenoprotein-phylogeny; Makefile variables SELENO_PREDICTIONS, GOLD): appendicularian
    sequences from the stage 2 final run; in OKI, OSKA and OdB3, predictions overlapping a gold-standard gene are
    replaced by the gold-standard protein (leaf <genome>_gold_<gene>; removed.tsv records it); KUM, AOM and Bar2 not
    in the trees; human MIEN1 added to the SELENOW references (RELATED_GENES). Best hits (make selenoprotein-rbh) from
    the final-run predictions for the genomes of that run and stage 1 predictions for the others, since prediction
    numbers differ between runs. selenoprotein_orthology.py: gold-standard leaves labelled "gold standard", named by
    gene, with the best hits of the predictions they replaced; the stage 1 naming of curated genes only with
    --annotation (first pass). First pass reproducible with SELENO_PREDICTIONS=work/stage1/pipeline/selenoprofiles GOLD=.
  - Held-out genomes (B. stygius, M. erythrocephalus, F. borealis) kept in the trees: the trees train nothing, and
    these genomes widen the sampling of appendicularians; they stay held out of the profiles and SECIS model.
  - AU tests (make selenoprotein-au, after selenoprotein-pmsf; machinery_lba.py analysis hypotheses;
    stage2/au_hypotheses.tsv): each hypothesis a constraint tree whose clades are expressions of leaf sets (app,
    appclade:, mrca:, group:, re:, joined by + and -), constrained search and AU test against the ML tree, under the
    ML model and LG+C60+F+G with the PMSF site frequencies. Hypotheses: SELENOF, SELENOT with the other tunicates and
    inside the chordates; SELENOO with the selenocysteine SELENOO of the other tunicates and inside the clade of human
    and Ciona selenocysteine SELENOO; SELENOW_PAR outside the smallest clade holding human SELENOW and SELENOV, with
    the appendicularian SELENOW of chromosome 1, with that of chromosome 2; DIO: the animal sequences of the smallest
    clade holding the human DIO and the gold-standard DIO as one clade (non-animal DIO outside).
- Step 8, round 2 trees and AU tests on Deigo (2026-10-30, commit 8e823ae; make selenoprotein-rbh, -phylogeny,
  -pmsf, -au, -orthology). p-AU, ML model / LG+C60+F+G PMSF (alignment sites, sequences):
  - SELENOO (724 sites, 76 sequences): the oikopleurid clade (10, gold OKI and OSKA) is now sister to the animal
    selenocysteine SELENOO (bootstrap 44 ML, 93 PMSF); the first-pass placement next to Phytophthora SelO is gone.
    Inside the human and Ciona selenocysteine SELENOO clade: in the ML tree (0.604 / 0.620); with the selenocysteine
    SELENOO of the other tunicates: 0.086 / 0.238.
  - SELENOF (176 sites, 50 sequences): appendicularian clade (10, nine genomes) sister to plant and green-algal Sep15
    (33 ML, 81 PMSF). With the other tunicates 0.284 / 0.069; inside the chordates 0.524 / 0.105. Neither group is a
    clade of the trees (insect SELENOF among them).
  - SELENOT (269 sites, 63 sequences): appendicularian clade (11, nine genomes) between a group of protist, plant and
    nematode SELT sequences (C. elegans SELT-2 among them) and the other animal SELENOT (90 ML, 86 PMSF). With the
    other tunicates (a clade of the trees) 0.126 / 0.183; inside the chordates 0.088 / 0.0069. The chordates are not
    a clade of the trees even without the appendicularians (nematode and insect sequences inside): the PMSF
    rejection may reject chordate monophyly rather than the appendicularian position; control added.
  - SELENOW (124 sites, 134 sequences): SELENOW_PAR clade (7: gold OKI, OSKA, OdB3; B. stygius, O. albicans,
    O. vanhoeffeni) sister to two profile sequences without species name (Se_user_SelW_1 and _2) and Euplotes
    (89) in the ML tree, sister to the mammalian SELENOV in the PMSF tree (92). Outside the smallest clade of human
    SELENOW and SELENOV: in the ML tree (0.464 / 0.601); with the appendicularian SELENOW of chromosome 1 (chr1a and the
    pair): 0.296 / 0.333; with that of chromosome 2: 0.023 / 0.050. SELENOV arose in placental mammals by duplication
    of SELENOW (Mariotti et al. 2012, doi:10.1371/journal.pone.0033066), so a sister relationship of SELENOW_PAR with
    SELENOV cannot be orthology; the PMSF placement is not required by the data (p-AU 0.601 for the alternative).
    The chromosome 2 clade (19 sequences, nine genomes) has MIEN1 as best human hit (14 of 19) and amphioxus
    sequences as sister (80 / 72). Vertebrate SELENOW is not a clade of these trees (short alignment).
  - DIO (290 sites, 153 sequences): the animal sequences of the smallest clade of the human and gold-standard DIO as
    one clade: 0.431 / 0.347 (not rejected). The ML tree holds DIO_PAR and DIO_chr1 in one appendicularian clade (26,
    support 77); the PMSF tree separates them, DIO_chr1 sister to two alveolate DIO (34).
  - Added (2026-10-30): controls where the reference group is not a clade of the trees (SELENOF other tunicates and
    chordates, SELENOO selenocysteine SELENOO of the other tunicates, SELENOT chordates) and SELENOW_PAR with the
    mammalian SELENOV (the PMSF position) under both models; machinery_lba.py redoes an AU test when its set of trees
    changes (the earlier code kept the old report). To run: make selenoprotein-au (finished constrained searches kept).
- Step 8, AU tests with controls (2026-10-30, commit 952221b; one AU test per family and model, all constrained
  trees together, so the earlier p-values shift slightly). p-AU, ML model / LG+C60+F+G PMSF:
  - SELENOT: with the other tunicates 0.134 / 0.197; inside the chordates 0.118 / 0.009; control (chordates a clade,
    appendicularians outside) 0.042 / 0.012. The data reject chordate monophyly itself (nematode and insect SELENOT
    among the chordate sequences), so the chordate test says nothing about the appendicularian position; the grouping
    with the other tunicates is not rejected.
  - SELENOO: with the selenocysteine SELENOO of the other tunicates 0.034 / 0.045 (rejected at 0.05, narrowly);
    control 0.469 / 0.591; inside the clade of human and Ciona selenocysteine SELENOO 0.544 / 0.594 (the ML tree).
    The oikopleurid clade (gold OKI, OSKA, OdB3) is sister to the animal selenocysteine SELENOO (tunicates, tick,
    amphioxus, vertebrates; 93 PMSF), outside the deuterostome cysteine SELENOO paralog (tunicates, amphioxus,
    turtle, chicken). Other appendicularian cysteine SELENOO copies (F. borealis, B. stygius, M. erythrocephalus)
    group with bacterial SelO (99 to 100): for the census (step 9, contamination and HGT screens).
  - SELENOF: no rejection (0.261 / 0.091 with the other tunicates, 0.412 / 0.125 inside the chordates; controls
    0.317 / 0.610 and 0.598 / 0.461): the position is unresolved.
  - SELENOW: SELENOW_PAR with the chromosome 2 SELENOW 0.031 / 0.054; with the chromosome 1 SELENOW 0.325 / 0.336;
    with the mammalian SELENOV 0.212 / 0.198; outside the human SELENOW-SELENOV clade 0.528 / 0.642. Only the grouping
    with the chromosome 2 (MIEN1-like) clade is rejected, under the ML model only; the PMSF placement next to SELENOV
    is one of several positions the data allow, and cannot be orthology (SELENOV is placental; Mariotti et al. 2012).
  - DIO: 0.431 / 0.347 (unchanged; one hypothesis).
  Conclusion of step 8: no appendicularian selenoprotein clade is placed outside the animals by a test; the
  gold-standard genes keep their family names. SELENOW_PAR stays named after SELENOW (its SELENOV-like position is not
  supported); the chromosome 2 SELENOW are MIEN1-like by best hit and placement but selenoproteins (gold standard),
  as in step 1. Next: step 9 (census from the final run, contamination screen, strict HGT screen), with the
  appendicularian cysteine SelO near bacteria and the SELENOS of O. longicauda and O. albicans to check there.
- Manuscript, step 8 (2026-10-30): Results section "Orthology of the selenoprotein genes" with Table SELENO-ORTHOLOGY
  (stage2/selenoprotein_table.py, make selenoprotein-table; one row per gold-standard gene) and Methods; step 1
  Methods note: with the extended taxonomy, 24 of the 111 selenocysteine predictions would have non-animal nearest
  relatives only (alveolate and amoebozoan DIO, among others), none bacterial, archaeal or fungal only (checked on
  the first-pass trees). SELENOW_chr2: best hit MIEN1 (no hit to SELENOW or SELENOV) and in a group (77 / 91) with
  selenocysteine SelW-like proteins of zebrafish, tilapia, Ciona, amphioxus and Artemia, and human MIEN1 (cysteine at
  the selenocysteine column); described as a SELENOW-like lineage, names kept. AU tests added for SELENOK and SELENOS
  (non-animal sister groups in round 2), with controls for SELENOK; to run (make selenoprotein-au).
- Step 9 started (2026-10-30): census code (Tiara part dropped later the same day, see below). stage1/collect.py now also writes the Selenoprofiles identifier of each
  codon and a predictions table (all labels, coding exons); make stage2-collect applies it to the final run and
  reruns the machinery census on the final-run predictions. Contamination screen: Tiara (scripts/run_tiara.sh, own
  image: Tiara needs Python < 3.10, outside environment.yml) and composition checks (stage2/census_contigs.py, make
  census-contigs: GC and tetranucleotide distance of a 2-10 kb window around the predictions against random windows
  of the same genome, as percentiles and robust z-scores; tested on a synthetic genome with a bacterial contig: robust
  z 97.5 against -1.1 for a host contig, percentile 98.7 diluted by the contaminant's own windows). stage2/census.py
  (make census): loci (codons within 30 nt or in one prediction merged), decisions fixed in advance (left out: Tiara
  non-eukaryotic, best human hit in another family, or composition outlier z > 5 with non-animal nearest relatives;
  gold-standard genes always counted), flags, counts against the other tunicates (stage 1, no screen), and the strict
  HGT screen (round 2 leaves with bacterial, archaeal or fungal nearest relatives; status contamination, not testable
  without intron, placement not stable under PMSF, or candidate). Tested on synthetic loci built from the round 2
  trees: the DIO cysteine clade of B. stygius, F. borealis and M. erythrocephalus (nearest relatives two bacterial
  DIO, 88%) and the cysteine SELENOO near bacteria appear in the screen as expected.
- Step 9, HGT screen revised (2026-10-30, question of Charles: what if a transfer is 500 My old, like that of
  cellulose synthase?). An old transfer would match its genome in composition (amelioration) and lie on a contig that
  Tiara calls eukaryotic, so the contig checks do not exclude it; but the first screen looked only at the nearest
  relatives outside the appendicularians, and missed a transfer older than the split from the other tunicates (other
  tunicates in between). The screen now also reads the nearest relatives outside the tunicates (two levels), and
  records the genomes of the appendicularian clade; an intronless candidate needs presence in other assemblies or
  species instead of being "not testable". First case found by the tunicate level, in the round 2 DIO tree: a clade
  of 41 appendicularian and 4 other tunicate (Pegea, Styela, Salpa, Thalia) sequences next to two flavobacterial
  sequences of the DIO profile (88 ML, 99 PMSF). These are thioredoxins (WCGPC/WCPPC; best human hit TXN, 63 to 100
  bits, and TXNL1, NHLRC2, DNAJC10): the DIO tree holds no animal thioredoxin, so their nearest relatives are the
  bacterial thioredoxin-like sequences of the profile. Status added: "other family" (best human hit outside the
  family of the tree) and "family not confirmed" (no human hit although the family has human members); a test of an
  old transfer needs a tree of the right family with broad eukaryotic and bacterial sampling.
- Decision (2026-10-30, Charles): Tiara is not run. Only old transfers are expected. Hypothesis of Charles: the
  ancestral tunicate had a symbiont, since lost; losses in parallel after the split of appendicularians and other
  tunicates would leave room for appendicularian-specific transfers (improbable, but transfer without symbiont is
  judged less probable). Old transfers match their genome in composition and lie on host contigs, so contig
  classification adds little; the evidence is presence at the same place of the tree in several species and
  assemblies, and introns shared between species (as for the tunicate GH6-containing genes: Li et al. 2020, Genes,
  doi:10.3390/genes11080937, shared splice sites in all tunicates examined; GH6-1 and the GH6 part of CesA form two
  independent monophyletic groups). scripts/run_tiara.sh removed; census_contigs.py still reads a Tiara output if
  one is placed in work/stage2/census/tiara/. The HGT screen lists the foreign sequences next to each candidate, to
  compare donor lineages: one lineage across candidates fits an ancestral symbiont, several lineages fit uptake of
  bacteria with food (appendicularians filter bacteria from the water). To consider for candidates: intron positions
  compared between species. The screen covers only the genes of the selenoprotein families.
- AU tests of SELENOK and SELENOS (2026-10-30, commit 1d630fb). SELENOS: with the other tunicates 0.19 / 0.51,
  inside the chordates 0.29 / 0.51 (not rejected). SELENOK: every constraint rejected (p ≤ 0.0016), controls
  included, at nearly the same cost (deltaL 70.6 against 70.3 for the other tunicates, 63.1 against 62.9 for the
  chordates, ML model): the SELENOK of the other tunicates and of the chordates are not clades in these data (181
  sites); uninformative for the appendicularians. One PMSF control p-AU printed as 7.6e-61 (deltaL 69.0, as its
  hypothesis): a numerical artefact of the AU fit, read as p < 0.001. Written into the manuscript.
- Final run collected (2026-10-30, make stage2-collect, commit 1d630fb): census preview without contig checks
  (make census-contigs not run yet): 247 loci in the twelve assemblies, 59 of them gold-standard genes, one left out
  (F. borealis DIO.14, best human hit GPX1). HGT screen: no candidate with introns; 27 intronless candidates whose
  family is confirmed by best hits, all in the three non-Oikopleura genomes: cysteine SELENOO of B. stygius,
  F. borealis and M. erythrocephalus next to Stenotrophomonas and Xylella SelO (100 / 99), and further F. borealis
  SELENOO next to a Pseudomonas protein; cysteine GPX of F. borealis (and two of B. stygius) next to Salmonella,
  Cellvibrio or Delftia GPX; one F. borealis MsrA. B. stygius may be contaminated by M. erythrocephalus (Pichon et al.
  2019), so their shared clade may come from one source. Contig checks needed before reading further.
- Introns and NMD (2026-10-30, question of Charles). Single-exon Selenoprofiles predictions in the final run:
  OKI 12/51, OSKA 18/56, OdB3 18/66, Olon 60/133, Oalb 32/72, Ovan 35/67, Bsty 39/55, Mery 30/55, Fbor 211/239. In
  F. borealis most splice junctions are non-canonical (Charles, 2026-10-24) and the prediction tools do not model
  them; B. stygius and M. erythrocephalus are fragmented: "no intron" in the HGT screen says little for these three
  genomes (to check with transcripts or a spliced alignment that allows non-canonical junctions). Gold standard: 5 of
  59 models intronless (SELENOW_PAR in the three assemblies, SELENOW_pair2 in OSKA2016 and OdB3); the selenocysteine
  codon lies more than 55 nt upstream of the last exon junction in 32 models (the mammalian rule for
  junction-dependent NMD; selenoprotein mRNAs are NMD substrates when the UGA is read as a stop: Moriarty et al. 1998,
  doi:10.1128/MCB.18.5.2932; Sun et al. 2001, doi:10.1091/mbc.12.4.1009) and in the last exon in 18.
- Census contig checks, first run (2026-10-30, commit d66521f). Assemblies: N50 1.3 kb (B. stygius), 1.9 kb
  (M. erythrocephalus), 2.4 kb (O. longicauda), 2.8 kb (F. borealis), 0.2 to 16 Mb for the others. Composition
  outliers (robust z > 5) among the contigs with predictions: OdB3 1, O. longicauda 1, M. erythrocephalus 4,
  F. borealis 9 (109 of its 222 contigs shorter than 2 kb, untested). Census decisions unchanged. HGT screen: 7 entries
  now "contamination likely (composition outlier)" (F. borealis cysteine DIO-like and GPX, one SELENOO), 25 intronless
  candidates left. Caveat found: in fragmented assemblies the 10-kb reference windows come only from the few long
  contigs; the B. stygius SELENOO candidates lie on 49- and 137-kb contigs and B. stygius GPX.1 on a 178-kb contig,
  long for an assembly of N50 1.3 kb, which is itself a reason to suspect a high-coverage contaminant, and they are
  compared with windows of such contigs (z -0.5, -0.2, -0.7). Added: comparisons on the 2-kb window against 2-kb
  reference windows; open reading frame density (fraction of bases in open frames of at least 100 codons, six frames)
  of the window and of the whole contig (tests: random 40% GC 0.03, random 65% GC 0.84, coding-dense 0.97; hence the
  rule only at GC up to 50%); per genome, the fraction of the assembly in contigs of each window length and the
  reference medians. New rule fixed before the rerun: ORF-dense contig (at least 0.6 of the contig in long open
  frames, GC up to 50%, contig of 5 kb or more) treated as a composition outlier. To rerun: make census-contigs.
- Census contig checks, second run (2026-10-30, commit 26ec103). Median open reading frame density of the 2-kb
  windows: O. dioica 0.41 to 0.44, O. longicauda 0.34, O. albicans 0.16, O. vanhoeffeni 0.00, M. erythrocephalus 0.15,
  B. stygius 0.21 (but 0.76 for 5-kb and 0.89 for 10-kb windows: its contigs of 10 kb or more, 2.9% of the assembly,
  are mostly ORF-dense, bacterial-like), F. borealis 0.80 at every window length. The F. borealis selenoprotein genes
  with SECIS hits lie on contigs of the same density (GPX.4, GPX.5, GPX.12, SELENOW.4: 0.70 to 0.77): either a
  compact genome with long open frames (most of its predictions are single-exon) or an assembly largely made of
  non-animal DNA; these checks cannot tell (Tiara or a database search could). Rule revised after seeing these data:
  an ORF-dense contig must also exceed the median of its genome by 0.3 (as first written, the rule left out
  F. borealis SELENOF.2, in DNA typical of its assembly). Census: 247 loci, 59 gold-standard genes, one left out
  (F. borealis DIO.14, best human hit GPX1). HGT screen: the B. stygius (6) and M. erythrocephalus (3) entries,
  cysteine SELENOO, DIO-like and GPX, lie on ORF-dense contigs of 26 to 324 kb (0.88 to 0.93) and are "contamination
  likely": the SELENOO and DIO-like clades shared by B. stygius, M. erythrocephalus and F. borealis are best explained
  by bacterial contamination of one lineage in the three samples. Left: 21 intronless F. borealis entries, 16 of them
  on contigs shorter than 3 kb, not testable with these data. No credible candidate of horizontal transfer in the
  selenoprotein families.
- Step 9, census (2026-10-30). census.py now assigns each locus to the O. dioica gene of its clade (gold-standard
  gene; lifted model in KUM, AOM, Bar2; else the gold-standard genes of the smallest appendicularian clade holding it
  in the round 2 ML tree) and adds the loci of the step 4 screens that the final run did not predict (22, without
  gene model; among them SELENOK of O. longicauda, O. albicans, O. vanhoeffeni and M. erythrocephalus, which the final
  run missed: Olon only as a 46-aa "unaligned" fragment, the others not at all; TXNRD of O. vanhoeffeni and
  M. erythrocephalus). Table CENSUS (results/tables/stage2_census_genes.md): O. dioica 20 genes in OKI2018 and
  OSKA2016, 19 in OdB3 (SELENOO exon missing), the same in the new assemblies; loci of 11, 17, 15 of the 20 genes in
  O. longicauda, O. albicans, O. vanhoeffeni, 8 and 7 in B. stygius and M. erythrocephalus, 5 in F. borealis (with
  SELENOL x2, SELENOU x1, 7 GPX outside the O. dioica clades; no KKSU). TXNL1.b only in O. dioica; SELENOW pair only
  in O. dioica (plus one O. albicans screen locus assigned by its partner). Not found in any appendicularian but in
  other tunicates: SELENON (8 of 10), FrnE (8), SELENOH (6), MSRB (6), SELENOI (C. intestinalis); SELENOL in 9 other
  tunicates and only F. borealis. No absence tests for these families (to consider: tblastn of the ascidian proteins,
  as for the machinery). Results and Methods written into the manuscript (section "The selenoproteomes of
  appendicularians", Tables CENSUS and CENSUS-FAMILIES; Methods "Census of the selenoprotein genes"). Next: absence
  tests for SELENON, SELENOH, MSRB, FrnE, SELENOL and SELENOI if wanted; step 10 (SECIS and operons).
- Step 9 follow-up 1, absence tests (2026-10-30, Charles: "one by one"). Families found in other tunicates but in no
  appendicularian: SELENON, SELENOH, MSRB, FrnE, SELENOL, SELENOI. In the final run: MSRB cysteine forms in every
  appendicularian, best human hit MSRB3 (or MSRB2), never MSRB1; SELENON only a cysteine form in F. borealis; FrnE
  fragments ("unaligned") in O. dioica, B. stygius, M. erythrocephalus and cysteine forms in F. borealis; SELENOL only
  F. borealis; SELENOH and SELENOI nothing. machinery_absence.py generalised (--labels, --no-profile,
  --outgroup-human); make selenoprotein-absence: selenocysteine forms of stage 1 as queries and references,
  outgroups CALM1 (SELENON EF-hand), SELENOW/SELENOV/SELENOT (SELENOH), MSRB2/MSRB3 (MSRB, references MSRB1-type only),
  GSTK1 (FrnE), TXN (SELENOL), CEPT1/CHPT1 (SELENOI). To run on Deigo.
- Absence tests, first run (2026-10-30, commit 4223f79). Positive controls: margins of 45 to 593 bits. SELENOH: no
  candidate in any genome or transcriptome. SELENON: none in the genomes; one 36-aa ORF of an OSKA2016 transcript
  (comp56540_c0_seq1), 77.8 bits to a human SELENON isoform of the profile, no outgroup hit: to inspect (a near
  identical human-like fragment would point to contamination of the transcriptome). MSRB: every candidate closer to
  human MSRB2 or MSRB3 than to the MSRB1-type selenocysteine forms (margins -20 to -56; F. borealis 0): no MSRB1.
  FrnE and SELENOL: only F. borealis (supported, 177 and 46 bits). SELENOI: supported candidates in every genome and
  in the three O. dioica transcriptomes (O. dioica transcripts: ORFs of about 390 aa, 285 bits to zebrafish SELENOI,
  175 bits closer than to CEPT1): an EPT1/SELENOI ortholog is present; whether it carries a selenocysteine is not
  known, since the final run predicted no SELENOI. machinery_absence.py now reports, for each candidate, the residues
  aligned to the selenocysteine of its best reference (TGA read as U; transcript ORFs read through TGA) and the
  protein. To rerun: make selenoprotein-absence.
- Absence tests, rerun (2026-10-30, commit d9b0868): SELENOI is a candidate selenoprotein of O. dioica. Transcripts of
  the three species (OKI2018 Trinity TRINITY_DN1283_c0_g1_i2 and _i5, OSKA2016 comp21960_c0_seq1, OdB3 EST
  KT0AAA233YG01RM2), read through TGA, encode an EPT1/SELENOI protein (394 aa in OKI2018, 285 bits to zebrafish
  SELENOI, 175 bits closer than to CEPT1) with an in-frame TGA aligned to the selenocysteine of zebrafish SELENOI, in
  a conserved C-terminal motif ...VFKIIPTKTDU-QI; other isoforms of the same genes end ...VFKIIPTS(I/M)SS(L/I)
  without it, apparently an alternative last exon. In OKI2018 a second in-frame TGA follows 3 codons later (OSKA2016
  H, OdB3 P). Missed by stage 1, the final run (no SELENOI prediction at all) and the step 4 screens (no SECIS hit
  above threshold in these transcripts). The genomic candidates of the other species stop before the selenocysteine
  region (short exon). SELENON: the OSKA2016 transcript fragment (36 aa, GEPLQFVFEEIKWQQELSWEEAARR..., 77.8 bits to a
  human SELENON isoform) is near identical to human SELENON: contamination of that transcriptome, not counted.
  SELENOI added to step 5's candidate genes (curation/candidates.faa; codons located by tblastn in the three
  assemblies, then the rule: TGA and SECIS of at least 20 bits in the 2 kb after the stop codon in the three
  assemblies; other genomes and proteomes as for TXNL1.b and KKSU). To run: make candidate-genes.
- make candidate-genes with SELENOI (2026-10-30, commit da275c9): only the proteome table has SELENOI. The genome
  searches were taken from the cache of the previous run, whose test (the cached candidates a subset of the current
  ones) let new candidates through unsearched: fixed (the candidate proteins searched are recorded beside the rows,
  and any change triggers a new search). Proteomes: the human, mouse and zebrafish SELENOI have their selenocysteine
  aligned to the candidate U (positions 387, 388, 387; 10 to 13 residues after it); the ascidian and amphioxus EPT1
  RefSeq proteins (374 to 392 aa) do not reach that position, except B. schlosseri (H, by extension). To rerun.
- make candidate-genes rerun (2026-10-30, commit d470695): SELENOI codon located by tblastn in OKI2018
  (NC_137292.1:6700889, -) and OdB3 (FN653023.1:889579, -), TGA in both, contexts KIIPTKTDUQI*KQK (OKI2018: the
  second TGA ends the frame, 2 codons after) and KIIPTKTDUQIPEKVQFDNHRSSKFE*IAE (OdB3: TAA, 17 codons after); no
  SECIS hit of either model (cmsearch, 0 bits) in the 2 kb after the stop codon in either assembly; OSKA2016 (S24): U
  position not aligned by tblastn (short last exon). Rule not met: SELENOI stays a candidate. Other genomes: U position
  not aligned (Olon, Oalb, Ovan, Bsty, Mery, ascidians, amphioxus), Fbor ATT (I, extended), Csav AAA (K, extended).
  To complete the evidence: stage2/candidate_transcripts.py places the SELENOI transcripts of the three assemblies on
  their assembly (minimap2 alignments of make gold-standard, junction placement of gold_standard.py), maps the codon
  aligned to the candidate U, and searches SECIS in the transcripts after the stop codon; candidate_genes.py
  (--transcripts) uses these positions for the rule (OSKA2016 included). Tested on a synthetic gene. To run: make
  candidate-genes (minutes; genome searches cached).
- make candidate-genes with candidate_transcripts.py (2026-10-30, commit 27e8a7c). SELENOI transcripts placed on the
  assemblies: OKI2018 TRINITY_DN1283_c0_g1_i2 and _i5 (NC_137292.1:6700889, -, exon 7 of 7; same position as
  tblastn), OSKA2016 comp21960_c0_seq1 (S24:93093, -, exon 7 of 7; KIIPTKTDUQIHQKSSIR*, TAA 9 codons after), OdB3
  EST KT0AAA233YG01RM2 (FN653023.1:889579, -, exon 5 of 5 of the EST; same as tblastn). TGA in the assembly in the
  three. No SECIS hit of either model in the transcripts after the stop codon (402, 255 and 544 nt) nor in the 2 kb
  of genome after it. Isoforms with the alternative end (VFKIIPTS(I/M)SS(L/I)): OKI2018 _i3, OSKA2016
  comp21960_c0_seq2, OdB3 three ESTs. Verdict: candidate (rule not met), not counted; the codons after the TGA differ
  between the three cryptic species while those before it are conserved. Other transcripts of the absence table
  (TRINITY_DN530, comp19304: CEPT1/CHPT1-like, TGT at the aligned position) are paralogs.
- Step 9 follow-up 1 closed (2026-10-30): absence tests and SELENOI written into the census section of the
  manuscript (Results: two paragraphs after the comparison with other tunicates; Methods: "Absence tests" extended to
  the selenoprotein families, "Candidate genes" with candidate_transcripts.py). SELENON, SELENOH, MSRB1 absent from
  appendicularians, FrnE and SELENOL from the Oikopleuridae (within the limits of tblastn on these assemblies).
  Next: follow-up 2 (screen-found SELENOK proteins into the SELENOK profile before the deposit), then step 10.
- Step 9 follow-up 2 started (2026-10-30): SELENOK of the other Oikopleura species into the SELENOK profile before the
  deposit. The screens found SELENOK in O. longicauda (SECIS scan, codon not placed; final run: a 46-aa "unaligned"
  fragment), O. albicans and O. vanhoeffeni (comparative screen, C-terminal open stretches of 34 and 21 aa ending
  GGUGR) and M. erythrocephalus (held out; 25 aa); the final run predicted none of them. The screen frames are
  C-terminal fragments, not usable as profile members. First, why the final SELENOK profile (24 sequences; filters:
  blast E < 0.01, p2g_refiltering identity with profile >= 0.12) misses them: make selenok-diagnose
  (diagnose_refiltering.py with SELENOK and the final build in OKI2018, the three Oikopleura species and the three
  held-out appendicularians; candidate predictions before filtering, now with their proteins). If the genes are
  found and filtered, their predicted proteins (Olon, Oalb, Ovan only; Mery stays held out) can become profile
  members; if not found, protein-to-genome alignment at the screen loci.
- make selenok-diagnose (2026-10-31, commit 5c7698f): the final SELENOK profile gives no candidate prediction at all,
  before filtering, in O. albicans, O. vanhoeffeni and M. erythrocephalus; in O. longicauda an "unaligned" fragment
  (SCLD01109698.1:5846-5983, -, coverage 0.55; the SECIS scan locus) and a "serine" fragment (SCLD01095635.1);
  in B. stygius an "unaligned" fragment, in F. borealis four "glycine" fragments; OKI2018 control found (coverage
  0.99). The genes are missed at the search, not by the filters (short exons, as for TXNL1 and KKSU). The profile
  members must therefore be built from the genomes: make selenok-regions (stage2/extract_regions.py,
  curation/selenok_loci.tsv) extracts 9 kb around each screen locus; gene models then built here by protein-to-genome
  alignment of the O. dioica SELENOK. results/sequences/ (small FASTA and TSV files) no longer ignored by git.
- SELENOI SECIS, wider search (2026-10-31, Charles: "search harder for a SECIS in / near SELENOI ... by just
  screening for stem-loops first"). make selenoi-regions (commit 91fc51e) gave the sequences; analysis
  stage2/selenoi_secis_search.py (make selenoi-secis; run here first, numbers below). Findings:
  - The TGA is also in O. albicans and O. vanhoeffeni (tblastn of the 45 residues before the U on their regions):
    ...KTN[TGA]QI in both; the five species read TGA-Q-I, then diverge (stop after 2 codons in OKI2018 and
    O. vanhoeffeni, 9 in OSKA2016, 17 in OdB3, 24 in O. albicans). Not found in O. longicauda (contig of 1.2 kb),
    B. stygius, M. erythrocephalus, F. borealis.
  - The TGA lies 9 nt inside a segment of 103-108 nt that the other transcripts splice out (OKI2018 _i3, OSKA2016
    seq2, three OdB3 ESTs; they end ...IIPT|SISSL*); the TGA transcripts retain it. The only placement of that
    junction consistent with the genome has AG...AG ends, not one of the O. dioica intron classes: to check.
  - A gene on the same strand, conserved in the three O. dioica species (MKLAVTLLSCSVVLA...), starts 306-348 nt
    after the TGA; the transcripts end 250-324 nt after it.
  - Covariance models without filters (cmsearch --max, 3 kb after the TGA): best Oikopleura model 9.5 to 13.5 bits,
    appendicularia model 10.4 to 15.4, RF00031 5.8 to 9.5; dinucleotide shuffles of the same sequences: medians 10.1,
    10.3, 7.8, 95th percentiles 15.4, 15.1, 11.0. Positive controls: the vertebrate SELENOI SECIS 1.1-1.3 kb after the
    codon (Oikopleura model 18.3-26.3, RF00031 17.4-22.5). C. intestinalis (stage 1 SELENOI site, miniprot,
    ...ENSD[TGA]HAA...): no SECIS either (11.8, 9.3 bits).
  - Stem-loop screen (first 600 nt): no stem-loop with the SECIS features in OKI2018, OSKA2016, O. albicans,
    O. vanhoeffeni; OdB3 437-594 (inside the downstream gene; z -0.93 to 1.32, at most 12.2 bits). Calibration on the
    gold-standard 3' sequences: 49 of 59 known SECIS in a stem-loop, 24 with the features; background 0.26 per kb.
  - Core patterns (no folding, no model): 0 to 3 per 600 nt (background 1.73 per kb in the shuffles; 34 of 59 known
    SECIS recovered, 1.36 per kb in the gold-standard background).
  - Five-species alignment of the 450 nt after the TGA (identity 38-90%): RNAalifold consensus of 8 base pairs,
    -1.72 kcal/mol, no conserved stem-loop.
  Reading (inference): no SECIS; the TGA most likely ends a truncated protein (...KTD) of the intron-retaining
  transcripts; SELENOI stays "candidate, not counted". To run on Deigo for the tables: make selenoi-secis.
- .gitignore (2026-10-31, Charles: "pay attention again to .gitignore for FASTA files"): the exceptions for
  results/sequences/ came before the global *.fa, *.sto patterns, which therefore won; moved after them (FASTA,
  protein FASTA, Stockholm and TSV files of results/sequences are tracked).
- SELENOK models (2026-10-31). make selenok-regions (commit d5ea80b). exonerate protein2genome found two exons in
  O. albicans and O. longicauda and nothing in O. vanhoeffeni and M. erythrocephalus; tblastn of the O. dioica
  protein finds only short pieces (the C-terminal half is Gly-rich). region_gene_models.py now uses genewise (best of
  the three gold-standard queries), joins the last exon to the exon of the screen codon (donor GT/GC/GA near the end
  of the model, acceptor AG within 400 nt upstream of the TGA, codon in frame, best global alignment score) and takes
  the start codon at most 10 codons upstream (or among the first 2 codons). Run here on the regions (O. longicauda
  offset corrected to the codon): O. longicauda 97 aa, three exons (introns of 94 and 48 nt), TGA at SCLD01109698.1
  5753-5755 (-), the codon given by the comparative screen as partner of M. erythrocephalus cand788, 76 nt before
  the SECIS of the scan (39.0 bits); O. albicans 98 aa, three exons, TGA at the screen codon; O. vanhoeffeni 66 aa,
  N-terminus incomplete (first exon not found), TGA at the screen codon; M. erythrocephalus 21-aa C-terminal fragment
  (held out). Identity to the O. dioica protein 51-55%. curation/selenok_loci.tsv: O. longicauda position now the
  codon (5753). New target stage2-release-profiles: the final profiles plus these models (M. erythrocephalus left
  out by the held-out rule of build_profiles.py), as a separate build ("release") for the deposit, so that the final
  run stays tied to the final profiles; test: make selenok-diagnose SELENOK_DIAG_BUILDS="final release" (held-out
  M. erythrocephalus, B. stygius and F. borealis are the real test).
- SELENOK release profiles (2026-10-31, commit 7e292d9). Models (Deigo): O. longicauda 97 aa (3 exons, complete),
  O. albicans 98 aa (3 exons, complete), O. vanhoeffeni 66 aa (2 exons, first exon missing), all with the TGA at the
  screen codon, 51-55% identical to the O. dioica protein; M. erythrocephalus 21-aa fragment, held out (not a member).
  Release SELENOK profile: 27 sequences (24 + 3). selenok-diagnose, final vs release: OKI2018 found by both; the three
  species of the models now predicted with their selenocysteine (coverage 0.99, 0.99, 0.68; expected, since they are
  members); held out: M. erythrocephalus none -> an "unaligned" N-terminal fragment on SCLF01722808.1
  (VDFVTGLFWGAVNFVTIFFSSMLNLDPP...; the screen codon is on SCLF01573378.1); B. stygius "unaligned" coverage 0.67 ->
  0.88 (N-terminus now complete), no selenocysteine exon; F. borealis only Gly/Pro-rich fragments with both. Side
  effect: more Gly/Pro-rich fragments without selenocysteine pass the SELENOK filters (OKI2018 2, O. longicauda 2,
  O. albicans 1), as the final profile already does in OSKA2016, AOM-5-5f and F. borealis. Decision: the deposit
  bundle carries the release profiles under the recommended name and the final profiles of the census
  (appendicularia_selenoprofiles_final.tar.gz). Written into the manuscript (census section, Methods "Stage 2
  profiles"; GeneWise cited, Birney et al. 2004, doi:10.1101/gr.1865504, via PubMed). Step 9 follow-up 2 closed.
  Next: step 10 (SECIS and operons).
- Step 10, SECIS and operons (2026-10-31, Charles: "let's do it"). Questions: which selenoprotein genes of the three
  O. dioica assemblies lie in operons and at which position; where their SECIS lies relative to the next gene (the
  first gene of an operon has a short 3' UTR); whether their mRNAs carry the spliced leader (SL); whether their
  neighbours are conserved between the cryptic species (Plessy 2024: operons mostly not conserved between Okinawa,
  Osaka and Barcelona); background of all genes. Operon definition: consecutive loci on the same strand with coding
  spans at most 500 nt apart (Plessy 2024: 500 bp intergenic; coding spans here because the GenBank OdB3 annotation
  has no UTRs). SL 3' end read in the transcripts of this study (TRINITY_DN1283_c0_g1_i2 and others start with
  CTCCTCCCATTTTTGAGTCCGATTTCGATTGTCTAACAG; OSKA2016 transcripts with TTTCGATTGTCTAACAG); the Ganot 2004 full text was
  not available through PubMed for the published sequence. stage2/operons.py (make operons), tested on a synthetic
  contig. Side finding: SELENOI mRNAs are SL trans-spliced (DN1283 isoforms start with the SL).
- SL sequence (2026-10-31, Charles, with the Ganot et al. 2004 PDF): ACTCATCCCATTTTTGAGTCCGATTTCGATTGTCTAACAG (40 nt,
  Fig. 1A). Ganot 2004 also: trans-splice acceptor consensus TTT(C/T/A)AG, most trans-spliced exons start with A
  (133 of 145); intercistronic regions of 23-30 nt (polyadenylation site of the upstream gene to the trans-splice site
  of the next) in the cyclin D3 cluster; only the last gene (cycD) had an AAUAAA signal; candidate operons defined by
  less than 300 bp between coding sequences (Fig. 7). operons.py: the published SL cited (search with its 3' 13 nt,
  since transcript assemblies can lack its 5' part), acceptor AG checked at each SL transcript start; gaps reported so
  that the 300-bp definition can be applied.
- make operons, first run (2026-10-31, commit cc53285). Background (coding spans <= 500 nt, same strand): 47% of the
  loci in operons in OKI2018 (15,364 loci), 40% in OSKA2016 (18,743), 49% in OdB3 (17,208). Selenoprotein genes in
  operons: 10 of 20, 13 of 20 (Fisher p = 0.037), 13 of 19; first position: 1, 5, 3 (not depleted; p = 0.10, 1,
  0.56). Every SECIS lies before the coding sequence of the next gene on the same strand, with at least 63 nt to
  spare (SELENOW_PAR, internal in operons of 7 and 11 genes in OKI2018 and OdB3: SECIS at 45-118 nt after the stop,
  next gene 182-183 nt after it), except SELENOW_pair2: in the three assemblies the only SECIS within 2 kb after its
  stop codon is that of SELENOW_pair1 (same score), beyond the start of SELENOW_pair1 (OKI2018: pair1 starts 197 nt
  after the pair2 stop, the SECIS at 584-648); transcripts run from pair2 into pair1 (6 in OKI2018, 1 in OSKA2016,
  2 in OdB3), so pair2 may depend on the pair1 SECIS in a bicistronic mRNA (to examine). Conserved neighbours: for
  example SELENOS internal in all three with neighbours on both sides homologous between the assemblies; SELENOK,
  SELENOW_chr2 and DIO_PAR last in all three with a conserved upstream neighbour. Bug: the acceptor check read the
  genome before the start of the alignment, but minimap2 aligns the AG at the end of the SL to the genomic AG, so no
  acceptor was found; fixed (the trans-splice site is placed from the end of the SL in the transcript). Added: the SL
  trans-splice sites of the next gene between the two coding sequences, and the space between the SECIS and the
  first of them. Tested on a synthetic contig (alignment starting inside the SL, antisense transcript). To rerun.
- make operons, second run (2026-10-31, commit 805135c): trans-splice sites now AG in nearly every case (for example
  OKI2018 SEPHS2 10 of 10, GPX_PAR 16 of 16); 5' UTRs of the SL-spliced selenoprotein mRNAs mostly 0-34 nt. SL found
  for 12 of 20 genes in OKI2018 and 15 of 20 in OSKA2016, 4 in OdB3 (ESTs with masked 5' ends). Space between the
  SECIS and the first trans-splice site of the next gene: 66 nt (OKI2018 SELENOW_PAR), 75 (OSKA2016 SEPHS2), 93
  (OSKA2016 TXNL1.b), 114-149 (SELENOF, SELENOS), more for the others. Correction of the first-run note on
  SELENOW_pair2: it has its own SECIS 23-92 nt after its stop codon in the three assemblies (42.7, 35.8, 36.6 bits),
  before SELENOW_pair1; operons.py took the best hit of the 2-kb window, which was the pair1 SECIS. Fixed: the SECIS
  of a gene is the best hit that ends before the next gene on the same strand (all hits listed). No bicistronic
  dependence. To rerun.
- make operons, third run (2026-10-31, commit dbefde6), written into the manuscript (Results "Selenoprotein genes in
  operons", Methods "Operons"). Numbers: selenoprotein genes in operons 10/20, 13/20, 13/19 (background 47%, 40%, 49%;
  P = 0.83, 0.037, 0.11); first position 1/10, 5/13, 3/13 (background 37%, 40%, 35% of operonic loci; P = 0.10, 1,
  0.56); SL transcripts for 12, 15 and 4 genes, 74 of 75 trans-splice sites AG, 55 of 75 within 40 nt of the start
  codon (median 2-6 nt); SECIS ends 82-341 nt after the stop codon; next coding sequence at least 63 nt after the
  SECIS (OdB3 SELENOW_PAR), next trans-splice site 66-516 nt after it (OKI2018 SELENOW_PAR 66); SELENOW_pair2 has its
  own SECIS (23-92 nt after the stop); in OdB3 an annotated CDS (CBY21828.1, neighbour not conserved) starts 47 nt
  after the pair2 stop, inside its SECIS. 19 genes in the three assemblies: 13 same operon status, 9 same position.
  Step 10 done. Next: step 11 (Discussion, Zenodo deposits), diary update.
- 3' UTRs of operon genes (2026-10-31, Charles: "a comparison of UTR length in operon genes and in operon
  selenoprotein genes"). operons.py now also compares, for the genes followed by another gene of their operon (first
  or internal), the distance from the stop codon to the next coding sequence and to the first trans-splice site of the
  next gene (SL transcripts of the whole transcript set), and the 3' UTR length for all genes (annotated exons beyond
  the coding sequence; gold-standard models for the selenoprotein genes, from transcripts; the OdB3 annotation has no
  UTRs): selenoprotein genes against the others, Mann-Whitney U test; distributions of the others in 25-nt bins
  (results/tables/stage2_operons_utr.tsv, stage2_operons_utr_bins.tsv) for a panel of Figure OPERONS. Tested on the
  synthetic contig. To rerun: make operons.
