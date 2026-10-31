# Stage 2

## Appendicularian SECIS model

`build_secis_model.py` bootstraps a SECIS covariance model from the stage 1 outputs: the 3′ sequences that
Selenoprofiles wrote after its predictions labelled selenocysteine, in seven oikopleurid genomes (OKI2018, OdB3,
*O. longicauda*, *O. albicans*, *O. vanhoeffeni*, *B. stygius*, *M. erythrocephalus*). OSKA2016 and *F. borealis* are
held out. Seeds of round 0 are the RF00031 hits of at least 15 bits; each round aligns the seeds to the current model,
rebuilds and calibrates it, and keeps as next seeds the training hits above the best hit in 200 random genomic
windows per training genome. The rounds stop when the seed set no longer changes.

Inputs: `work/stage1/pipeline/selenoprofiles/*.three_prime.fa` (stage 1), `data/genomes/` and
`data/references/RF00031.cm` (`make fetch`).

### Reproduce (on a Slurm cluster, in the image built from `environment.yml`)

```bash
sbatch [--partition=<partition>] --cpus-per-task=8 scripts/run_in_container.sh make secis
```

`make secis` runs three targets, which can also be run separately:

| Target | Output | Time on 1 CPU (measured) |
|---|---|---|
| `secis-stage1` | RF00031 scores, `results/tables/stage1_secis_scores.tsv`, `figures/stage1_secis_scores.*` | 4 min |
| `secis-model` | `work/stage2/secis_model/`; the final model, seed alignment, seeds and rounds copied to `results/models/` and `results/tables/stage2_secis_model_rounds.tsv` | about 50 min (5 rounds, then a full calibration of about 25 min) |
| `secis-stage2` | scores of the new model, `results/tables/stage2_secis_scores.tsv`, `figures/stage2_secis_scores.*` | 4 min |

The scripts use as many CPUs as Slurm allocates (`SLURM_CPUS_PER_TASK`). A rebuild with the same inputs should give
the same rounds as `results/tables/stage2_secis_model_rounds.tsv` and the same seed sequences as
`results/models/appendicularia_secis.seeds.fa`: the background windows come from fixed random seeds, and Infernal's
calibration uses its fixed default seed.

### Replicate with other background windows

```bash
sbatch --cpus-per-task=8 scripts/run_in_container.sh \
    make secis-model BACKGROUND_SEED=200 SECIS_OUT=work/stage2/secis_model_seed200
```

Replicates write only to their `SECIS_OUT` directory; the files in `results/` are those of the default run
(`BACKGROUND_SEED=100`). Compare their `rounds.tsv` and `final.seeds.fa` with the default ones.

Other options of `build_secis_model.py` (`--help`): training genomes, seed threshold of round 0, number of background
windows, maximum number of rounds, calibration length of the intermediate rounds.

## Appendicularian profiles of the Sec machinery

`build_profiles.py` (here with the options of `make machinery-profiles`) extends each standard Selenoprofiles profile
of the machinery set (SEPHS2, pstk, SecS,
eEFsec, SBP2, secp43) with appendicularian proteins: the stage 1 predictions labelled "homologue" or "selenocysteine"
that are at least half the median profile length, and the curated transcript proteins of
`curation/machinery_proteins.faa`. They are added with `mafft --add --keeplength`, so the columns and the standard
configuration remain valid, and the profile is rebuilt with `selenoprofiles build`. OdB3 contributes no sequence and
serves as the test genome. The profiles keep the standard names, so that they can replace the standard ones.

`test_machinery_profiles.py` runs Selenoprofiles with these profiles on the nine genomes and compares the
predictions, by label, with stage 1 (`results/tables/stage2_machinery_profiles_test.tsv`).

```bash
sbatch [--partition=<partition>] --cpus-per-task=4 scripts/run_in_container.sh make machinery
```

| Target | Output | Time on 1 CPU (measured) |
|---|---|---|
| `machinery-profiles` | `work/stage2/machinery_profiles/` (profiles; alignments before the build in `alignments/`); copied to `results/models/machinery_profiles/` and `results/models/appendicularia_machinery_profiles.tar.gz` | 2 min |
| `machinery-test` | `work/stage2/machinery_test/` (resumable per genome); `results/tables/stage2_machinery_profiles_test.tsv` | about 75 min for the nine genomes |

| `machinery-fragments` | PSTK predictions discarded by Selenoprofiles in the genomes without a PSTK prediction, classified by blastp against the PSTK found elsewhere and KTI12 (`curation/machinery_outgroups.faa`); `results/tables/stage2_machinery_fragments.tsv` | about 5 min per genome |
| `machinery-census` | census per gene and genome, `results/tables/stage2_machinery_census.tsv` and `.md` (the table of the manuscript) | seconds |

The committed test table was produced with the same profile alignments under the names `<family>_app`, before the
profiles took the standard names; a rerun of `make machinery-test` regenerates it with the current files.

### Absence tests (`machinery_absence.py`, `make machinery-absence`)

For a machinery gene missing from appendicularian genomes, tests whether it is absent or only divergent: the stage 1
predictions of the gene in other genomes are searched with tblastn in each appendicularian genome (and the
*O. dioica* transcriptomes), and each candidate is compared by blastp with the gene's proteins and with related
proteins of other families (outgroups). The margin (best gene hit minus best outgroup hit, in bits) classifies it as
supported (at least 30), ambiguous (10 to 30, to be resolved by phylogeny) or a related protein. Positive controls are
the stage 1 predictions themselves. SECp43 uses TIA1-like proteins as outgroups; PSTK uses KTI12. The scans of the
largest genomes take longest; run them on the cluster (`sbatch --cpus-per-task=4 scripts/run_in_container.sh make
machinery-absence`, from a writable clone).

### Reciprocal best hits (`machinery_rbh.py`, `make machinery-rbh`)

Each machinery prediction (and SECp43) is searched against the human RefSeq proteome (blastp; genes from the RefSeq
annotation), and the best human protein is searched back against the genome (tblastn). Status: RBH (expected human
gene, and the reverse search returns to the prediction), forward (expected gene, but the human protein prefers another
locus), other gene (for example SECISBP2L for SECISBP2) or no hit; margins below 10 bits to the next human gene are
flagged as ambiguous, for phylogeny. Needs `data/genomes/Hsap/protein.faa.gz` (fetched with the human genome).

### Phylogenies: sequences

`machinery_phylogeny.py` combines the profile sequences (at most three per taxonomic group with `--max-per-group 3`;
groups from `taxonomy.py`, below),
the SPS1/SPS2-labelled sequences of Mariotti et al. 2015 (SEPHS), the predictions of the other tunicates,
cephalochordates and vertebrates (stage 1) and of the appendicularians (stage 2), the curated transcript proteins,
and, for PSTK and the SECIS-binding proteins, the RefSeq proteins of the genomes that have them (KTI12 proteins join
the outgroup). `--collect-only` writes the sequence sets without aligning. `scripts/proteins_from_gff.py` builds a
protein file from a genome and its RefSeq GFF3 when the proteins were not downloaded.

### Phylogenies: curation, alignment and trees (`make machinery-phylogeny`)

Before the trees, the sequence sets are curated and every removal is recorded in `results/trees/<family>.removed.tsv`
with its reason: the human stage 1 predictions are left out (the human RefSeq proteins represent human); identical
sequences are kept once, under the genome or RefSeq name; after a first alignment, sequences covering less than half
of the core columns (columns where at least half of the sequences have a residue; `--min-core`) are removed and the
rest realigned, except appendicularian, curated and tunicate sequences, which are kept even when partial.
Alignment: MAFFT L-INS-i. Trimming: ClipKIT gappy mode, removing the columns with at least 90 % gaps (insertions found
in a few sequences, such as retained introns or rare isoforms); 50 % for the SECIS-binding proteins, whose tunicate
predictions often cover only the C-terminal domain. Trees: IQ-TREE 2 with ModelFinder restricted to LG, WAG and JTT
and 1,000 ultrafast bootstraps; PSTK and SEPHS also under LG+C20+F+G (`results/trees/C20/`, from the same trimmed
alignments). Outputs per family in `results/trees/`: `.faa`, `.groups.tsv`, `.removed.tsv`, `.aln`, `.trimmed`,
`.treefile`, `.iqtree`; figures `figures/phylogeny_<family>.*` (full trees) and `_collapsed.*` (display trees,
`figures/phylogeny_collapsed.py`). About 4 hours on 8 CPUs, most of it the C20 tree of SEPHS.

Taxonomic groups (`taxonomy.py`, shared by the tree, long-branch, orthology and figure scripts) are read from the
sequence names: genome code or genus at the start of the name, Ensembl species prefix of a protein identifier and
UniProt species mnemonic of an entry name (`config/taxon_codes.tsv`), then any token that is a known genus
(`config/taxon_groups.tsv`, case-insensitive), then descriptive bacterial names; the rest is "unclassified". The
machinery trees of `results/trees/` (2026-10-19 to 10-24) were built before this module, when only the first token was
looked up and sequences left unclassified escaped the limit of three per group: a rerun at a later commit subsamples
the profile sequences differently (the earlier commits reproduce them).

### Long-branch analyses (`machinery_lba.py`, `make machinery-lba`)

Appendicularian proteins evolve fast, so their position is exposed to long-branch attraction. For each family:
LG+C60+F+G with PMSF (guide: the maximum-likelihood tree); removal of the fastest 10 to 40 % of the sites (IQ-TREE site
rates), each reduced alignment reanalysed; removal of named long branches (`--exclude`, by default the insect-gut fungi
in PSTK); an AU test of the grouping of the appendicularian clade with the other tunicates of the same form; for SEPHS,
two more AU tests (appendicularian selenocysteine clade inside the clade of all tunicate SEPHS, and at the position
expected for an ortholog), under the selected model and under LG+C60+F+G with the PMSF site frequencies. Leaves are
assigned to taxa by `taxonomy.py`; the appendicularian clade is the largest one in the maximum-likelihood tree.
Outputs: `results/trees/lba/` (trees, reports, site rates, `<family>.summary.tsv`), `results/tables/stage2_machinery_lba.tsv`
(all families), and, through `machinery-lba-figures`, the display trees of the PMSF trees
(`figures/phylogeny_<family>_collapsed_pmsf.*`). About 2 hours on 8 CPUs. `make machinery-lba` starts afresh; to
resume after an interruption, or to refresh the summaries after a change of the script, run
`python3 stage2/machinery_lba.py` directly: finished analyses (IQ-TREE report present) are kept, interrupted ones
redone.

### Domain analysis (`machinery_domains.py`, `make machinery-domains`)

Tests whether the two domains of each protein support the same tree. Each trimmed alignment is split at a domain
boundary of the human protein (`curation/machinery_domains.tsv`: UniProt coordinates from InterPro entries of
CATH-Gene3D, SUPERFAMILY and Pfam; UniProt sequences in `curation/machinery_domain_references.faa`, carried to the
human RefSeq sequence of the alignment by a global pairwise alignment): PSTK kinase domain / C-terminal region, SEPHS
N-terminal / C-terminal halves, SECIS-binding protein region upstream of the L7Ae domain / L7Ae domain. Sequences
covering less than 30 % of the columns of either part are left out. For each part: a tree (IQ-TREE, as for the
full-length trees), the supports of the appendicularian clade and of its grouping with the other tunicates, the
constrained trees of `machinery_lba.py`, and AU tests of both domain trees, of the full-length tree (pruned) and of
the constrained trees against the data of that part. Outputs: `results/trees/domains/` (alignments, trees, reports,
constraints), `results/tables/stage2_machinery_domains.tsv` (one row per family, part and tested tree),
`results/tables/stage2_machinery_domain_conflicts.tsv` (clades supported by at least 95 or 90 % in one domain tree that
are incompatible with a clade so supported in the other) and, through `machinery-domain-figures`, display trees
(`figures/phylogeny_<family>_collapsed_domain_<domain>.pdf`). About 2.5 hours on 2 CPUs. As for `machinery_lba.py`,
running the script directly keeps the finished IQ-TREE analyses and only rewrites the tables.

`make machinery-domains-pmsf` repeats the analysis under LG+C60+F+G with PMSF (`--pmsf`): site frequencies estimated
on the maximum-likelihood tree of each part, then the tree, the constrained trees and the AU tests with these
frequencies; the full-length PMSF tree of `results/trees/lba/` is added to the candidate trees. Outputs
`results/trees/domains/<family>.<domain>.pmsf.*`, `results/tables/stage2_machinery_domains_pmsf.tsv`,
`results/tables/stage2_machinery_domain_conflicts_pmsf.tsv` and `figures/phylogeny_<family>_collapsed_domain_<domain>_pmsf.pdf`.
It keeps finished analyses, so it can be resubmitted after an interruption.

### Table MACHINERY (`machinery_table.py`, `make machinery-table`)

Combines the census (`stage2_machinery_census.tsv`), the reciprocal best hits (`stage2_machinery_rbh.tsv`), the
long-branch table and the SEPHS trees into the table of the manuscript, with an orthology column and a row for the
second, non-selenocysteine SEPHS copy (`results/tables/stage2_machinery_table.{md,tsv}`). It reads only `results/` and
runs in seconds; `make machinery-lba` runs it at the end.

## SECIS search in the validated hits (`secis_validated_hits.py`, `make secis-hits`)

Step 1 of the plan of 2026-10-24. The stage 1 predictions with a selenocysteine codon in the six *Oikopleura*
assemblies, left out only on evidence against them from hit validation round 1 (best human hit in an unrelated
family, counting MIEN1 for SELENOW and SEPHS1 for SEPHS2 as related; non-animal nearest relatives in the family tree);
the curated *O. dioica* genes always kept; *B. stygius*, *M. erythrocephalus* and *F. borealis* reported apart. No
new search: the RF00031 and current-model scores of `results/tables/stage1_secis_scores.tsv` and
`stage2_secis_scores.tsv`, with the random windows of each genome and the cysteine homologs as controls, and the
length of each 3' sequence (shorter than 6,000 nt at contig ends). Outputs `results/tables/stage2_secis_hits*.tsv`.

## SECIS model trained on the six *Oikopleura* assemblies (`secis_model_training.py`, `make secis-v2`)

Step 3 of the plan of 2026-10-24. Positive windows: the 2,000 nucleotides after the stop codon of each selenoprotein
gene, shortened at the next selenoprotein coding sequence on the same strand; stop codons from the gold-standard
models of step 2 (OKI2018, OSKA2016, OdB3) and from the stage 1 predictions kept by step 1, extended to the next stop
codon (*O. longicauda*, *O. albicans*, *O. vanhoeffeni*). Negative windows: random genomic windows of the six
assemblies, kept away from the positive windows, and of four non-target tunicates (separate random seeds for training
and evaluation), dinucleotide-shuffled
positive windows, and windows after the stage 1 cysteine homologs. Training bootstraps from RF00031 with one SECIS
element per locus and a 90% identity filter on the seeds; a merged model (RF00031 seed sequences with the same
appendicularian seeds) is built in each fold. Folds: one per family and one per species (the three *O. dioica*
assemblies together); each is scored with its new and merged models, RF00031 and the current model; the evaluation
negatives above the training threshold are listed with their coordinates (candidate SECIS elements). Outputs:
`work/stage2/secis_v2/` (windows, one directory per fold, final model), `results/tables/stage2_secis_v2_*.tsv` and
`results/models/oikopleura_secis.{cm,sto,seeds.fa}`. Needs Infernal and Biopython (container image); about 16 jobs
of training, four at a time (`--jobs`). Run of 2026-10-25: 84 of 111 held-out windows detected in the species folds
at 0.32 evaluation negatives per Mb (RF00031: 27); the final model replaces the first appendicularian model for the
SECIS of genes, with a threshold of 20 bits in the 3′ windows (`docs/plan.md`).

## UGA/SECIS screens (`uga_secis_screen.py`, `make secis-scan`, `make uga-screen`)

Step 4 of the plan of 2026-10-24, independent of protein families. `scan`: genome-wide SECIS scans with
`results/models/oikopleura_secis.cm` (both strands, hits of 15 bits or more) of the six *Oikopleura* assemblies, the
three appendicularian genomes held out of stage 2 (*B. stygius*, *M. erythrocephalus*, *F. borealis*; computed for
the final census) and four non-target tunicates, and of copies shuffled in windows of 1 kb (dinucleotides
preserved) as background;
SECIS searches of the transcript sets: annotated transcripts (extended by 1 kb of genome after their 3′ end, and
300 nt before CDS-only models), sense strand; transcript assemblies (OKI2018 Trinity, OSKA2016 transcripts, OdB3
ESTs), both strands; human and *C. intestinalis* RefSeq mRNAs as controls, also with RF00031. `screen`: the
genome-scan threshold from the shuffled copies (at most one expected false hit in the largest *Oikopleura*
assembly); transcript screen (a frame through a TGA, from an ATG to a stop codon before a SECIS hit of 20 bits or
more; frames overlapping a longer frame in another reading frame left out); comparative screen (TGA/TGA pairs in
tblastx alignments between species of the regions upstream of genome-scan hits, as SGP2; the held-out genomes are
compared with all species but are not partners of the *Oikopleura* assemblies); blastp of the candidates
against the human RefSeq proteome and the gold-standard proteins; known genes from the gold standard, the stage 1
predictions and the annotated selenocysteines of the controls. Tested on a reduced OKI2018 and OdB3 set (seven
gold-standard genes: all found by the transcript screen; the comparative test, OdB3 standing in for another
species, found the five genes present in both).

## Description of the SECIS model (`secis_model_description.py`, `make secis-description`)

Reads repository files only (the model, its seeds, the gold-standard windows after the stop codons, the training
window table). Seeds realigned with `cmalign`: identities over the consensus positions and single-linkage groups at
0.7; per consensus position, residue counts; per base pair, canonical pairs, pair types and mutual information with
a permutation p-value over one seed per group. The same pair statistics for decoys (best hits of the model in 20
dinucleotide-shuffled copies of each gold-standard window, `cmsearch --max`, the 76 best), because a structural
alignment gives complementary bases to sequences that are not SECIS elements. Changes between the orthologous SECIS
hits of the three *O. dioica* assemblies, aligned by sequence to the OKI2018 hit (`--reference-genome`), counted by
three-taxon parsimony and compared with the fraction of new bases that would keep a canonical pair (transition rate
from the other unpaired positions). The SECIS hits of TXNL1.b and KKSU (not in training) are described as the seeds.
cmsearch tables are cached in `work/stage2/secis_model_description/` (the decoy search takes about 6 minutes with
8 CPUs). Figures: `figures/secis_model_figures.py` (also classes the genome-scan hits at the genome threshold,
`results/tables/stage2_secis_scan_classes.tsv`).

## Hits without a selenoprotein gene upstream (`secis_orphan_hits.py`, `make secis-orphans`)

The genome-scan hits at the genome threshold classed "nothing found upstream" or "stage 1 prediction without
selenocysteine" (`results/tables/stage2_secis_scan_classes.tsv`): annotated transcripts at the hit and upstream on
the same strand (genomes with an annotation), blastx of the 5 kb upstream against the human RefSeq proteome and the
gold-standard proteins (best human and best gold-standard hit, each with the stop codons inside the alignment and the residue aligned to each), copies of the hit in its
own genome (blastn, 80% identity over 80% of the hit and 50 nt flanks), and conservation in the other *Oikopleura*
assemblies (best HSP, and a genome-scan hit nearby). Output `results/tables/stage2_secis_orphan_hits.tsv`. Tested on
OKI2018 chromosome 2 and two OdB3 contigs (with the gold-standard proteins standing in for the human proteome).

## Candidate genes (`candidate_genes.py`, `make candidate-genes`)

Step 5 of the plan of 2026-10-24. Candidates: `curation/candidates.tsv` and `curation/candidates.faa` (TXNL1.b; the
KKSU family of step 4). The rule fixed in advance for TXNL1.b, applied to every candidate: TGA at the curated position
and a SECIS of at least 20 bits (Oikopleura model) in the 2 kb after the next stop codon, in the three *O. dioica*
assemblies; then tblastn against the appendicularian, tunicate and amphioxus genomes (codon aligned to the
selenocysteine, codons to the next stop, SECIS after it with both models) and blastp against the RefSeq proteomes.
Tested on OKI2018 chromosome 2 and the OdB3 TXNL1.b contig (contexts of `candidates.tsv` reproduced).

SELENOI (`curation/candidates.faa`), found by the absence tests of step 9, has no curated row: its codons are located
by `candidate_transcripts.py` (run first by `make candidate-genes`). The transcripts that encode it in the three
*O. dioica* transcript sets (rows of `results/tables/stage2_absence_SELENOI.candidates.tsv`) are placed on their
assembly with the minimap2 alignments of `make gold-standard` (`work/stage2/gold_standard/<genome>.transcripts.paf`;
aligned again, these transcripts only, when that file is missing) and the junction placement of
`curation/gold_standard.py`; the codon aligned to the first U of the candidate protein is mapped to the assembly
(`results/tables/stage2_candidate_transcripts.tsv`: codon in the transcript and in the assembly, exon, stop codon,
transcript length after it and SECIS hits there; `stage2_candidate_transcripts_positions.tsv`: the positions used by
the rule). Without such a position, the codon aligned by tblastn is used. Tested on a synthetic gene (five exons,
minus strand, transcript stored in both orientations, SECIS seed in the 3' UTR): codon position and exon reproduced.

## Operons and SECIS elements (`operons.py`, `make operons`)

Step 10. In each *O. dioica* assembly, the loci of the public annotation (coding spans; transcripts with overlapping
coding spans merged; loci overlapping a gold-standard selenoprotein gene replaced by it) are chained into operons
(consecutive loci on the same strand, coding spans at most 500 nt apart; Plessy et al. 2024 used 500 nt between genes).
For each gold-standard gene: neighbours, operon members and position, the SECIS element of the Oikopleura model after
the stop codon and the space left before the next gene on the same strand, transcripts starting with the spliced
leader (3' end GATTGTCTAACAG, read in the transcripts of this study) and where they start, transcripts that run from
the upstream neighbour into the gene, and whether the neighbours are homologous between the three assemblies
(blastp). Background: the same definition for all loci (Fisher's exact test). Tested on a synthetic contig (operon of
three loci, annotated model replaced by the gold-standard gene, SL read in both transcript orientations).
Outputs: `results/tables/stage2_operons.tsv`, `stage2_operons_summary.tsv`.

## Operon partners (`operon_partners.py`, `make operon-proteins`, `make operon-search`, `make operon-partners`)

Do the genes that share an operon with a selenoprotein gene have to do with redox or selenium? Loci and operons as in
`operons.py`; the protein of each annotated locus (longest coding sequence, translated from the genome; `make
operon-proteins`) is searched by DIAMOND (`--more-sensitive`, E <= 1e-5) in the Swiss-Prot database of the OIST
Bioinformatics User Group (release 2020_06, https://github.com/oist/BioinfoUgrp/blob/master/DataBases.md), human
entries only (`--taxonlist 9606`; `scripts/operon_search.sh`, outside the container). Flags of the best hit: from its
protein name (oxidoreductase-type names; seleno-), and, with the UniProt table of `make human-uniprot` (optional, on a
node with internet access), its UniProt keywords (Oxidoreductase, Redox-active center, Antioxidant, Selenium,
Selenocysteine) and directly annotated GO terms (cell redox homeostasis, response to oxidative stress, cellular oxidant
detoxification, antioxidant and peroxidase activity, selenium compound and selenocysteine metabolism). Partners and
immediate neighbours of the selenoprotein genes are compared with the other loci in operons (Fisher's exact test), per
assembly. A legacy blastp fallback (`--human-fasta`) exists for use without the database. Tested on a synthetic contig
(operon of four loci with a gold-standard gene, a second operon, two single loci; four human proteins; mock DIAMOND
output; Fisher's test checked against scipy). Outputs: `results/tables/stage2_operon_partners.tsv`,
`stage2_operon_partners_summary.tsv`; best hits of all loci in `work/stage2/operon_partners/<genome>.loci.tsv`.

## miniprot at the selenocysteine codon (`miniprot_junctions.py`, `make miniprot-junctions`)

The manuscript suggests that miniprot misses SELENOT and SEPHS2 because their codons lie 2 and 1 nt from an intron.
The gold-standard proteins of the three assemblies (selenocysteine as U, and as C) are aligned by miniprot to each
assembly as in the pipeline (`-I --gff` on the indexed genome); for each gene and query, the best alignment over the
gene is compared with the gold-standard model: the pipeline's rule for calling the codon (U aligned to TGA, or U next
after the alignment end with TGA as next codon; code copied from the pipeline's `bin/sec_check.py`), whether the codon
lies in one CDS feature and in frame, and the shift of the alignment's intron that overlaps the gold-standard intron
nearest to the codon. The container has miniprot 0.18; the pipeline's module pins 0.11. Tested on a synthetic gene on
both strands with a mock miniprot output. Output: `results/tables/stage2_miniprot_junctions.tsv`.

## SELENOI junction (`selenoi_junction.py`, `make selenoi-junction`)

The transcripts of the three *O. dioica* assemblies (from `make selenoi-regions`, in `results/sequences/`) are aligned
to the codon regions with minimap2 (`-x splice:hq -un`); for each junction, all placements giving the same transcript,
their splice sites, canonical placements within 15 nt with 1 or 2 mismatches, and the bases before each end of the
removed segment. Output: `results/tables/stage2_selenoi_junction.tsv`.

## SELENOK of the other Oikopleura species and release profiles (`extract_regions.py`, `region_gene_models.py`)

The final SELENOK profile does not find the SELENOK genes that the UGA/SECIS screens found in *O. longicauda*,
*O. albicans*, *O. vanhoeffeni* and *M. erythrocephalus* (`make selenok-diagnose`: no candidate prediction before the
filters, or a fragment). `make selenok-regions` extracts 9 kb around the screen codons (`curation/selenok_loci.tsv`)
to `results/sequences/selenok_regions.fa`; `make selenok-models` (`region_gene_models.py`) builds gene models there:
GeneWise with each gold-standard SELENOK protein, the last exon joined to the exon of the screen codon (donor
GT/GC/GA, acceptor AG, TGA in frame, best global alignment to the query), start codon at most 10 codons upstream
(`results/sequences/selenok_models.faa`, `results/tables/stage2_selenok_models.tsv`). `make stage2-release-profiles`
builds the final profiles plus these proteins (the *M. erythrocephalus* model is held out), the profiles of the
deposit (`make stage2-bundle`, which also carries the final profiles); `make selenok-diagnose
SELENOK_DIAG_BUILDS="final release"` compares the two SELENOK profiles.

## SECIS element of SELENOI (`selenoi_regions.py`, `selenoi_secis_search.py`, `secis_stemloop_screen.py`)

`make selenoi-regions` writes the regions around SELENOI (the three *O. dioica* assemblies, the tblastn loci of the
other appendicularians, tunicates and amphioxus, the stage 1 sites of vertebrates and *C. intestinalis*) and the
*O. dioica* SELENOI transcripts to `results/sequences/selenoi_regions.fa` (with `selenoi_regions.tsv`), so that the
search runs without the genomes. `make selenoi-secis` (`selenoi_secis_search.py`) locates the candidate codon in
each appendicularian region (tblastn of the 45 residues before the U), places the transcripts on their region
(retention or splicing of the segment that holds the codon), and searches the 3' sequences (3 kb; 8 kb for
vertebrates) with the covariance models without filters (`cmsearch --max`, against dinucleotide shuffles), with the
stem-loop screen of `secis_stemloop_screen.py` (local minimum free energy stem-loops, then SECIS features) and core
patterns without folding, and compares five *Oikopleura* species (MAFFT, RNAalifold). The stem-loop screen and the
core patterns are calibrated on the gold-standard 3' sequences (known SECIS elements against nt 1001-2000).

## TXNL1 homologs (`txnl1_homologs.py`, `make txnl1-homologs`)

Whether *O. dioica* TXNL1.b (168 residues, PITH region only, selenocysteine) and human TXNL1 are reciprocal best hits,
what the next hits of TXNL1.b are in its own assembly, and where the thioredoxin part of human TXNL1 (residues 1-117) and
its PITH part (118-289) hit. In each of the three *O. dioica* assemblies: blastp of TXNL1.b against the human
RefSeq proteome (hits grouped by gene, margin to the next gene); blastp of human TXNL1 and of its two parts against
the proteins of the assembly (annotated coding sequences translated from the genome, and the gold-standard proteins)
and tblastn against its genome (colinear HSPs grouped into loci); the same with TXNL1.b itself; each top
subject placed relative to the TXNL1.b gold-standard locus and back-searched against the human proteome (blastp of the
protein, blastx of the locus). Outputs `results/tables/stage2_txnl1_homologs.tsv` and
`stage2_txnl1_homologs_summary.tsv`. Tested on synthetic data (a split TXNL1-like gene, a thioredoxin gene next to
it, a PITH paralog and an unannotated copy).

## TXNL1.a and TXNL1.b (`txnl1_paralogs.py`, `make txnl1-paralogs`)

After `make txnl1-homologs`: *O. dioica* has a full-length TXNL1, TXNL1.a (thioredoxin domain and PITH region;
OKI2018 XP_081633344.1), the reciprocal best hit of human TXNL1, and the selenocysteine copy TXNL1.b, with the PITH
region only (names after the tunicate nomenclature guidelines, Stolfi et al. 2015). With
the OKI2018 copies as queries: (1) in the three assemblies, the residue of the full-length copy and of human TXNL1
aligned to the selenocysteine; (2) tblastn against the appendicularian, tunicate and amphioxus genomes: each site
(colinear HSPs of either query) with its structure (full-length, PITH only), the codon aligned to the selenocysteine,
the codons to the next stop, SECIS hits after it, the annotation and the best human gene; (3) blastp against the
RefSeq proteomes and the annotated *O. dioica* proteins; (4) a tree of the PITH regions (TXNL1 and PITHD1 proteins,
and genome sites without an annotated protein; MAFFT, ClipKIT, IQ-TREE), rooted on PITHD1; (5) the neighbourhood of
the full-length copy (a GPR107-like gene: merged model in OSKA2016; transcripts aligned to both). Outputs
`results/tables/stage2_txnl1_paralogs_{queries,sites,proteins,tree,neighbours}.tsv`, `results/trees/txnl1_pith/`,
`figures/txnl1_pith_tree.{png,pdf}`, and the TXNL1.a models and proteins drawn with TXNL1.b in Figure CANDIDATES
(`figures/inputs/reannotation_TXNL1.b.paralogs.{tsv,fa}`: model, codon aligned to the selenocysteine, drawn span). Tested on synthetic data (two *O. dioica*-like assemblies, one with a merged
model, a genome without annotation holding both copies, and a proteome with TXNL1 and PITHD1).

## Introns of TXNL1.a and TXNL1.b (`txnl1_introns.py`, `make txnl1-introns`)

Coding introns of TXNL1.b (gold standard) and TXNL1.a (public model, from `figures/inputs/reannotation_TXNL1.b.paralogs.*`)
in the three *O. dioica* assemblies: residue after each intron and phase, the TXNL1.b positions carried onto TXNL1.a
through a local alignment; shared when the other gene has an intron of the same phase at the same aligned residue.
Reads repository files only. Output `results/tables/stage2_txnl1_introns.tsv`.

## Position of the selenocysteine codon (`sec_codon_positions.py`, `make sec-positions`)

In each gold-standard model: coding exon of the selenocysteine codon, distance to the nearest coding exon junction and
upstream of the last one, codons to the stop codon; each compared with the sense codons of the same gene (mid-rank
fraction), and summed over the genes of OSKA2016, with all genes and with one gene per family. Reads repository files
only. Outputs `results/tables/stage2_sec_codon_positions.tsv` and `_summary.tsv`.

## Tandem stop codons and SECIS-like hits (`tandem_stops.py`, `make tandem-stops`)

A path to a C-terminal selenoprotein: a gene ending with a TGA followed, within one or two codons, by another stop codon
would gain a selenocysteine and at most two more residues if a SECIS element appeared downstream (as KKSU). For each
annotated genome (one coding sequence per gene): the primary stop codon and the next in-frame stop codons (at +1, +2,
+3, against the per-codon frequency of stop codons at +4 to +33); the best hit of the genome scan of step 4 (and of its
shuffled copy, same coordinates) in the 2,000 nt after the stop codon, by class (TGA with a stop at +1 or +2, TGA
otherwise, TAA or TAG, known selenoproteins); candidates of the path class with a hit at the genome threshold. Tested on
a synthetic genome (both strands, CDS rows with and without the stop codon, a known selenoprotein). Outputs
`results/tables/stage2_tandem_stops_{codons,secis,candidates}.tsv`; per gene in `work/stage2/tandem_stops/`.

## Orthology of the selenoproteins

The approach of the machinery applied to the selenoprotein families that have appendicularian sequences
(`SELENO_TREES` in the Makefile: DIO, GPX, MsrA, SELENOF, SELENOK, SELENOL, SELENOO, SELENOS, SELENOT, SELENOU,
SELENOW, TXNRD), to give a first global view before family-by-family analyses.

Two rounds, set by the Makefile variables `SELENO_PREDICTIONS` and `GOLD`. First pass (2026-10-24): stage 1
predictions of the appendicularians (`SELENO_PREDICTIONS=work/stage1/pipeline/selenoprofiles GOLD=`). Round 2 (step 8,
the default from 2026-10-29): predictions of the stage 2 final run (`work/stage2/pipeline/final`), in which the
predictions of OKI, OSKA and OdB3 that overlap a gold-standard gene of the family are replaced by its gold-standard
protein (leaf `<genome>_gold_<gene>`, recorded in `<family>.removed.tsv`); the test assemblies KUM, AOM and Bar2 are
not in the trees.

1. `make selenoprotein-rbh`: `machinery_rbh.py --families selenoproteins` on the predictions of `SELENO_PREDICTIONS`
   for the genomes it holds and on the stage 1 predictions of the other genomes.
   Each family has its human genes (`SELENOPROTEIN_GENES`: for example GPX1 to GPX8, DIO1 to DIO3, SELENOW and
   SELENOV), possibly none (SELENOL, SELENOU: status "best hit only"); the margin is the score of the best expected
   gene minus that of the best other gene, so that a small margin flags paralogs that best hits do not separate.
   Output `results/tables/stage2_selenoprotein_rbh.tsv`.
2. `make selenoprotein-phylogeny`: `machinery_phylogeny.py` with the same settings as the machinery trees
   (profile sequences, at most three per taxonomic group; stage 1 predictions of the other tunicates,
   cephalochordates and vertebrates; human RefSeq proteins of the family, of a related family as outgroup for
   TXNRD (GSR) and SELENOF (SELENOM), and of a related gene as reference for SELENOW (MIEN1, the cysteine homolog
   that is the best human hit of the SelW-like proteins of tunicates and amphioxus); curation, MAFFT, ClipKIT,
   IQ-TREE), with the appendicularian sequences of the round and `--min-fraction 0.3`, since the profiles of some
   families mix short and long proteins.
   Families without appendicularian sequence are collected but not aligned. Outputs `results/trees/selenoproteins/`
   and `figures/selenoproteins/phylogeny_<family>.pdf`, `phylogeny_<family>_collapsed.pdf`.
3. `make selenoprotein-pmsf`: `machinery_lba.py --analyses pmsf` on these trees (LG+C60+F+G with PMSF, guide: the
   maximum-likelihood tree). Outputs `results/trees/selenoproteins/lba/`, `results/tables/stage2_selenoprotein_lba.tsv`
   (the supports of the largest appendicularian clade only) and `figures/selenoproteins/phylogeny_<family>_collapsed_pmsf.pdf`.
4. `make selenoprotein-orthology`: `selenoprotein_orthology.py` splits the appendicularian sequences of each tree
   into appendicularian clades and reads, for each, its support, its sister group, the human genes of the nearest
   side that contains vertebrate sequences (human RefSeq proteins and vertebrate predictions with their best human
   hit), and the best human hits of its members. Outputs `results/tables/stage2_selenoprotein_orthology.{tsv,md}`.
   Gold-standard proteins are labelled "gold standard", named by their gene, and take the best human hits of the
   predictions they replaced. With `GOLD=` (first pass), the curated genes among stage 1 predictions are named from
   `results/tables/stage1_annotation.tsv`. This is a summary for review, not an assignment: the families that it
   leaves open are then analysed as the machinery was (AU tests, site removal, domains).
5. `make selenoprotein-au` (after `selenoprotein-pmsf`): `machinery_lba.py --analyses pmsf hypotheses` on the families
   of `SELENO_AU` (DIO, SELENOF, SELENOK, SELENOO, SELENOS, SELENOT, SELENOW), with the positions of `au_hypotheses.tsv`: for each, a
   constraint tree (clades given as expressions of leaf sets: the appendicularian clade, appendicularian clades
   holding named leaves, the smallest clade holding named leaves, taxonomic groups, regular expressions; see the
   script), a constrained search, and an AU test of the constrained trees against the unconstrained one, under the
   model of the maximum-likelihood tree and under LG+C60+F+G with the PMSF site frequencies. Tested positions: the
   appendicularian SELENOF, SELENOK, SELENOO, SELENOS and SELENOT with the other tunicates and inside the chordate
   (SELENOF, SELENOK, SELENOS, SELENOT) or animal (SELENOO) clade; SELENOW_PAR outside the smallest clade of the human SELENOW and SELENOV, and with the
   appendicularian SELENOW of chromosome 1, of chromosome 2, or with the mammalian SELENOV; the animal sequences of
   the smallest clade of the human and gold-standard DIO as one clade, without the non-animal DIO. Where the reference
   group of a hypothesis is not a clade of the trees, a control places the group as a clade with the appendicularian
   clade outside it: a rejection of both means that the data reject the monophyly of the group, not the position of
   the appendicularians. All constrained trees of a family enter one AU test (redone when the set changes). Outputs `results/trees/selenoproteins/lba/<family>.hypotheses.tsv`
   (resolved leaf sets, and whether each clade is in the maximum-likelihood tree), `.hyp_<id>.*` (constrained
   trees), `.hypotheses_au.iqtree` (AU tests), and rows "AU test: ..." in `results/tables/stage2_selenoprotein_lba.tsv`.
6. `make selenoprotein-table`: `selenoprotein_table.py`, Table SELENO-ORTHOLOGY of the manuscript (one row per
   gold-standard gene: best human hits of the predictions it replaced, its appendicularian clade, sister group and
   nearest vertebrate genes). Outputs `results/tables/stage2_selenoprotein_table.{tsv,md}`.
7. `make selenoprotein-absence` (Deigo): `machinery_absence.py` for the families found in other tunicates but in no
   appendicularian (SELENON, SELENOH, MSRB, FrnE, SELENOL, SELENOI): the selenocysteine forms of stage 1 (other
   genomes) searched in each appendicularian genome and transcriptome by tblastn, each candidate scored against them
   and against human proteins of related families (CALM1 for the EF-hand of SELENON; SELENOW, SELENOV, SELENOT for
   SELENOH; MSRB2 and MSRB3 for MSRB, whose references are the MSRB1-type selenocysteine forms only; GSTK1 for FrnE;
   TXN for SELENOL; CEPT1 and CHPT1 for SELENOI); for each candidate, the residues aligned to the selenocysteine of
   its best reference (TGA read as U; transcripts read through TGA codons) and that protein. Outputs
   `results/tables/stage2_absence_<family>.{summary,candidates,controls}.tsv`.

SELENOK, missed by the stage 1 settings in every appendicularian genome, is in the round 2 trees (stage 2 final run
and gold-standard proteins).

## Stage 2: profiles, runs and benchmark (steps 6 and 7)

Goal: settings and data files of the pipeline for annotating the selenoproteins of new appendicularian assemblies,
measured against the gold standard on assemblies the profiles have not seen.

### Inputs (`make stage2-inputs`, `make stage2-truth`)

`build_profiles.py` builds Selenoprofiles profiles from the standard ones and appendicularian proteins (sources and
rules in `--help`): gold-standard proteins (the profile of a gene is its name up to the first `_`, without dots:
TXNL1.b goes to the profile TXNL1b, because Selenoprofiles reads the family of a prediction up to the first dot),
stage 1 machinery predictions, validated stage 1 selenoprotein predictions and curated proteins. Standard families
keep their columns (`mafft --add --keeplength`); a member whose selenocysteine falls outside the selenocysteine
column of the standard alignment is removed. KKSU and TXNL1b have no standard profile: their members are aligned
with MAFFT L-INS-i and their configuration replaces the AWSI filter, which needs more sequences, by
`x.label != 'pseudo'`. *B. stygius*, *M. erythrocephalus* and *F. borealis* are held out of every build (decision of
2026-10-24). Configuration keys listed in `stage2/profile_config_overrides.tsv` replace those of the standard
configuration, with the reason: the SELENOS filter no longer requires the vertebrate motif (5 of 9 residues) for
predictions labelled selenocysteine, because it removed the *O. dioica* SELENOS (`make stage2-diagnose`).

| Build (`work/stage2/profiles/`) | Selenoprotein members | Machinery members (stage 1) | Used by runs |
|---|---|---|---|
| `E0` | gold standard of OKI2018, OSKA2016, OdB3 | the same three assemblies | `eval`, `eval_miniprot` |
| `LOO_<g>` | the same without assembly `<g>` | the same without `<g>` | `loo_<g>` |
| `final` | E0 + validated stage 1 predictions of *O. longicauda*, *O. albicans*, *O. vanhoeffeni* + *O. albicans* KKSU | the six *Oikopleura* assemblies | `final` |

The final profiles are copied to `results/models/stage2_profiles/` (alignments and configurations) and archived with
their `.profile_data` in `work/stage2/appendicularia_selenoprofiles.tar.gz` (for Zenodo); the members of every
build are listed in `results/tables/stage2_profile_members.tsv`. `cm_thresholds.py` writes the *Oikopleura* SECIS
model with a gathering threshold of 20 bits (`work/stage2/models/oikopleura_secis.cm`), so that the pipeline
searches it with `--cut_ga` (reporting and inclusion threshold), without `--secis_min_bits`.

`truth_sets.py` writes the reference sets (`results/tables/stage2_truth_sets.tsv`, one table per genome in
`work/stage2/truth/`): the gold standard of the three *O. dioica* assemblies; its models lifted, with minimap2, to
three new assemblies of the same cryptic species from Plessy et al. 2024 (KUM-M3-7f, Okinawa; AOM-5-5f, Osaka;
Bar2_p4, Barcelona; `config/genomes.tsv`, role `test`), with the status of each lifted model and the codon read
at the lifted selenocysteine position (a segment end that falls in an indel is placed from the nearest aligned base
within 30 nt; regions and alignments are kept in `work/stage2/truth/`); and, for the other *Oikopleura* species, the
stage 1 predictions kept by step 1 and the KKSU sites of the comparative screen. The benchmark counts every
reference whose codon is a TGA on the assembly, including lifted models disrupted by an indel of the assembly.

### Pipeline runs (`stage2/run/`)

Each run has a params file `stage2/run/<run>.yaml` and a samplesheet; `stage2/run/run.sh <run>` launches the
pipeline (release 0.2.0) with it. Stage 2 settings (`stage2/run/stage2.config`, added by `run.sh`): Selenoprofiles
3′ sequences of 2,000 nt (the window length on which the 20-bit threshold was set), the *Oikopleura* SECIS model at
its gathering threshold, the appendicularian profiles (`--selenoprofiles_overrides`), and the profiles column
`"eukarya,machinery,KKSU,TXNL1b"` (quoted: the field contains commas).

| Run | Profiles | Genomes | Purpose |
|---|---|---|---|
| `default` | standard; Rfam RF00031 | KUM, AOM, Bar2 | stage 1 settings on the new assemblies (stage 1 itself covers the others) |
| `eval` | E0 | 3 gold, 3 new, Olon, Oalb, Ovan, Bsty, Mery, Fbor | stage 2 settings |
| `eval_miniprot` | E0 | as `eval` | miniprot `-C 0.5` (`stage2/run/miniprot_C05.config`): exons at GA-AG introns |
| `loo_OKI`, `loo_OSKA`, `loo_OdB3` | LOO_<g> | `<g>` and the new assembly of its species | profiles without the species tested |
| `final` | final | as `eval` | census (after the benchmark) |

### Benchmark (`make stage2-benchmark`)

`benchmark.py` compares each run present with the reference sets: per reference gene, the Sec site at the reference
codon (Selenoprofiles label, miniprot support), its family, the labels of predictions overlapping the gene when the
codon is missed, the exon structure of the prediction (identical introns; coding bases shared with the reference)
and the SECIS score; per run and genome, the counts, and the Sec sites outside the reference set. Outputs:
`results/tables/stage2_benchmark_{genes,summary,other_sites}.tsv`. Stage 1 is read from `work/stage1/pipeline`.
The target also copies the parameters, software versions and execution trace of each run to
`results/pipeline_info/stage2_<run>_*`.
The fair tests are the gold-standard assemblies in the `loo_*` runs, and the new assemblies in all runs; the gold
assemblies in `eval` are training data (upper bound).

### Filter diagnosis (`make stage2-diagnose`)

`diagnose_refiltering.py` runs Selenoprofiles with one profile on one genome (jobs `genome:family:build`, build =
`standard` or a directory of `work/stage2/profiles/`) and, through a pre-filtering action (`refiltering_hook.py`),
reports for every candidate prediction the measures used by the filters (coverage, identity with the profile, AWSI
z-scores, selenocysteine alignment, number of profile sequences) and the value of `p2g_filtering`, `p2g_refiltering`
and each of its `and` clauses: `results/tables/stage2_refiltering_diagnosis.tsv` (column
`refiltering_clauses_false`: the clauses that remove the prediction).

### Reproduce

```bash
ml singularity
sbatch [--partition=<partition>] scripts/run_in_container.sh scripts/fetch.sh    # adds KUM, AOM, Bar2 (Zenodo)
sbatch [--partition=<partition>] --cpus-per-task=4 scripts/run_in_container.sh make stage2-inputs stage2-truth
# one after the other, from the repository root (eval_miniprot right after eval reuses its Selenoprofiles tasks)
for run in default eval eval_miniprot loo_OKI loo_OSKA loo_OdB3; do
    NXF_WORK=<scratch>/stage2 PROFILE=<profile> stage2/run/run.sh $run || break
done
scripts/run_in_container.sh make stage2-benchmark
```

### Recommended settings (`stage2/recommended/`, `make stage2-bundle`)

`stage2/recommended/README.md` describes the settings chosen from the benchmark for new appendicularian assemblies
(final profiles, SECIS model with its threshold, 2,000-nt 3′ sequences, the profiles column; miniprot with the
pipeline's options). `make stage2-bundle` assembles them with the data files into
`work/stage2/appendicularia_selenoprotein_annotation.tar.gz`, for Zenodo.

The expected selenoproteome of *O. dioica* for the pipeline's expectation check, with the new profile names, is in
`stage2/run/expectation_oikopleura_dioica.csv`.


## Census (step 9)

Census of the selenoprotein genes of the nine appendicularian genomes (and of the three new *O. dioica* assemblies),
from the stage 2 final run, with a contamination screen and the list of candidates for the strict screen of
horizontal gene transfer (decision of 2026-10-24 in `docs/plan.md`: contamination is the default explanation of
foreign-looking sequences; one HGT screen, with no false positive allowed).

1. `make stage2-collect` (Deigo, minutes): `stage1/collect.py --prefix stage2_final` on `work/stage2/pipeline/final`
   gives `results/tables/stage2_final_{sec_sites,genes_summary,trnasec,annotation,predictions}.tsv` (the
   selenocysteine codons now carry the Selenoprofiles prediction identifier, and `predictions` lists every prediction
   with its coding exons); `machinery_census.py` on the final-run predictions gives
   `results/tables/stage2_final_machinery_census.{tsv,md}` (fragments: those of the machinery analysis).
2. `make census-contigs` (Deigo): `census_contigs.py`, for each contig carrying a prediction or a selenocysteine
   codon, its length, GC content and open reading frame density (fraction of bases in open reading frames of at
   least 100 codons, six frames: bacterial DNA is mostly coding); on a window of 2, 5 or 10 kb centred on its
   predictions, GC content, distance of the tetranucleotide composition (centred log-ratio of canonical
   tetranucleotide counts) to the median of random windows of the same genome and length, and open reading frame
   density, as percentiles and robust z-scores; the same on the 2-kb window (in fragmented assemblies the longer
   reference windows come only from the few long contigs, which may themselves be contaminants). Tiara classes when
   a Tiara output is present (not run: decision of 2026-10-30 in `docs/plan.md`). Per genome: length, contigs, N50,
   and for each window length the fraction of the assembly in contigs at least that long and the medians of the
   reference windows. Outputs `results/tables/stage2_census_contigs.tsv`, `stage2_census_genomes.tsv`.
3. `make census` (seconds, from `results/` only): `census.py` merges the selenocysteine codons into loci, decides for
   each locus whether it is counted (rules in the script: left out only on evidence against it, contig classified
   non-eukaryotic by Tiara, best human hit in another family, or composition outlier or ORF-dense contig with
   non-animal nearest relatives; gold-standard genes always counted), and records the flags (composition and GC outliers, nearest
   relatives in the round 2 ML and PMSF trees, no SECIS hit, codon near the contig end, miniprot only). The counts
   are compared with the other tunicates, the cephalochordates and the vertebrates (stage 1 codons, same locus
   definition, no screen). The strict HGT screen lists every appendicularian leaf of the round 2 trees whose nearest
   relatives outside the appendicularians, or outside the tunicates (a transfer before the split of the tunicate
   lineages, like that of cellulose synthase), are all bacterial, archaeal or fungal, with its contig checks and a
   status: contamination (composition, ORF density; Tiara if run), other family (best human hit elsewhere: the tree lacks the animal
   members of its family, as for the thioredoxins found by the DIO profile), family not confirmed (no human hit),
   placement not stable under PMSF, candidate without intron, or candidate for the tests. An old transfer is
   expected to match its genome in composition; its evidence is its presence at the same place of the tree in
   several species and assemblies; the column of foreign sequences shows the donor lineage, to compare between
   candidates (an ancestral symbiont predicts one donor lineage; bacteria taken up with food predict several). Outputs
   `results/tables/stage2_census_loci.tsv`, `stage2_census.{tsv,md}`, `stage2_hgt_screen.tsv`. Each locus is also
   assigned to the O. dioica gene of its clade (gold-standard gene, lifted model in KUM, AOM and Bar2, or the
   gold-standard genes of the smallest appendicularian clade holding it in the round 2 tree), and the loci of the
   step 4 screens that the final run did not predict are added without gene model (SELENOK of O. longicauda,
   O. albicans, O. vanhoeffeni and M. erythrocephalus, for example): `stage2_census_genes.{tsv,md}` (Table CENSUS).
