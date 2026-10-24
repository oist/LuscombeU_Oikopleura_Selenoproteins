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

`build_machinery_profiles.py` extends each standard Selenoprofiles profile of the machinery set (SEPHS2, pstk, SecS,
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
| `machinery-profiles` | `work/stage2/machinery_profiles/`; copied to `results/models/machinery_profiles/` and `results/models/appendicularia_machinery_profiles.tar.gz` | 2 min |
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

`machinery_phylogeny.py` combines the profile sequences (at most three per taxonomic group with `--max-per-group 3`),
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

### Long-branch analyses (`machinery_lba.py`, `make machinery-lba`)

Appendicularian proteins evolve fast, so their position is exposed to long-branch attraction. For each family:
LG+C60+F+G with PMSF (guide: the maximum-likelihood tree); removal of the fastest 10 to 40 % of the sites (IQ-TREE site
rates), each reduced alignment reanalysed; removal of named long branches (`--exclude`, by default the insect-gut fungi
in PSTK); an AU test of the grouping of the appendicularian clade with the other tunicates of the same form; for SEPHS,
two more AU tests (appendicularian selenocysteine clade inside the clade of all tunicate SEPHS, and at the position
expected for an ortholog), under the selected model and under LG+C60+F+G with the PMSF site frequencies. Leaves are
assigned to taxa by genus or genome code; the appendicularian clade is the largest one in the maximum-likelihood tree.
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
`results/models/oikopleura_secis.cm` (both strands, hits of 15 bits or more) of the six *Oikopleura* assemblies and
four non-target tunicates, and of copies shuffled in windows of 1 kb (dinucleotides preserved) as background;
SECIS searches of the transcript sets: annotated transcripts (extended by 1 kb of genome after their 3′ end, and
300 nt before CDS-only models), sense strand; transcript assemblies (OKI2018 Trinity, OSKA2016 transcripts, OdB3
ESTs), both strands; human and *C. intestinalis* RefSeq mRNAs as controls, also with RF00031. `screen`: the
genome-scan threshold from the shuffled copies (at most one expected false hit in the largest *Oikopleura*
assembly); transcript screen (a frame through a TGA, from an ATG to a stop codon before a SECIS hit of 20 bits or
more; frames overlapping a longer frame in another reading frame left out); comparative screen (TGA/TGA pairs in
tblastx alignments between species of the regions upstream of genome-scan hits, as SGP2); blastp of the candidates
against the human RefSeq proteome and the gold-standard proteins; known genes from the gold standard, the stage 1
predictions and the annotated selenocysteines of the controls. Tested on a reduced OKI2018 and OdB3 set (seven
gold-standard genes: all found by the transcript screen; the comparative test, OdB3 standing in for another
species, found the five genes present in both).

## Orthology of the selenoproteins

The approach of the machinery applied to the selenoprotein families that have appendicularian sequences
(`SELENO_TREES` in the Makefile: DIO, GPX, MsrA, SELENOF, SELENOK, SELENOL, SELENOO, SELENOS, SELENOT, SELENOU,
SELENOW, TXNRD), to give a first global view before family-by-family analyses.

1. `make selenoprotein-rbh`: `machinery_rbh.py --families selenoproteins` on the stage 1 predictions of all genomes.
   Each family has its human genes (`SELENOPROTEIN_GENES`: for example GPX1 to GPX8, DIO1 to DIO3, SELENOW and
   SELENOV), possibly none (SELENOL, SELENOU: status "best hit only"); the margin is the score of the best expected
   gene minus that of the best other gene, so that a small margin flags paralogs that best hits do not separate.
   Output `results/tables/stage2_selenoprotein_rbh.tsv`.
2. `make selenoprotein-phylogeny`: `machinery_phylogeny.py` with the same settings as the machinery trees
   (profile sequences, at most three per taxonomic group; stage 1 predictions of the other tunicates,
   cephalochordates and vertebrates; human RefSeq proteins of the family, and of a related family as outgroup for
   TXNRD (GSR) and SELENOF (SELENOM); curation, MAFFT, ClipKIT, IQ-TREE), with the stage 1 predictions of the
   appendicularians and `--min-fraction 0.3`, since the profiles of some families mix short and long proteins.
   Families without appendicularian sequence are collected but not aligned. Outputs `results/trees/selenoproteins/`
   and `figures/selenoproteins/phylogeny_<family>.pdf`, `phylogeny_<family>_collapsed.pdf`.
3. `make selenoprotein-pmsf`: `machinery_lba.py --analyses pmsf` on these trees (LG+C60+F+G with PMSF, guide: the
   maximum-likelihood tree). Outputs `results/trees/selenoproteins/lba/`, `results/tables/stage2_selenoprotein_lba.tsv`
   (the supports of the largest appendicularian clade only) and `figures/selenoproteins/phylogeny_<family>_collapsed_pmsf.pdf`.
4. `make selenoprotein-orthology`: `selenoprotein_orthology.py` splits the appendicularian sequences of each tree
   into appendicularian clades and reads, for each, its support, its sister group, the human genes of the nearest
   side that contains vertebrate sequences (human RefSeq proteins and vertebrate predictions with their best human
   hit), and the best human hits of its members. Outputs `results/tables/stage2_selenoprotein_orthology.{tsv,md}`.
   This is a summary for review, not an assignment: the families that it leaves open are then analysed as the
   machinery was (AU tests, site removal, domains).

Known gap: SELENOK was not predicted by the stage 1 settings in any appendicularian genome, so that its tree, if
any, lacks the *O. dioica* protein until the stage 2 re-calling or a curated protein from transcripts.
