# Curation

Baseline of the three *O. dioica* assemblies, written as data:

| File | Contents |
|---|---|
| `baseline_sec_genes.tsv` | 20 Sec genes per assembly (18 curated before stage 2, and TXNL1.b and KKSU added on 2026-10-26 after step 5; TXNL1.b, the selenocysteine paralog of TXNL1, was called TXNL1 until 2026-10-27): Sec codon (lowest coordinate of the TGA), strand, Sec context of the protein alignment, status (`sec`, `uncertain`), notes on assembly errors |
| `machinery_pstk.tsv` | PSTK protein alignment spans (RefSeq model truncated in OKI2018) |
| `machinery_proteins.faa` | curated *O. dioica* machinery proteins from transcripts, added to the appendicularian Selenoprofiles profiles (header: family, genome, source): PSTK from the OKI2018 Trinity transcriptome and from the OSKA2016 transcript assembly; SEPHS2, with its Sec, from the OSKA2016 transcript assembly |
| `candidates.tsv`, `candidates.faa` | candidate selenoprotein genes examined in step 5: TXNL1.b (a paralog of thioredoxin-like protein 1, TXNL1, limited to its PITH region) and the family found by the UGA/SECIS screens of step 4 (provisional name KKSU, after its C-terminus), with their codons in each assembly (both added to the baseline on 2026-10-26); their proteins (U for the candidate selenocysteine), queries of `make candidate-genes` |

The Sec codons were placed by aligning curated reference proteins and checking the Sec context; transcripts and
cross-assembly comparison support them (see the stage 1 tables). Assembly errors are recorded, not corrected.

`check_baseline.py` (run by `make curation`) reports, for each curated Sec codon, which stage 1 tools found it, the
RF00031 score and the transcript support: `results/tables/baseline_vs_stage1.tsv` and `_summary.tsv`.

`mariotti2015_sps.faa`: the selenophosphate synthetase sequences of Mariotti et al. 2015 (Genome Research,
doi:10.1101/gr.190538.115, Supplemental Material S8: summary set, vertebrates, insects), named by species and type
(SPS1 or SPS2) as in the paper; used as labelled references in the SEPHS tree (stage2/machinery_phylogeny.py). The
supplement gives codon alignments: they are translated codon by codon (gap codons skipped, stop codons, including
the selenocysteine TGA, read as X).
Header: `SEPHS2|<species>.<type>|<source>`.

`gold_standard.py` (run by `make gold-standard`; step 2 of the plan of 2026-10-24) builds gene models of the curated
selenoprotein genes of OKI2018, OSKA2016 and OdB3. The transcripts of each assembly are aligned with minimap2
(`-x splice:hq -un`) and their junctions placed by the script (fewest mismatches, plus a splice-site cost from the
intron classes of *O. dioica*: GT-AG 87.6%, GA-AG 9%, GC-AG 2.7%, GG-AG 0.7%; Denoeud et al. 2010). Candidate
models are exon chains from the supporting transcripts, the annotated mRNAs and the stage 1 prediction; each is read
on the genome sequence (frame through the curated TGA read as selenocysteine), so that GFF3 and protein agree with
the assembly. The model is chosen by the curated context, length, completeness and transcript support; among
candidates equal on these criteria except length, by similarity to the same gene in the other assemblies; an
annotated mRNA whose coding sequence is another protein's (operons) is not preferred as annotated chain; models
split at the TGA by the annotation are also joined into one candidate; a chain contradicted by most supporting
transcripts is ranked last; 3' UTR exons of read-through transcripts into a downstream gene are left out, and UTRs stop at the coding
sequence of another selenoprotein gene of the table (one transcript holding two genes). When a
transcript differs from the assembly by indels, its protein is used only if it is more similar to the same gene in
the other assemblies; the model is then mapped from the transcript and marked with the GFF3 attribute `exception`.
The table reports splice sites (non-canonical ones with all equivalent placements), differences between transcripts
and assembly, UTRs, the treatment of the gene by the public annotation and the change it needs. Outputs:
`curation/gold_standard/<genome>.gff3` (with `transl_except` for the selenocysteine), `<genome>.candidate_models.tsv`
(all candidate models and transcript frames, and the annotated models near the codon, with splice sites and annotated
CDS rows, for review), `<genome>.proteins.faa`,
`<genome>.utr3.fa`, `<genome>.after_stop.fa`, and `results/tables/gold_standard.tsv`.

`machinery_anchors.py`, `machinery_projections.py` and `gold_standard.py --anchor sense` (run by
`make gold-standard-machinery`; decisions of 2026-11-01) build gene models of the machinery genes (PSTK, SEPSECS,
EEFSEC, the SECIS-binding protein SBP2, and SEPHS_b, the SEPHS copy without selenocysteine; SEPHS2 is in the
selenoprotein gold standard) in OKI2018, OSKA2016, OdB3, *O. longicauda*, *O. albicans* and *O. vanhoeffeni*, so that
the machinery trees use final protein sequences. `machinery_anchors.py` writes `machinery_genes.tsv` (one row per
machinery prediction of the final stage 2 run, with an anchor codon in the middle of its coding sequence; predictions
labelled pseudo, fragments shorter than half the longest prediction of the gene, and predictions overlapping a
selenoprotein gold-standard gene left out) and `machinery_references.faa` (human RefSeq proteins of each gene). A
first pass builds the *O. dioica* models around the anchors, from transcripts, annotations and predictions, ranking
the candidates by the sum of their local alignment scores to the human proteins and, in the final pass, to the
first-pass proteins of the other assemblies (a fragment scores less, a fusion no more): candidates within 5% of the best
score come first, then those without unsupported splice sites, complete ones, those with transcript support, those
not contradicted by a transcript, and higher score. `machinery_projections.py` aligns the human
proteins and the first-pass proteins to the six genomes with miniprot; the final pass adds these alignments as
candidate models and the first-pass proteins of the other genomes as references. Outputs:
`gold_standard_machinery/<genome>.gff3` (no `transl_except`), `<genome>.proteins.faa`,
`<genome>.candidate_models.tsv` (since 2026-11-03 with the insertions, deletions and mismatches of each supporting
transcript against the assembly, drawn in the machinery panels), and `results/tables/gold_standard_machinery.tsv` (column `anchor_codon` instead of
`sec_codon`; `reference_score`: the score and the reference it was reached with).

`machinery_overrides.tsv`: exon chains of machinery genes set by the curator, with the evidence, where no candidate
of the automatic procedure holds the right chain (`gold_standard.py --overrides`; the model is read on the genome along
this chain, and the notes of the gold-standard table say so; its UTRs come from the transcripts that agree with its
first or last two exons, even when they differ elsewhere, for example by a retained intron). OdB3 PSTK (2026-11-03): the GenBank exons 1 to 3, with the
38-nt third exon that the prediction skipped, and the last intron of the ESTs (AT-AG; first set at a GG-AG site 3 nt
away, which added a residue absent from OKI2018 and OSKA2016; corrected when the transcript tracks of Figure PSTK
showed the difference). *O. albicans* PSTK (2026-11-03): the prediction read into its last intron (2 kb) and missed the
last exon, found by tblastn in `results/gene_regions/PSTK` (no transcripts).

`gene_regions.py`: genomic regions of one gene in several assemblies (gold-standard models, regions aligned with
minimap2 to assemblies of other individuals, an annotated protein, contigs given by hand), written to
`results/gene_regions/<gene>/` for checks by hand (`make pstk-regions`).
