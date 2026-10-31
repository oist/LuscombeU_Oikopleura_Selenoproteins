# RefSeq update for the OKI2018 selenoprotein genes: plan

Draft of 2026-10-31; nothing started. After the preprint, we send NCBI RefSeq an evidence package that corrects the
models of the 20 selenoprotein genes of OKI2018 (GCF_907165135.1), and we point the next annotation run to the
appendicularian profiles and SECIS model. The OSKA2016 and OdB3 annotations are optional.

## Why

- Of the 20 gold-standard genes of OKI2018, three RefSeq models declare the selenocysteine (DIO_PAR, GPX_PAR, SELENOF)
  and six genes have no RefSeq model at all (manuscript, "Gold-standard gene models in three *O. dioica*
  assemblies"). Analyses that start from the annotation miss or truncate most of them.
- RefSeq has a dedicated selenoprotein procedure, which annotates 77% of vertebrate selenoprotein genes correctly and
  fails mostly in families with a C-terminal selenocysteine (Ticó et al. 2026). In OKI2018 it declared 3 of 20.
- Ticó et al. (2026) argue that selenoproteins must be corrected in the public annotations that researchers use, and
  propose to feed selenocysteine codons predicted by Selenoprofiles to the annotation pipelines.

## What needs to change

Changes to the current RefSeq models, from the `annotation_fix` column of `results/tables/gold_standard.tsv`:

| Change | Genes | Number |
| --- | --- | --- |
| none: selenocysteine declared, model as the gold standard | GPX_PAR, SELENOF | 2 |
| declare the selenocysteine (`transl_except`) | GPX_chr1a | 1 |
| move an intron and change the coding sequence | DIO_PAR (also a conserved N-terminal extension of 58 residues), DIO_chr1 | 2 |
| change the coding sequence on the annotated exon chain | GPX_chr1b, GPX_chr2, SELENOS (exon chain extended), SELENOT, SELENOW_chr2, SEPHS2, TXNL1.b | 7 |
| new exon chain over an existing model | SELENOO, TXNRD | 2 |
| new gene model | KKSU, SELENOK, SELENOW_PAR, SELENOW_chr1a, SELENOW_pair1, SELENOW_pair2 | 6 |

Three genes cannot be encoded on the assembly as it is: DIO_chr1 (7 bases of the transcript coding sequence absent from
the assembly, 6 assembly bases absent from the transcript), SELENOO (5 bases absent) and SELENOS (13 bases absent).
RefSeq can represent such genes as corrected models ("PREDICTED: LOW QUALITY PROTEIN", bases inserted relative to the
genome; RefSeq FAQ), as it already does for DIO_chr1, whose current corrected model starts 4 codons after the TGA.
The alternative, a new version of the OKI2018 assembly with these bases corrected, would be submitted by the assembly
submitters and would start a new annotation; it is not part of this plan.

## Evidence package

One package for all 20 genes, generated from the repository (script to write: `stage2/refseq_package.py`,
`make refseq-package`; outputs in `results/refseq/`):

- a table with one row per gene: RefSeq gene and protein accessions, change needed, selenocysteine codon (RefSeq
  chromosome accession and position), `transl_except` location, protein length, supporting OKI2018 transcripts,
  SECIS element (position after the stop codon, score of the *Oikopleura* model), best human hit and orthology, and
  conservation in OSKA2016 and OdB3;
- the gold-standard models in GFF3 on the RefSeq chromosome accessions (NC_137291.1 to NC_137294.1 for these
  genes), and their proteins with U at the selenocysteine;
- the corresponding pages of Figures REANNOTATION and CANDIDATES, which draw each RefSeq model against the gold
  standard;
- links to the preprint, to the OKI2018 transcriptome (Zenodo, doi:10.5281/zenodo.23472675) and to the deposits of
  the profiles and SECIS model.

## Sending it

- Channel: the Gene and RefSeq Feedback form, which RefSeq names for problems in the sequence or annotation of a
  record; "suspected errors in mRNA or protein sequences may be corrected, when supported by evidence" (RefSeq FAQ).
- First message: a summary (20 genes, 18 needing a change, the table above), the package (Zenodo deposit) and the
  preprint, with the question of how RefSeq prefers to receive the corrections (one report per gene, or the batch).
- The model records (XM_, XP_) come from the NCBI Eukaryotic Genome Annotation Pipeline, which re-annotates an
  organism "when new evidence is available (e.g. RNA-Seq) or when a new assembly is released" (EGAP description).
  The message therefore also offers the appendicularian profiles and SECIS model for the next run, as Ticó et al.
  (2026) propose for Selenoprofiles predictions.

## Evidence for the next annotation run

- The RNA-seq reads behind the OKI2018 transcriptome (three stages, ENA ERX4506646 to ERX4506648; Bliznina et al.
  2021) are public, so the pipeline probably had transcript evidence already; the models fail at the UGA codon, not for
  lack of transcripts. To check in the annotation report of GCF_907165135.1 (evidence list and annotation release).
- The unit's developmental RNA-seq time course of Okinawa animals would add evidence for the RefSeq annotation if it
  is deposited in the SRA; not planned here.
- A Third Party Annotation (TPA) submission of the gold-standard models to the INSDC is possible only with the
  published paper; to decide after peer review.

## Other annotations (optional)

- OSKA2016 (ANISEED): send the GFF3 of the gold-standard models to the ANISEED curators.
- OdB3 (GenBank/ENA records of the submitters): low priority; the gold-standard models stay available in the
  repository and on Zenodo.

## Steps

1. Write `stage2/refseq_package.py` and build the package; check each row against Figure REANNOTATION.
2. Deposit the package on Zenodo with the other deposits of the paper (Data and code availability), at the release
   of the preprint.
3. Send the first message through the Gene and RefSeq Feedback form; follow the procedure RefSeq asks for.
4. Offer the profiles and SECIS model for the next annotation run (same message, or the contact RefSeq gives).
5. When a new RefSeq annotation of OKI2018 appears, rerun the annotation check of the pipeline on it and count the
   genes with the selenocysteine declared in a correct model (3 of 20 today).

## Questions to settle

- [ ] Who sends the reports, and under which name.
- [ ] The three genes with assembly errors: corrected models in RefSeq, or a new version of the OKI2018 assembly?
- [ ] Include OSKA2016 (ANISEED) and OdB3?
- [ ] Send before or after peer review (the plan assumes after the preprint).

## Sources

- Ticó M et al. 2026. Overcoming the widespread flaws in the annotation of vertebrate selenoprotein genes in public
  databases. PLoS Comput Biol. [doi:10.1371/journal.pcbi.1013885](https://doi.org/10.1371/journal.pcbi.1013885)
- RefSeq Frequently Asked Questions, NCBI Bookshelf [NBK50679](https://www.ncbi.nlm.nih.gov/books/NBK50679/)
  (feedback form, corrected models, model records, selenoproteins).
- The NCBI Eukaryotic Genome Annotation Pipeline,
  [process description](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process) (evidence, re-annotation).
- NCBI [Contact Us](https://www.ncbi.nlm.nih.gov/home/about/contact/) (Gene and RefSeq Feedback form).
- `results/tables/gold_standard.tsv` (columns `annotation_status`, `annotated_proteins`, `annotation_fix`, `notes`).
