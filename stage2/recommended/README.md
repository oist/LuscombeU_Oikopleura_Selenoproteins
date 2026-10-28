# Selenoprotein annotation of appendicularian genome assemblies: recommended settings

Settings and data files for [oist/LuscombeU_Selenoproteins_Pipeline](https://github.com/oist/LuscombeU_Selenoproteins_Pipeline)
release 0.2.0, chosen in stage 2 of this repository (`stage2/README.md`, `results/tables/stage2_benchmark_*.tsv`).
The deposit archive (`make stage2-bundle`) holds this file and:

| File | Use |
|---|---|
| `appendicularia_selenoprofiles.tar.gz` | `--selenoprofiles_overrides`: Selenoprofiles profiles of the selenoproteins and of the machinery, extended with appendicularian proteins (standard names, replacing the standard profiles), and two new profiles, `KKSU` and `TXNL1b` |
| `oikopleura_secis.cm` | `--secis_cm`: SECIS covariance model trained on *Oikopleura* genes, with a gathering threshold of 20 bits (used by the pipeline's `--cut_ga`) |
| `appendicularia.config` | Nextflow configuration (`-c`): Selenoprofiles writes 2,000 nt of 3′ sequence after each prediction, the window length on which the 20-bit threshold was set |
| `params.yaml` | parameters file (`-params-file`) to edit |
| `samplesheet.csv` | example samplesheet: the profiles column must name the new profiles, in quotes |
| `expectation_oikopleura_dioica.csv` | expected selenoproteome of *O. dioica* (20 genes), for `--expectation` |
| `profile_members.tsv` | the proteins added to each profile |

## Run

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r 0.2.0 -profile <profile> \
    -params-file params.yaml -c appendicularia.config
```

with, in the samplesheet, `"eukarya,machinery,KKSU,TXNL1b"` in the `profiles` column of each appendicularian genome.
For *O. dioica*, uncomment the expectation lines of `params.yaml` (the check counts selenocysteine sites per family;
TXNRD, SELENOT and SELENOK are on the X-specific region).

## What to expect

Measured on *O. dioica* (results of 2026-10-29): with profiles built without the gold standard of the species tested
(leave-one-species-out), the pipeline found the selenocysteine codon of every reference gene: 20 of 20 in OKI2018
and OSKA2016, 19 of 19 in OdB3, and 19, 18 and 19 in the new assemblies of the same three species KUM-M3-7f,
AOM-5-5f and Bar2_p4 (the genes whose codon could be placed on them from the gold standard), against 11 to 14 with
the default settings. Each codon found by Selenoprofiles had a SECIS hit of
at least 20 bits in the 2,000 nt after the gene. Sites outside the reference genes were other copies of selenoprotein
genes (haplotypes, duplicates) or genes the references missed, and miniprot calls 3 nt from a true codon.

Read the selenocysteine sites in `compare/<sample>.sec_sites.tsv` (or `compare/<sample>.sec_codons.gff3`), from both
tools: in *O. dioica*, SELENOT and SEPHS2 are found by Selenoprofiles only (codons 1 and 2 nt from an intron), and
some genes by miniprot only. Selenoprofiles gene models: median coding-base sensitivity and precision 1.00, introns
identical to the reference in about half of the genes; use them as evidence for gene annotation rather than as final
models.

Changes to the standard profiles: the SELENOS filter no longer requires a motif conserved in vertebrates for
predictions labelled selenocysteine (the *O. dioica* SELENOS lacks it); the GPX profile keeps its Gene Ontology terms.
Limits: fragments shorter than 61 residues covering 40% of the profile or less are removed by Selenoprofiles' default
filter (two in the *O. albicans* assembly); the tests outside *O. dioica* rest on stage 1 predictions.
