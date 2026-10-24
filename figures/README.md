# Figures

Figure scripts read `results/tables/`, `results/trees/` and, where noted, the pipeline outputs in `work/`. Run them with
`make` (targets below); they need Infernal, matplotlib and Biopython, available in the image.

| Script | Make target | Figure |
|---|---|---|
| `stage1_secis_scores.py` | `secis-stage1` | RF00031 scores of the Selenoprofiles 3′ sequences and of random genomic windows, per genome (`stage1_secis_scores.*`; table `results/tables/stage1_secis_scores.tsv`) |
| `stage1_secis_scores.py` with the stage 2 model | `secis-stage2` | the same scores with the appendicularian SECIS model, with the curated *O. dioica* Sec genes and new random windows (`stage2_secis_scores.*`) |
| `stage1_overview.py` | `stage1-figures` | stage 1 per genome (18 genomes): selenocysteine predictions and their RF00031 support, fraction in correctly annotated genes, tRNA-Sec loci and machinery genes (`stage1_overview.*`) |
| `stage1_secis_scores.py --all` | `stage1-figures` | RF00031 scores of the 3′ sequences of all eighteen genomes, with random windows where the genome is available (`stage1_secis_scores_all.*`) |
| `phylogeny_collapsed.py` | `machinery-phylogeny` | display versions of the machinery trees: nodes with ultrafast bootstrap below 50 collapsed, clades of one taxonomic group (`config/taxon_groups.tsv`) drawn as one tip, tunicate sequences kept (`phylogeny_<family>_collapsed.*`) |
| `stage2/machinery_phylogeny.py` | `machinery-phylogeny` | full maximum-likelihood trees of PSTK, SEPHS and the SECIS-binding proteins, with protein lengths (`phylogeny_<family>.*`); the C20 versions in `results/trees/C20/` |
| `phylogeny_collapsed.py --trees results/trees/lba --tree-suffix .pmsf` | `machinery-lba-figures` | display versions of the LG+C60+F+G (PMSF) trees (`phylogeny_<family>_collapsed_pmsf.*`), proposed as the main tree figures |
| `phylogeny_collapsed.py --trees results/trees/domains --formats pdf` | `machinery-domain-figures` | display versions of the trees of each domain (`phylogeny_<family>_collapsed_domain_<domain>.pdf`; PDF only, to keep the repository small); with `DOMAIN_SUFFIX=.pmsf` (run by `machinery-domains-pmsf`), the PMSF trees (`..._pmsf.pdf`) |
| `stage2/machinery_phylogeny.py`, `phylogeny_collapsed.py` (selenoprotein families) | `selenoprotein-phylogeny`, `selenoprotein-pmsf` | trees of the selenoprotein families in `figures/selenoproteins/`: full (`phylogeny_<family>.pdf`), collapsed (`_collapsed.pdf`) and collapsed PMSF (`_collapsed_pmsf.pdf`); PDF only |
| `machinery_scheme.py` | `machinery-scheme` | schematic of the selenocysteine machinery, drawn from scratch with SBGN-style glyphs; gene names as in Table MACHINERY (`machinery_scheme.pdf`, and `machinery_scheme.svg` for editing, for example with Inkscape) |
| `reannotation_inputs.py`, `reannotation_figure.py` | (run directly; inputs need www.ebi.ac.uk and the genome sequence) | reannotation of one gene in the three *O. dioica* assemblies: public models against gold-standard models with SECIS elements, aligned on the selenocysteine TGA, introns coloured by splice site (A), gold-standard and reference proteins on the coordinates of a reference with its Pfam domains (B); `reannotation_<gene>.*` (SEPHS2, DIO_PAR, drawn with the previous layout); inputs with the source of each value in `figures/inputs/reannotation_<gene>.{tsv,fa}` |
| `reannotation_page.py --genes ...` | (run directly, after `reannotation_inputs.py` for each gene) | A4 portrait page of reannotated genes, three per row: gene models (gene structures half as thick as in the one-gene figure, public models named in captions) and proteins of each gene; assembly names and gold-standard protein labels once per row, key at the bottom of the page; three rows (nine genes) fit on a page (`reannotation_page.*`, or `--name`) |

For `stage1_secis_scores.py`, `--max` turns off the `cmsearch` filters, so that every sequence gets its best score
instead of "no hit" (slower), and `--replot` redraws a figure from its table. In the tree figures, bands mark the
outgroup family, the paralogs of the gene that names the family and the vertebrate clade of that gene; a column gives
the protein lengths (median and range for collapsed tips).

The schematic (`machinery_scheme.*`) is Figure SCHEME of the manuscript, in the introduction, with its legend.
