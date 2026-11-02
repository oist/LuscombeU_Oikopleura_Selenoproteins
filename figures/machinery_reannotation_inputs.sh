#!/bin/sh
# Inputs of the machinery reannotation panels (figures/reannotation_page.py --origin start --no-secis), from the
# machinery gold standard (curation/gold_standard_machinery). Needs www.ebi.ac.uk and the OKI2018 and OdB3 genome
# sequences (FASTA, contig names as in the gold standard: NC_* for OKI2018, FN* for OdB3):
#   OKI_GENOME=... ODB3_GENOME=... sh figures/machinery_reannotation_inputs.sh
# Reference proteins: the human protein (UniProt) and the Ciona intestinalis RefSeq protein of the machinery tree.
set -e
GOLD=${GOLD:-curation/gold_standard_machinery}
OUT=${OUT:-figures/inputs}
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --table results/tables/gold_standard_machinery.tsv \
          --genome OKI="$OKI_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" "$@" > /dev/null; }

run --gene PSTK --embl FN653020.1 --reference "human=Q8IV42" \
    --reference-fasta "C. intestinalis=results/trees/pstk.faa:Cint_RefSeq_PSTK_XP_002120985.2:XP_002120985.2"
