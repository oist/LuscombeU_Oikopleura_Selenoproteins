#!/bin/sh
# Inputs of the machinery reannotation panels (figures/reannotation_page.py --origin start --no-secis), from the
# machinery gold standard (curation/gold_standard_machinery). Needs www.ebi.ac.uk and the OKI2018 and OdB3 genome
# sequences (FASTA, contig names as in the gold standard: NC_* for OKI2018, FN* for OdB3):
#   OKI_GENOME=... ODB3_GENOME=... sh figures/machinery_reannotation_inputs.sh
# Reference proteins: the human protein (UniProt) and the Ciona intestinalis RefSeq protein of the machinery tree.
set -e
GOLD=${GOLD:-curation/gold_standard_machinery}
OUT=${OUT:-figures/inputs}
TRANSCRIPTS=${TRANSCRIPTS:-4}   # supporting transcripts drawn per assembly
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --table results/tables/gold_standard_machinery.tsv \
          --genome OKI="$OKI_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" --transcripts "$TRANSCRIPTS" "$@" > /dev/null; }

run --gene PSTK --embl FN653020.1 --reference "human=Q8IV42" \
    --reference-fasta "C. intestinalis=results/trees/pstk.faa:Cint_RefSeq_PSTK_XP_002120985.2:XP_002120985.2" \
    --note "OKI:XP_081642556.1=first three exons only; reads on into the third intron (110 codons, 98 of them as in the gold standard)" \
           "OSKA:Oidioi.CG.OSKA2016.S10.g09825.01.t=starts in the fourth exon and continues past the stop codon into seven exons downstream (604 codons)" \
           "OdB3:CBY23135.1=same first three exons; other first site of the fourth exon, and no last exon (178 codons)" \
    --gold-note "OKI=two splice sites outside the intron classes, as in both supporting transcripts (possible assembly errors)" \
                "OSKA=two splice sites outside the intron classes, as in the supporting transcript" \
                "OdB3=exon chain set by the curator: the GenBank exons 1 to 3 and the last intron of the ESTs"
