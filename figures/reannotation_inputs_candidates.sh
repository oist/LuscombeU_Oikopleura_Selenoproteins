#!/bin/sh
# Inputs of the reannotation plots of the candidate genes (step 5: TXNL1, KKSU), from their gold-standard models
# (make gold-candidates: curation/gold_candidates, results/tables/gold_candidates.tsv). Same needs as
# figures/reannotation_inputs.sh (www.ebi.ac.uk; OKI2018 and OdB3 contigs with the names of the gold standard):
#   OKI_GENOME=... ODB3_GENOME=... sh figures/reannotation_inputs_candidates.sh
# KKSU has no homolog outside the Oikopleuridae: its reference is the O. albicans open stretch of the comparative screen
# (curation/candidates.faa), and the page draws its proteins on the OKI2018 gold-standard protein.
set -e
GOLD=${GOLD:-curation/gold_candidates}
OUT=${OUT:-figures/inputs}
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --table results/tables/gold_candidates.tsv --genome OKI="$OKI_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" "$@" > /dev/null; }

run --gene TXNL1 --embl FN653019.1 --reference "human=O43396" \
    --reference-fasta "C. intestinalis=figures/inputs/reference_TXNL1_Cint.fa:XP_078490413.1:XP_078490413.1" \
    --note "OKI:XP_081630352.1=same exon chain; TGA read as stop" \
           "OSKA:Oidioi.CG.OSKA2016.S27.g15032.01.t=same coding exons; TGA read as stop" \
           "OdB3:CBY22487.1=same exon chain; TGA read as stop" \
    --gold-note "OSKA=5′ UTR exon from the transcripts" \
    --reference-note "human=the O. dioica protein aligns with the PITH domain only (human residues 120-283)"

run --gene KKSU --embl FN653017.1 \
    --reference-fasta "O. albicans=curation/candidates.faa:KKSU_Oalb:SCLG01000111.1 (genome; open stretch of the comparative screen)" \
    --note "OKI:XP_081632581.1=mRNA of another protein (173 aa) with this gene in its 3′ UTR" \
           "OSKA:Oidioi.CG.OSKA2016.S22.g13919.01.t=same exon chain; TGA read as stop" \
           "OdB3:CBY21026.1=same exon chain; TGA read as stop"
