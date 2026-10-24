#!/bin/sh
# Inputs of the reannotation figures (figures/reannotation_figure.py, figures/reannotation_page.py), with the
# curator's notes. Needs www.ebi.ac.uk and the OKI2018 and OdB3 genome sequences (FASTA, contig names as in the gold
# standard: NC_* for OKI2018, FN* for OdB3):
#   OKI_GENOME=... ODB3_GENOME=... sh figures/reannotation_inputs.sh
# GOLD (default curation/gold_standard) and OUT (default figures/inputs) can be set too.
set -e
GOLD=${GOLD:-curation/gold_standard}
OUT=${OUT:-figures/inputs}
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --genome OKI="$OKI_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" "$@" > /dev/null; }

run --gene SEPHS2 --embl FN653032.1 --reference "human SEPHS2=Q99611" \
    --reference-fasta "Ciona intestinalis SEPHS=results/models/machinery_profiles/SEPHS2.fa:198432475:XP_002129526.1" \
    --note "OKI:XP_081635792.1=coding sequence starts 79 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S25.g14713.01.t=TGA exon skipped; fused with the upstream gene" \
           "OdB3:CBY08599.1=other frame" "OdB3:CBY08598.1=starts at an internal ATG"

run --gene DIO_PAR --embl FN653278.1 --reference "human DIO2=Q92813" "Halocynthia roretzi deiodinase=Q6U6H1" \
    --note "OKI:XP_081636384.1=selenocysteine declared; one intron 11 nt away from the GA-AG intron of the transcripts" \
           "OSKA:Oidioi.CG.OSKA2016.S2.g03083.01.t=TGA in an intron of the annotated coding sequence" \
           "OdB3:CBY14356.1=ends before the TGA" "OdB3:CBY14355.1=other frame" "OdB3:CBY14354.1=starts after the TGA" \
    --gold-note "OdB3=protein from the ESTs (one extra base in the assembly)"

run --gene TXNRD --embl FN653072.1 --reference "human TXNRD1=Q16881" "Ciona intestinalis TXNRD=H2Y2A4" \
    --reference-note "Ciona intestinalis TXNRD=model ends at the TGA" \
    --note "OKI:XP_081643508.1=N-terminus from other exons (no glutaredoxin domain); TGA read as stop" \
           "OSKA:Oidioi.CG.OSKA2016.S8.g08759.01.t=extra exon in the GA-AG intron; ends in an intron before the TGA" \
           "OdB3:CBY19735.1=one intron differs from the ESTs; TGA read as stop" \
    --gold-note "OdB3=exons beyond the ESTs from the GenBank model"

run --gene SELENOT --embl FN653020.1 --reference "human SELENOT=P62341" \
    --reference-fasta "Ciona intestinalis SELENOT=results/trees/selenoproteins/SELENOT.faa:Cint_SELENOT.2.selenocysteine:Cint_SELENOT.2 (stage 1 prediction)" \
    --note "OKI:XP_081643256.1=same exon chain; coding sequence starts at an internal ATG, 62 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S10.g09589.01.t=TGA in an intron of the coding sequence; extra 5′ exon; ends in the GA-AG intron" \
           "OdB3:CBY23255.1=same exon chain; coding sequence starts at an internal ATG, 62 codons after the TGA"

run --gene DIO_chr1 --embl FN653016.1 --reference "human DIO3=P55073" "Halocynthia roretzi deiodinase=Q6U6H1" \
    --note "OKI:XP_081632420.1=coding sequence starts 4 codons after the TGA; model corrected for the assembly (low-quality protein)" \
           "OSKA:Oidioi.CG.OSKA2016.S37.g16254.01.t=TGA exon skipped; joined to 5 exons upstream of the gene" \
           "OdB3:CBY06720.1=ends at the TGA (read as stop)" "OdB3:CBY06719.1=starts 4 codons after the TGA" \
    --gold-note "OKI=mapped from the transcript (the assembly lacks 7 bases of the coding sequence and has 6 others)"

run --gene SELENOW_PAR --embl FN653024.1 --reference "human SELENOW=P63302" \
    --reference-fasta "Ciona intestinalis SELENOW=results/trees/selenoproteins/SELENOW.faa:Cint_SELENOW.1.selenocysteine:Cint_SELENOW.1 (stage 1 prediction)" \
    --note "OKI:XP_081636423.1=mRNA of another protein (592 aa), with the gene in its 3′ UTR" \
           "OdB3:CBY23864.1=mRNA of another protein (589 aa), with the gene in its 3′ UTR" \
    --gold-note "OSKA=mapped from the transcript (the assembly lacks 2 bases of the coding sequence)"
