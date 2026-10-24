#!/bin/sh
# Inputs of the reannotation figures (figures/reannotation_figure.py, figures/reannotation_page.py), with the
# curator's notes. Needs www.ebi.ac.uk and the OKI2018 and OdB3 genome sequences (FASTA, contig names as in the gold
# standard: NC_* for OKI2018, FN* for OdB3):
#   OKI_GENOME=... ODB3_GENOME=... sh figures/reannotation_inputs.sh
# GOLD (default curation/gold_standard) and OUT (default figures/inputs) can be set too.
# Reference proteins: the human protein (UniProt; best reciprocal hit of stage 2) and the Ciona intestinalis protein
# of the stage 1 predictions on the Ciona genome (results/trees/; the prediction with the best local alignment score
# to the gold-standard protein, or one starting with a methionine when the scores are close), labelled by species so
# that the page draws one row per species.
set -e
GOLD=${GOLD:-curation/gold_standard}
OUT=${OUT:-figures/inputs}
T=results/trees/selenoproteins
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --genome OKI="$OKI_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" "$@" > /dev/null; }

run --gene SEPHS2 --embl FN653032.1 --reference "human=Q99611" \
    --reference-fasta "C. intestinalis=results/trees/SEPHS2.faa:Cint_SEPHS2.5.glycine:Cint_SEPHS2.5" \
    --note "OKI:XP_081635792.1=coding sequence starts 79 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S25.g14713.01.t=TGA exon skipped; fused with the upstream gene" \
           "OdB3:CBY08599.1=other frame" "OdB3:CBY08598.1=starts at an internal ATG"

run --gene DIO_PAR --embl FN653278.1 --reference "human=Q92813" \
    --reference-fasta "C. intestinalis=$T/DIO.faa:Cint_DIO.2.selenocysteine:Cint_DIO.2" \
    --note "OKI:XP_081636384.1=selenocysteine declared; one intron 11 nt away from the GA-AG intron of the transcripts" \
           "OSKA:Oidioi.CG.OSKA2016.S2.g03083.01.t=TGA in an intron of the annotated coding sequence" \
           "OdB3:CBY14356.1=ends before the TGA" "OdB3:CBY14355.1=other frame" "OdB3:CBY14354.1=starts after the TGA" \
    --gold-note "OdB3=protein from the ESTs (one extra base in the assembly)"

run --gene TXNRD --embl FN653072.1 --reference "human=Q16881" \
    --reference-fasta "C. intestinalis=$T/TXNRD.faa:Cint_TXNRD.14.selenocysteine:Cint_TXNRD.14" \
    --note "OKI:XP_081643508.1=N-terminus from other exons (no glutaredoxin domain); TGA read as stop" \
           "OSKA:Oidioi.CG.OSKA2016.S8.g08759.01.t=extra exon in the GA-AG intron; ends in an intron before the TGA" \
           "OdB3:CBY19735.1=one intron differs from the ESTs; TGA read as stop" \
    --gold-note "OdB3=exons beyond the ESTs from the GenBank model"

run --gene SELENOT --embl FN653020.1 --reference "human=P62341" \
    --reference-fasta "C. intestinalis=$T/SELENOT.faa:Cint_SELENOT.2.selenocysteine:Cint_SELENOT.2" \
    --note "OKI:XP_081643256.1=same exon chain; coding sequence starts at an internal ATG, 62 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S10.g09589.01.t=TGA in an intron of the coding sequence; extra 5′ exon; ends in the GA-AG intron" \
           "OdB3:CBY23255.1=same exon chain; coding sequence starts at an internal ATG, 62 codons after the TGA"

run --gene DIO_chr1 --embl FN653016.1 --reference "human=P55073" \
    --reference-fasta "C. intestinalis=$T/DIO.faa:Cint_DIO.2.selenocysteine:Cint_DIO.2" \
    --note "OKI:XP_081632420.1=coding sequence starts 4 codons after the TGA; model corrected for the assembly (low-quality protein)" \
           "OSKA:Oidioi.CG.OSKA2016.S37.g16254.01.t=TGA exon skipped; joined to 5 exons upstream of the gene" \
           "OdB3:CBY06720.1=ends at the TGA (read as stop)" "OdB3:CBY06719.1=starts 4 codons after the TGA" \
    --gold-note "OKI=mapped from the transcript (the assembly lacks 7 bases of the coding sequence and has 6 others)"

run --gene SELENOW_PAR --embl FN653024.1 --reference "human=P63302" \
    --reference-fasta "C. intestinalis=$T/SELENOW.faa:Cint_SELENOW.1.selenocysteine:Cint_SELENOW.1" \
    --note "OKI:XP_081636423.1=mRNA of another protein (592 aa), with the gene in its 3′ UTR" \
           "OdB3:CBY23864.1=mRNA of another protein (589 aa), with the gene in its 3′ UTR" \
    --gold-note "OSKA=mapped from the transcript (the assembly lacks 2 bases of the coding sequence)"

run --gene GPX_PAR --embl FN653291.1 --reference "human=P36969" \
    --reference-fasta "C. intestinalis=$T/GPX.faa:Cint_GPX.5.selenocysteine:Cint_GPX.5" \
    --note "OKI:rna-XR_014134294.1=non-coding RNA that starts upstream and ends in the TGA exon" \
           "OKI:XP_081635498.1=selenocysteine declared; same model as the gold standard" \
           "OSKA:Oidioi.CG.OSKA2016.S25.g14813.01.t=TGA exon skipped; starts in the second exon" \
           "OdB3:CBY25160.1=starts 14 nt before the TGA, in another frame for 41 codons" \
    --gold-note "OdB3=mapped from the transcript (one extra base in the assembly)"

run --gene SELENOK --embl FN653034.1 --reference "human=Q9Y6D0" \
    --reference-fasta "C. intestinalis=$T/SELENOK.faa:Cint_SELENOK.2.selenocysteine:Cint_SELENOK.2" \
    --note "OdB3:CBY09056.1=same exon chain; TGA read as stop"

run --gene SELENOS --embl FN653033.1 --reference "human=Q9BQE4" \
    --reference-fasta "C. intestinalis=$T/SELENOS.faa:Cint_SELENOS.1.selenocysteine:Cint_SELENOS.1" \
    --note "OKI:XP_081639144.1=same exon chain; frame lost where the assembly lacks 13 bases; ends before the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S13.g10968.01.t=starts 9 codons upstream; TGA read as stop" \
           "OdB3:CBY08859.1=N-terminal part; ends in an intron" "OdB3:CBY08858.1=other first exon; TGA read as stop" \
    --gold-note "OKI=mapped from the transcript (the assembly lacks 13 bases of the coding sequence)"

# page 2
run --gene GPX_chr1a --embl FN653021.1 --reference "human=P36969" \
    --reference-fasta "C. intestinalis=$T/GPX.faa:Cint_GPX.5.selenocysteine:Cint_GPX.5" \
    --note "OKI:XP_081637767.1=same model; TGA translated as X (selenocysteine not declared)" \
           "OSKA:Oidioi.CG.OSKA2016.S6.g07306.01.t=other 5′ exons; TGA read as stop" \
           "OSKA:Oidioi.CG.OSKA2016.S6.g07307.01.t=starts 13 codons after the TGA" \
           "OdB3:CBY23346.1=TGA read as stop" "OdB3:CBY23347.1=starts 13 codons after the TGA" \
    --gold-note "OdB3=the two GenBank models joined through the TGA (no transcript)"

run --gene GPX_chr1b --embl FN653042.1 --reference "human=P36969" \
    --reference-fasta "C. intestinalis=$T/GPX.faa:Cint_GPX.5.selenocysteine:Cint_GPX.5" \
    --note "OKI:XP_081636615.1=same exon chain; coding sequence starts at an internal ATG, 17 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S3.g04759.01.t=starts 17 codons after the TGA; reads through one intron, places another differently" \
           "OdB3:CBY19528.1=TGA read as stop" "OdB3:CBY19527.1=starts 17 codons after the TGA"

run --gene GPX_chr2 --embl FN653019.1 --reference "human=P36969" \
    --reference-fasta "C. intestinalis=$T/GPX.faa:Cint_GPX.1.selenocysteine:Cint_GPX.1" \
    --note "OKI:XP_081630045.1=same exon chain; TGA translated as X (selenocysteine not declared)" \
           "OSKA:Oidioi.CG.OSKA2016.S27.g15046.01.t=starts in the second exon; TGA in an intron of the coding sequence" \
           "OdB3:CBY22474.1=TGA in an extra intron of the coding sequence"

run --gene SELENOF --embl FN653029.1 --reference "human=O60613" \
    --reference-fasta "C. intestinalis=$T/SELENOF.faa:Cint_SELENOF.1.selenocysteine:Cint_SELENOF.1" \
    --note "OKI:XP_081628534.1=selenocysteine declared; same model as the gold standard" \
           "OSKA:Oidioi.CG.OSKA2016.S16.g12238.01.t=TGA in an intron of the coding sequence" \
           "OdB3:CBY08124.1=same exon chain; TGA read as stop"

run --gene SELENOO --embl FN653064.1 --reference "human=Q9BVL4" \
    --reference-fasta "C. intestinalis=$T/SELENOO.faa:Cint_SELENOO.1.selenocysteine:Cint_SELENOO.1" \
    --note "OKI:XP_081634508.1=TGA read as stop; GT-AG or GA-AG introns where the gold standard, mapped from the transcript, has other sites" \
           "OSKA:Oidioi.CG.OSKA2016.S11.g10191.01.t=N-terminal part" "OSKA:Oidioi.CG.OSKA2016.S11.g10192.01.t=central part" \
           "OSKA:Oidioi.CG.OSKA2016.S11.g10193.01.t=C-terminal part, TGA read as stop, skips one exon" \
    --gold-note "OKI=mapped from the transcript (the assembly lacks 5 bases of the coding sequence)" \
                "OdB3=the exon of the TGA is not in the assembly"

run --gene SELENOW_chr1a --embl FN653018.1 --reference "human=P63302" \
    --reference-fasta "C. intestinalis=$T/SELENOW.faa:Cint_SELENOW.1.selenocysteine:Cint_SELENOW.1" \
    --note "OSKA:Oidioi.CG.OSKA2016.S7.g08360.01.t=TGA in an intron of the coding sequence"

run --gene SELENOW_chr2 --embl FN653019.1 --reference "human=Q9BRT3" --reference-note "human=MIEN1, best human hit, a cysteine homolog" \
    --reference-fasta "C. intestinalis=$T/SELENOW.faa:Cint_SELENOW.3.selenocysteine:Cint_SELENOW.3" \
    --note "OKI:XP_081633592.1=same exon chain; coding sequence starts at an internal ATG, 24 codons after the TGA" \
           "OSKA:Oidioi.CG.OSKA2016.S40.g16478.01.t=starts 24 codons after the TGA; one intron placed differently" \
           "OdB3:CBY22624.1=same exon chain; coding sequence starts at an internal ATG, 24 codons after the TGA" \
    --gold-note "OSKA=3′ UTR exons of a read-through transcript into the next gene left out"

run --gene SELENOW_pair1 --embl FN653018.1 --reference "human=P63302" \
    --reference-fasta "C. intestinalis=$T/SELENOW.faa:Cint_SELENOW.1.selenocysteine:Cint_SELENOW.1" \
    --note "OKI:XP_081628722.1=mRNA of another protein (674 aa), with SELENOW_pair2 and this gene in its 3′ UTR" \
           "OSKA:Oidioi.CG.OSKA2016.S5.g07175.01.t=one model across SELENOW_pair2 and this gene; both TGAs in introns" \
           "OdB3:CBY21827.1=TGA in an intron of the coding sequence"

run --gene SELENOW_pair2 --embl FN653018.1 --reference "human=P63302" \
    --reference-fasta "C. intestinalis=$T/SELENOW.faa:Cint_SELENOW.1.selenocysteine:Cint_SELENOW.1" \
    --note "OKI:XP_081628722.1=mRNA of another protein (674 aa), with this gene and SELENOW_pair1 in its 3′ UTR" \
           "OSKA:Oidioi.CG.OSKA2016.S5.g07175.01.t=one model across this gene and SELENOW_pair1; both TGAs in introns" \
           "OdB3:CBY21828.1=mRNA of another protein (94 aa) further downstream, with this gene in its 5′ UTR"
