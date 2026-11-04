#!/bin/sh
# Inputs of the machinery reannotation panels (figures/reannotation_page.py --origin start --no-secis), from the
# machinery gold standard (curation/gold_standard_machinery). Needs www.ebi.ac.uk and the OKI2018, OSKA2016 and OdB3
# genome sequences (FASTA, contig names as in the gold standard: NC_* for OKI2018, S* for OSKA2016, FN* for OdB3; the
# splice sites of the models are read from them, the curator's chains included):
#   OKI_GENOME=... OSKA_GENOME=... ODB3_GENOME=... sh figures/machinery_reannotation_inputs.sh
# Reference proteins: the human protein (UniProt) and the Ciona intestinalis RefSeq protein of the machinery tree.
set -e
GOLD=${GOLD:-curation/gold_standard_machinery}
OUT=${OUT:-figures/inputs}
TRANSCRIPTS=${TRANSCRIPTS:-4}   # supporting transcripts drawn per assembly
run() { python3 figures/reannotation_inputs.py --gold "$GOLD" --table results/tables/gold_standard_machinery.tsv \
          --genome OKI="$OKI_GENOME" OSKA="$OSKA_GENOME" OdB3="$ODB3_GENOME" --out "$OUT" --transcripts "$TRANSCRIPTS" "$@" > /dev/null; }

run --gene PSTK --embl FN653020.1 --reference "human=Q8IV42" \
    --reference-fasta "C. intestinalis=results/trees/pstk.faa:Cint_RefSeq_PSTK_XP_002120985.2:XP_002120985.2" \
    --note "OKI:XP_081642556.1=first three exons only; reads on into the third intron (110 codons, 98 of them as in the gold standard)" \
           "OSKA:Oidioi.CG.OSKA2016.S10.g09825.01.t=starts in the fourth exon and continues past the stop codon into seven exons downstream (604 codons)" \
           "OdB3:CBY23135.1=same first three exons; other first site of the fourth exon, and no last exon (178 codons)" \
    --gold-note "OKI=third and fourth introns with splice sites outside the intron classes (TG-AG, AG-AA), as in the three supporting transcripts" \
                "OSKA=third and fourth introns with splice sites outside the intron classes (AA-AG, TT-GT), as in the supporting transcripts; the third intron is retained in two of them" \
                "OdB3=exon chain set by the curator: the GenBank exons 1 to 3, the acceptor of the third intron (GA-AG; retained in both ESTs) from the prediction, and the last intron of the ESTs (AT-AG)"

# Genes under review (2026-11-03): public models, gold-standard models and supporting transcripts; notes to add after
# the review. Reference proteins: human (UniProt) and a Ciona protein (C. intestinalis from the machinery trees, or the
# C. robusta RefSeq protein of the Selenoprofiles profile).
run --gene SBP2 --embl FN653040.1 --reference "human=Q93073" \
    --reference-fasta "C. intestinalis=results/trees/SBP2.faa:Cint_RefSeq_SECISBP2L_XP_078482116.1:XP_078482116.1"
run --gene SEPHS_b --reference "human=P49903" \
    --reference-fasta "C. intestinalis=results/trees/SEPHS2.faa:Cint_SEPHS2.5.glycine:Selenoprofiles prediction, glycine form" \
    --gold-note "OSKA=no Selenoprofiles prediction; anchor from an alignment of the OKI2018 protein (curation/machinery_extra_anchors.tsv); exon chain of the transcript" \
                "OdB3=no Selenoprofiles prediction; anchor from an alignment of the OKI2018 protein; locus split by two assembly gaps (100 and 675 N): 34 residues read from the first gap (X), no EST; protein left out of the trees"
run --gene SEPSECS --embl FN653034.1 --reference "human=Q9HD40" \
    --reference-fasta "C. robusta=results/models/machinery_profiles/SecS.fa:XP_002120570:XP_002120570.1"
run --gene EEFSEC --embl FN653026.1 --reference "human=P57772" \
    --reference-fasta "C. robusta=results/models/machinery_profiles/eEFsec.fa:XP_002124067:XP_002124067.1"

# Panel C of Figure PSTK: equivalent placements of the third and fourth introns, in the gold-standard assemblies and
# in the assemblies of other individuals (results/gene_regions/PSTK, make pstk-regions; needs minimap2)
python3 figures/intron_placements.py --gene PSTK --assembly OKI=gold KUM=lift:OKI OSKA=gold AOM=lift:OSKA OdB3=gold \
    Bar2=lift:OdB3 --introns 3 4 --gold-dir "$GOLD" --out "$OUT/intron_placements_PSTK.tsv" > /dev/null
