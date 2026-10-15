#!/usr/bin/env bash
# Stage 1 for one genome: Selenoprofiles, miniprot Sec codons, Rfam RF00031 SECIS, Secmarker tRNA-Sec,
# and (when an annotation exists) the annotation check. Each step is skipped when its output exists.
# usage: run_genome.sh <genome id>     environment: CPUS (default 1); SELENOPROFILES_CACHE (optional directory
#        with <id>.selenoprofiles.gtf from an earlier run of the same version and profiles, to avoid re-running it)
set -euo pipefail
id=$1
R=$(cd "$(dirname "$0")/.." && pwd)
B=$R/stage1/bin
W=$R/work/stage1/$id; mkdir -p "$W"; cd "$W"
G=$R/data/genomes/$id/genome.fna.gz
Q=$R/work/stage1/queries.faa
DATA=$R/data/references/selenoprofiles_data
CPUS=${CPUS:-1}
CHUNK=120000000      # miniprot index of at most 120 Mb per run (about 3 GB of RAM)

# 1. Selenoprofiles 4.5.8, profiles eukarya,machinery
if [[ ! -s $id.selenoprofiles.gtf ]]; then
    if [[ -n ${SELENOPROFILES_CACHE:-} && -s $SELENOPROFILES_CACHE/$id.selenoprofiles.gtf ]]; then
        cp "$SELENOPROFILES_CACHE/$id.selenoprofiles.gtf" .
    else
        export HOME=$W; selenoprofiles -setup > /dev/null
        gunzip -c "$G" > target.fa; mkdir -p tmp
        selenoprofiles -o sp_results -t target.fa -s "$id" -p eukarya,machinery -selenoprofiles_data_dir "$DATA" \
            -temp tmp -ncpus "$CPUS" -output_gtf_file "$id.selenoprofiles.gtf" -output_fasta_file "$id.selenoprofiles.faa"
        rm -rf tmp sp_results target.fa target.fa.n?? formatdb.log
    fi
fi

# 2. miniprot alignment of the Sec-containing profile sequences, by genome chunks; Sec codons with sec_check.py
if [[ ! -s $id.sec_check.tsv ]]; then
    seqkit fx2tab -n -l -i "$G" | awk -v M=$CHUNK '{s+=$2; if (s>M) {k++; s=$2}; print $1"\t"k+0}' > chunks.tsv
    : > miniprot.gff.part
    for k in $(cut -f2 chunks.tsv | sort -un); do
        awk -v k="$k" '$2==k{print $1}' chunks.tsv > ids.txt
        seqkit grep -f ids.txt "$G" > chunk.fa 2> /dev/null
        miniprot -t "$CPUS" -I --gff chunk.fa "$Q" 2> /dev/null >> miniprot.gff.part
    done
    rm -f chunk.fa ids.txt; mv miniprot.gff.part miniprot.gff
    python3 "$B/sec_check.py" --gff miniprot.gff --queries "$Q" --genome "$G" --output "$id.sec_check.tsv"
fi

# 3. Sec sites (Selenoprofiles + miniprot) and Rfam RF00031 on the 3 kb downstream of each Sec codon
if [[ ! -s $id.secis_rf00031.tbl ]]; then
    python3 "$B/secis_windows.py" --genome "$G" --gtf "$id.selenoprofiles.gtf" --sec-check "$id.sec_check.tsv" --prefix "$id"
    cmsearch --toponly --nohmmonly -T 10 --cpu "$CPUS" --tblout "$id.secis_rf00031.tbl.part" "$R/data/references/RF00031.cm" "$id.windows.fa" > /dev/null
    mv "$id.secis_rf00031.tbl.part" "$id.secis_rf00031.tbl"
fi

# 4. tRNA-Sec with Secmarker 1.0.0
if [[ ! -d secmarker ]]; then
    gunzip -c "$G" > target.fa
    python -m secmarker.secmarker_lib -t target.fa -o secmarker.part -cpu "$CPUS" > secmarker.log 2>&1 || true
    mv secmarker.part secmarker; rm -f target.fa secmarker/genome_link.fa*
fi

# 5. Annotation check (Tico et al. 2026 classes), when the genome has an annotation
A=$(ls "$R/data/genomes/$id"/annotation.gff3* 2>/dev/null | head -1 || true)
if [[ -n $A && ! -s $id.sec_annotation.tsv ]]; then
    python3 "$B/sec_annotation_check.py" --sample "$id" --selenoprofiles-gtf "$id.selenoprofiles.gtf" --annotation "$A" --prefix "$id"
fi
# 6. Transcript evidence (genomes with rows in config/transcripts.tsv): spliced alignment, Sec codon in an exon
labels=(); alns=()
while IFS=$'\t' read -r tid tgenome rest; do
    [[ $tgenome == "$id" && -s $R/data/transcripts/$tid.fa.gz ]] || continue
    [[ -s $tid.alignments.tsv ]] || { python3 "$B/map_transcripts.py" --genome "$G" --transcripts "$R/data/transcripts/$tid.fa.gz" --output "$tid.alignments.tsv.part" && mv "$tid.alignments.tsv.part" "$tid.alignments.tsv"; }
    labels+=("$tid"); alns+=("$tid.alignments.tsv")
done < <(tail -n +2 "$R/config/transcripts.tsv")
if (( ${#labels[@]} )) && [[ ! -s $id.transcript_support.tsv ]]; then
    python3 "$B/transcript_support.py" --sites "$id.sec_sites.tsv" --alignments "${alns[@]}" --labels "${labels[@]}" --output "$id.transcript_support.tsv"
fi
echo "stage 1 done: $id"
