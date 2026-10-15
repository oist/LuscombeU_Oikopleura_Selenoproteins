#!/usr/bin/env bash
# Download the public data listed in config/*.tsv into data/ and check md5 sums when given.
# Requirements: NCBI datasets CLI (for "ncbi" rows), curl, unzip, md5sum.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/genomes data/transcripts data/references
fail=0
fetch_url () {  # url destination md5
    local url=$1 dest=$2 md5=${3:-}
    if [[ $url == PLACEHOLDER* ]]; then echo "SKIP $dest: source URL not yet filled in config ($url)"; fail=1; return; fi
    if [[ ! -s $dest ]]; then curl -fL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"; fi
    if [[ -n $md5 ]]; then echo "$md5  $dest" | md5sum -c --quiet - || { echo "MD5 mismatch: $dest"; fail=1; }; fi
}
# genomes (and annotations)
tail -n +2 config/genomes.tsv | while IFS=$'\t' read -r id species assembly source acc annot md5 role notes; do
    out=data/genomes/$id
    mkdir -p "$out"
    if [[ $source == ncbi ]]; then
        if [[ ! -s $out/genome.fna.gz ]]; then
            datasets download genome accession "$acc" --include genome,gff3 --filename "$out/ncbi.zip"
            unzip -o -q "$out/ncbi.zip" -d "$out/ncbi"
            find "$out/ncbi" -name '*_genomic.fna' -exec cat {} + | gzip > "$out/genome.fna.gz"
            gff=$(find "$out/ncbi" -name 'genomic.gff' | head -1 || true)
            [[ -n $gff ]] && gzip -c "$gff" > "$out/annotation.gff3.gz"
            rm -rf "$out/ncbi" "$out/ncbi.zip"
        fi
    else
        fetch_url "$acc" "$out/genome.fna.gz" "$md5"
        url=${annot%% *}; [[ -n $url ]] && fetch_url "$url" "$out/annotation.gff3" ""
    fi
done
# transcripts
tail -n +2 config/transcripts.tsv | while IFS=$'\t' read -r id genome desc source url md5 notes; do
    fetch_url "$url" "data/transcripts/$id.fa.gz" "$md5"
done
# references (Rfam model and seed; Selenoprofiles profiles)
fetch_url https://rfam.org/family/RF00031/cm data/references/RF00031.cm ""
fetch_url https://rfam.org/family/RF00031/alignment data/references/RF00031.sto ""
[[ -d data/references/selenoprotein_profiles ]] || git clone -q --depth 1 --branch v1.3.0 https://github.com/marco-mariotti/selenoprotein_profiles data/references/selenoprotein_profiles
md5sum data/genomes/*/genome.fna.gz data/transcripts/*.fa.gz 2>/dev/null > data/MD5SUMS || true
echo "Done. Checksums of downloaded files: data/MD5SUMS"
exit $fail
