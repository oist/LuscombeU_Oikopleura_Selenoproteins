#!/usr/bin/env bash
#SBATCH --job-name=fetch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.log
# Download the public data listed in config/*.tsv into data/ and check md5 sums when given.
# Requirements: NCBI datasets command-line tool (for "ncbi" rows), curl, unzip, tar, gzip, md5sum (all in the
# image built from environment.yml; see "Before running" in README.md).
# Run from the repository root: scripts/fetch.sh, or sbatch [--partition=<data transfer partition>] scripts/fetch.sh
set -euo pipefail
# Repository root: under Slurm the script runs from a spool copy, so use the submission directory.
if [[ -n ${SLURM_JOB_ID:-} ]]; then R=$SLURM_SUBMIT_DIR; else R=$(cd "$(dirname "$0")/.." && pwd); fi
[[ -f $R/config/genomes.tsv ]] || { echo "Submit or run from the repository root ($R is not)" >&2; exit 1; }
cd "$R"
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
# Selenoprofiles profiles at the commit used by the pipeline's test profile
if [[ ! -d data/references/selenoprofiles_data ]]; then
    mkdir -p data/references/selenoprofiles_data
    curl -fsL https://codeload.github.com/mariottigenomicslab/selenoprotein_profiles/tar.gz/4829bc969570587c222f1ce4721038fd75c804a0 \
        | tar -xz --strip-components=1 -C data/references/selenoprofiles_data
fi
md5sum data/genomes/*/genome.fna.gz data/transcripts/*.fa.gz 2>/dev/null > data/MD5SUMS || true
echo "Done. Checksums of downloaded files: data/MD5SUMS"
exit $fail
