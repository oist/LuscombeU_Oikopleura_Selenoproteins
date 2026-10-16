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
# Rows of a TSV file without its header. Empty fields become "-" because read(1) merges consecutive tabs;
# blank() turns them back into empty strings.
rows () { tail -n +2 "$1" | awk -F'\t' -v OFS='\t' '{for (i = 1; i <= NF; i++) if ($i == "") $i = "-"; print}'; }
blank () { local v; for v in "$@"; do [[ ${!v} == - ]] && printf -v "$v" ''; done; return 0; }
# Download url to dest (gzip-compressed, whatever the source compression) unless dest exists; then check the md5 of
# the uncompressed content when given.
decompress () { case $1 in *.bz2) if command -v lbzip2 > /dev/null; then lbzip2 -dc; else bzip2 -dc; fi ;; *.gz) gzip -dc ;; *) cat ;; esac; }
get_file () {  # url dest md5_uncompressed
    local url=$1 dest=$2 md5=${3:-}
    if [[ ! -s $dest ]]; then
        mkdir -p "$(dirname "$dest")"
        if ! { curl -fL --retry 3 -o "$dest.download" "$url" && decompress "$url" < "$dest.download" | gzip -n > "$dest.part"; }; then
            rm -f "$dest.download" "$dest.part"; echo "MISSING $dest: download failed ($url)"; fail=1; return
        fi
        rm -f "$dest.download"; mv "$dest.part" "$dest"
    fi
    if [[ -n $md5 ]] && ! gzip -dcf "$dest" | md5sum | grep -q "^$md5 "; then
        echo "MD5 mismatch: $dest (expected $md5 for the uncompressed content)"; fail=1
    fi
}
fetch_url () {  # url destination md5
    local url=$1 dest=$2 md5=${3:-}
    if [[ $url == PLACEHOLDER* ]]; then echo "SKIP $dest: source URL not yet filled in config ($url)"; fail=1; return; fi
    if [[ ! -s $dest ]]; then curl -fL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"; fi
    if [[ -n $md5 ]]; then echo "$md5  $dest" | md5sum -c --quiet - || { echo "MD5 mismatch: $dest"; fail=1; }; fi
}
# files distributed inside archives (config/archives.tsv). A copy placed by hand in data/downloads/ is used when
# the download fails; the archive and the extracted file are checked against their md5 sums.
mkdir -p data/downloads
while IFS=$'\t' read -r archive url md5 member dest member_md5 notes; do
    blank archive url md5 member dest member_md5 notes
    if [[ -s $dest ]]; then              # extracted earlier, or placed by hand: check the uncompressed content
        if [[ -n $member_md5 ]] && ! gzip -dcf "$dest" | md5sum | grep -q "^$member_md5 "; then
            echo "MD5 mismatch: $dest (expected $member_md5 for the uncompressed content)"; fail=1
        fi
        continue
    fi
    z=data/downloads/$archive
    if [[ ! -s $z ]]; then
        if ! { curl -fL --retry 3 -o "$z.part" "$url" && mv "$z.part" "$z"; }; then
            rm -f "$z.part"; echo "MISSING $dest: download of $archive failed; place a copy of the archive in $z, or the file itself in $dest"; fail=1; continue
        fi
    fi
    if [[ -n $md5 ]] && ! echo "$md5  $z" | md5sum -c --quiet -; then
        echo "MD5 mismatch: $z (an error page instead of the archive?)"; fail=1; continue
    fi
    [[ $member == PLACEHOLDER* ]] && { echo "SKIP $dest: archive member not yet filled in config/archives.tsv"; fail=1; continue; }
    mkdir -p "$(dirname "$dest")"
    unzip -p "$z" "$member" > "$dest.tmp"
    if [[ -n $member_md5 ]] && ! echo "$member_md5  $dest.tmp" | md5sum -c --quiet -; then
        echo "MD5 mismatch: $member in $z"; rm -f "$dest.tmp"; fail=1; continue
    fi
    if [[ $dest == *.gz ]]; then gzip -n -c "$dest.tmp" > "$dest"; rm "$dest.tmp"; else mv "$dest.tmp" "$dest"; fi
done < <(rows config/archives.tsv)

# genomes (and annotations)
rows config/genomes.tsv | while IFS=$'\t' read -r id species assembly source acc annot md5 role notes; do
    blank id species assembly source acc annot md5 role notes
    out=data/genomes/$id
    mkdir -p "$out"
    [[ $source == archive ]] && continue      # fetched from config/archives.tsv above
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
# transcripts: downloaded, or placed by hand when there is no public source (source "manual"); in both cases the
# file is checked by the md5 of its uncompressed content
rows config/transcripts.tsv | while IFS=$'\t' read -r id genome desc source url md5 notes; do
    blank id genome desc source url md5 notes
    dest=data/transcripts/$id.fa.gz
    if [[ ! -s $dest && ( $source == manual || $url == PLACEHOLDER* ) ]]; then
        echo "MISSING $dest: no download source yet; place a copy there (plain or gzipped)"; fail=1; continue
    fi
    get_file "$url" "$dest" "$md5"
done

# other annotation files (config/extra_files.tsv)
rows config/extra_files.tsv | while IFS=$'\t' read -r id genome desc url md5 dest notes; do
    blank id genome desc url md5 dest notes
    get_file "$url" "$dest" "$md5"
done
# references (Rfam model and seed; Selenoprofiles profiles)
fetch_url https://rfam.org/family/RF00031/cm data/references/RF00031.cm ""
fetch_url https://rfam.org/family/RF00031/alignment data/references/RF00031.sto ""
fetch_url https://rfam.org/family/RF01852/cm data/references/RF01852.cm ""
# Selenoprofiles profiles at the commit used by the pipeline's test profile
if [[ ! -d data/references/selenoprofiles_data ]]; then
    mkdir -p data/references/selenoprofiles_data
    curl -fsL https://codeload.github.com/mariottigenomicslab/selenoprotein_profiles/tar.gz/4829bc969570587c222f1ce4721038fd75c804a0 \
        | tar -xz --strip-components=1 -C data/references/selenoprofiles_data
fi
md5sum data/genomes/*/genome.fna.gz data/transcripts/*.fa.gz 2>/dev/null > data/MD5SUMS || true
echo "Done. Checksums of downloaded files: data/MD5SUMS"
exit $fail
