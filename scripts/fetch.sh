#!/usr/bin/env bash
#SBATCH --job-name=selenoproteins-fetch
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
        if ! { curl -fsSL --retry 3 -o "$dest.download" "$url" && decompress "$url" < "$dest.download" | gzip -n > "$dest.part"; }; then
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
    if [[ ! -s $dest ]]; then
        if ! { curl -fsSL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"; }; then
            rm -f "$dest.part"; echo "MISSING $dest: download failed ($url)"; fail=1; return
        fi
    fi
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
        if ! { curl -fsSL --retry 3 -o "$z.part" "$url" && mv "$z.part" "$z"; }; then
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
while IFS=$'\t' read -r id species assembly source acc annot md5 role notes; do
    blank id species assembly source acc annot md5 role notes
    out=data/genomes/$id
    mkdir -p "$out"
    [[ $source == archive ]] && continue      # fetched from config/archives.tsv above
    if [[ $source == ncbi ]]; then
        if [[ ! -s $out/genome.fna.gz ]]; then     # written last: its presence means the download completed
            rm -rf "$out/ncbi" "$out/ncbi.zip"
            if ! { datasets download genome accession "$acc" --include genome,gff3 --filename "$out/ncbi.zip" \
                   && unzip -o -q "$out/ncbi.zip" -d "$out/ncbi"; }; then
                echo "MISSING $out/genome.fna.gz: download of $acc failed"; fail=1; rm -rf "$out/ncbi" "$out/ncbi.zip"; continue
            fi
            gff=$(find "$out/ncbi" -name 'genomic.gff' | head -1)
            if [[ -n $gff ]]; then gzip -n -c "$gff" > "$out/annotation.gff3.gz.part" && mv "$out/annotation.gff3.gz.part" "$out/annotation.gff3.gz"; fi
            find "$out/ncbi" -name '*_genomic.fna' | sort | xargs cat | gzip -n > "$out/genome.fna.gz.part"
            mv "$out/genome.fna.gz.part" "$out/genome.fna.gz"
            rm -rf "$out/ncbi" "$out/ncbi.zip"
        fi
    else
        fetch_url "$acc" "$out/genome.fna.gz" "$md5"
        url=${annot%% *}; if [[ -n $url ]]; then fetch_url "$url" "$out/annotation.gff3" ""; fi
    fi
done < <(rows config/genomes.tsv)
# transcripts: downloaded, or placed by hand when there is no public source (source "manual"); in both cases the
# file is checked by the md5 of its uncompressed content
while IFS=$'\t' read -r id genome desc source url md5 notes; do
    blank id genome desc source url md5 notes
    dest=data/transcripts/$id.fa.gz
    if [[ ! -s $dest && ( $source == manual || $url == PLACEHOLDER* ) ]]; then
        echo "MISSING $dest: no download source yet; place a copy there (plain or gzipped)"; fail=1; continue
    fi
    get_file "$url" "$dest" "$md5"
done < <(rows config/transcripts.tsv)

# other annotation files (config/extra_files.tsv)
while IFS=$'\t' read -r id genome desc url md5 dest notes; do
    blank id genome desc url md5 dest notes
    get_file "$url" "$dest" "$md5"
done < <(rows config/extra_files.tsv)
# references (Rfam model and seed; Selenoprofiles profiles)
fetch_url https://rfam.org/family/RF00031/cm data/references/RF00031.cm ""
fetch_url https://rfam.org/family/RF00031/alignment data/references/RF00031.sto ""
fetch_url https://rfam.org/family/RF01852/cm data/references/RF01852.cm ""
# Selenoprofiles data, as "selenoprofiles -download" assembles it: profiles, Gene Ontology tables and the UniRef50
# tag database formatted for BLAST (formatdb, from blast-legacy). The three archives are kept in
# data/downloads/selenoprofiles/ and checked against config/references.tsv when a checksum is recorded there.
SP=data/references/selenoprofiles_data
SP_URL=${SELENOPROFILES_DATA_URL:-https://mariottigenomicslab.bio.ub.edu/selenoprofiles_data}
if [[ ! -s $SP/uniref50_seleno_sim.fasta.psq ]]; then
    mkdir -p data/downloads/selenoprofiles
    ok=1
    for f in selenoprotein_profiles.tar.gz gene_ontology.tar.gz uniref50_seleno_sim.fasta.gz; do
        fetch_url "$SP_URL/$f" "data/downloads/selenoprofiles/$f" "$(awk -F'\t' -v f="$f" '$4 ~ "/"f"$" {print $6}' config/references.tsv)"
        [[ -s data/downloads/selenoprofiles/$f ]] || ok=0
    done
    if (( ok )) && ! command -v formatdb > /dev/null; then
        echo "MISSING $SP: formatdb (blast-legacy) is needed to format the UniRef50 tag database; rerun inside the image:"
        echo "        scripts/run_in_container.sh scripts/fetch.sh"
        ok=0; fail=1
    fi
    if (( ok )); then
        rm -rf "$SP.part"; mkdir -p "$SP.part"
        tar -xzf data/downloads/selenoprofiles/selenoprotein_profiles.tar.gz -C "$SP.part"
        tar -xzf data/downloads/selenoprofiles/gene_ontology.tar.gz -C "$SP.part"
        gzip -dc data/downloads/selenoprofiles/uniref50_seleno_sim.fasta.gz > "$SP.part/uniref50_seleno_sim.fasta"
        (cd "$SP.part" && formatdb -i uniref50_seleno_sim.fasta -p T -l formatdb.log)
        rm -rf "$SP"; mv "$SP.part" "$SP"
    else
        fail=1
    fi
fi
md5sum data/genomes/*/genome.fna.gz data/transcripts/*.fa.gz data/downloads/*/* data/downloads/*.zip 2>/dev/null > data/MD5SUMS || true
echo "Done. Checksums of downloaded files: data/MD5SUMS"
exit $fail
