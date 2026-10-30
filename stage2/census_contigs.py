#!/usr/bin/env python3
"""Contig checks for the census (step 9): composition of the contigs that carry Selenoprofiles predictions, compared
with the rest of the genome, and the Tiara classification of these contigs when it was run (scripts/run_tiara.sh).

For each genome of --genomes (data/genomes/<genome>/genome.fna.gz):
  - reference windows: up to --windows windows of each length of WINDOWS (2, 5 and 10 kb), at random positions on the
    contigs at least that long (fixed seed); for each, GC content and canonical tetranucleotide frequencies;
    the composition profile of a window is the centred log-ratio of its 136 canonical tetranucleotide counts (plus
    one), and its distance is the Euclidean distance to the component-wise median of the reference windows of the
    same length;
  - each contig carrying a prediction (--predictions, or a selenocysteine codon of --sites): its length, and for the
    window of the largest reference length that fits in it (centred on the middle of its predictions), GC content,
    composition distance, their percentiles among the reference windows of that length, and their robust z-scores
    ((x - median) / (1.4826 MAD) of the reference windows; contaminant windows among the reference windows raise
    the percentiles of contaminant contigs less than the robust z-scores); contigs shorter than the shortest window
    have no composition test;
  - the same GC and composition measures on the 2-kb window centred on the predictions, against the 2-kb reference
    windows: in fragmented assemblies, the longer reference windows come only from the few long contigs, which may
    themselves be contaminants, whereas the 2-kb windows sample most of the assembly;
  - the open reading frame density: the fraction of bases in stop-to-stop open reading frames of at least
    ORF_CODONS codons (six frames, N as a stop), on the window (with its robust z-score against the reference windows
    of the same length) and on the whole contig. Bacterial DNA is mostly coding (about 85 to 90%), animal genomic DNA
    much less; random sequence of 30 to 45% GC has few long open frames (GC-rich sequence has more);
  - Tiara classes (first and second stage) of the contig, read from <tiara-dir>/<genome>.tiara.tsv when present
    (Tiara classifies sequences of at least 3,000 nt; not run, decision of 2026-10-30).
Per genome (--genome-table): assembly length, contigs, N50; for each reference length, the fraction of the assembly
in contigs at least that long, the median GC and open reading frame density of the reference windows, and the 99th
percentile of their composition distances; the bases of the contigs of each Tiara first-stage class.
Outputs: --table (one row per genome and contig) and --genome-table.
"""
import argparse, collections, csv, gzip, os, random, sys

import numpy as np

GENOMES = ["OKI", "OSKA", "OdB3", "KUM", "AOM", "Bar2", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
WINDOWS = [2000, 5000, 10000]
ORF_CODONS = 100
CODE = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(b"ACGT"): CODE[b] = i; CODE[b + 32] = i     # upper and lower case


def canonical_index():
    """Tetranucleotide code (0-255) -> index of its canonical form (0-135)."""
    def rc(c):
        out = 0
        for _ in range(4):
            out = out * 4 + (3 - (c & 3)); c >>= 2
        return out
    canon = sorted({min(c, rc(c)) for c in range(256)})
    pos = {c: i for i, c in enumerate(canon)}
    return np.array([pos[min(c, rc(c))] for c in range(256)], dtype=np.int64)


CANON = canonical_index()


def composition(seq):
    """GC fraction (of A, C, G, T) and canonical tetranucleotide counts of a sequence (bytes)."""
    a = CODE[np.frombuffer(seq, dtype=np.uint8)]
    valid = a < 4
    n = int(valid.sum())
    gc = float(((a == 1) | (a == 2)).sum()) / n if n else float("nan")
    if len(a) < 4: return gc, np.zeros(136)
    v = valid[:-3] & valid[1:-2] & valid[2:-1] & valid[3:]
    a4 = a.astype(np.int64)
    code = a4[:-3] * 64 + a4[1:-2] * 16 + a4[2:-1] * 4 + a4[3:]
    return gc, np.bincount(CANON[code[v]], minlength=136).astype(float)


def orf_fraction(seq, min_codons=ORF_CODONS):
    """Fraction of the bases of a sequence (bytes) in stop-to-stop open reading frames of at least min_codons codons,
    in any of the six frames (TAA, TAG, TGA, and codons with a base other than A, C, G, T, end a frame)."""
    a = CODE[np.frombuffer(seq, dtype=np.uint8)]
    L = len(a)
    if L < 3 * min_codons: return 0.0
    cov = np.zeros(L, dtype=bool)
    rc = np.where(a < 4, 3 - a, 4)[::-1]
    for strand, s_ in ((0, a), (1, rc)):
        for f in range(3):
            n = (L - f) // 3
            c = s_[f:f + 3 * n].reshape(n, 3)
            stop = (c[:, 0] == 3) & (((c[:, 1] == 0) & ((c[:, 2] == 0) | (c[:, 2] == 2))) |
                                     ((c[:, 1] == 2) & (c[:, 2] == 0)))
            stop |= (c >= 4).any(axis=1)
            bounds = np.concatenate(([-1], np.flatnonzero(stop), [n]))
            lengths = np.diff(bounds) - 1
            for k in np.flatnonzero(lengths >= min_codons):
                st, en = f + 3 * (bounds[k] + 1), f + 3 * bounds[k + 1]
                if strand == 0: cov[st:en] = True
                else: cov[L - en:L - st] = True
    return float(cov.mean())


def clr(counts):
    x = np.log(counts + 1.0)
    return x - x.mean()


def read_fasta(path, keep=None):
    """Yield (name, sequence bytes) of a (gzipped) FASTA file; with keep, only those names (others: length only)."""
    opener = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None: yield name, b"".join(buf)
                name, buf = line[1:].split()[0].decode(), []
            else:
                buf.append(line.strip())
        if name is not None: yield name, b"".join(buf)


def percentile(values, x):
    """Percentage of values at or below x."""
    if not len(values) or x != x: return ""
    return f"{100.0 * np.searchsorted(np.sort(values), x, side='right') / len(values):.1f}"


def robust_z(values, x):
    """(x - median) / (1.4826 MAD) of the reference values: unlike percentiles, not diluted by contaminant windows
    among the reference windows."""
    if not len(values) or x != x: return ""
    med = float(np.median(values))
    mad = 1.4826 * float(np.median(np.abs(np.asarray(values) - med)))
    return f"{(x - med) / mad:.1f}" if mad > 0 else ""


def read_tiara(path):
    out = {}
    if not os.path.exists(path): return out
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0 and f[0].lower().startswith("sequence"): continue
            if len(f) >= 2: out[f[0].split()[0]] = (f[1], f[2] if len(f) > 2 else "")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", nargs="+", default=GENOMES)
    ap.add_argument("--genomes-dir", default="data/genomes")
    ap.add_argument("--predictions", default="results/tables/stage2_final_predictions.tsv")
    ap.add_argument("--sites", default="results/tables/stage2_final_sec_sites.tsv")
    ap.add_argument("--tiara-dir", default="work/stage2/census/tiara")
    ap.add_argument("--windows", type=int, default=2000, help="reference windows per length and genome")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--table", default="results/tables/stage2_census_contigs.tsv")
    ap.add_argument("--genome-table", default="results/tables/stage2_census_genomes.tsv")
    a = ap.parse_args()
    where = collections.defaultdict(lambda: collections.defaultdict(list))   # genome -> contig -> positions
    for path, cols in ((a.predictions, ("contig", "start", "end")), (a.sites, ("contig", "sec_codon_start", "sec_codon_end"))):
        if not os.path.exists(path): continue
        for r in csv.DictReader(open(path), delimiter="\t"):
            where[r["genome"]][r[cols[0]]].append((int(r[cols[1]]) + int(r[cols[2]])) // 2)
    rows, genome_rows = [], []
    for g in a.genomes:
        fasta = f"{a.genomes_dir}/{g}/genome.fna.gz"
        if not os.path.exists(fasta):
            print(f"{g}: {fasta} not found; skipped", file=sys.stderr); continue
        lengths = {n: len(s) for n, s in read_fasta(fasta)}
        rng = random.Random(f"{a.seed}-{g}")
        wanted = collections.defaultdict(list)                  # contig -> [(window length, start)]
        for w in WINDOWS:
            long_ = [(n, L) for n, L in lengths.items() if L >= w]
            if not long_: continue
            ends = np.cumsum([L - w + 1 for _, L in long_])
            for _ in range(a.windows):
                x = rng.randrange(int(ends[-1]))
                i = int(np.searchsorted(ends, x, side="right"))
                wanted[long_[i][0]].append((w, x - (int(ends[i - 1]) if i else 0)))
        ref = {w: {"gc": [], "clr": [], "orf": []} for w in WINDOWS}
        contigs = {}
        for n, s in read_fasta(fasta):
            for w, x in wanted.get(n, ()):
                gc, k = composition(s[x:x + w])
                if gc == gc:
                    ref[w]["gc"].append(gc); ref[w]["clr"].append(clr(k)); ref[w]["orf"].append(orf_fraction(s[x:x + w]))
            if n in where[g]:
                L = len(s)
                fit = [w for w in WINDOWS if w <= L]
                gc_all, _ = composition(s)
                rec = {"length": L, "gc": gc_all, "window": fit[-1] if fit else None, "orf": orf_fraction(s), "w": {}}
                mid = int(np.median(where[g][n]))
                for w in fit:
                    x = min(max(0, mid - w // 2), L - w)
                    wgc, k = composition(s[x:x + w])
                    rec["w"][w] = (wgc, clr(k), orf_fraction(s[x:x + w]))
                contigs[n] = rec
        stats = {}
        for w in WINDOWS:
            if not ref[w]["clr"]: continue
            m = np.median(np.array(ref[w]["clr"]), axis=0)
            d = np.array([np.linalg.norm(v - m) for v in ref[w]["clr"]])
            stats[w] = (np.array(ref[w]["gc"]), m, d, np.array(ref[w]["orf"]))
        tiara = read_tiara(f"{a.tiara_dir}/{g}.tiara.tsv")
        for n in sorted(contigs):
            rec = contigs[n]
            w = rec["window"]
            row = {"genome": g, "contig": n, "length": rec["length"], "gc": f"{rec['gc']:.4f}",
                   "predictions_and_sites": len(where[g][n]), "window": w or "", "window_gc": "",
                   "gc_percentile": "", "gc_robust_z": "", "composition_distance": "", "composition_percentile": "",
                   "composition_robust_z": "", "orf_fraction_window": "", "orf_robust_z": "",
                   "contig_orf_fraction": f"{rec['orf']:.3f}", "gc_robust_z_2kb": "", "composition_robust_z_2kb": "",
                   "orf_robust_z_2kb": "",
                   "tiara_first_stage": tiara.get(n, ("", ""))[0] if tiara else "not run",
                   "tiara_second_stage": tiara.get(n, ("", ""))[1] if tiara else ""}
            if w and w in stats:
                gcs, m, d, orfs = stats[w]
                wgc, wclr, worf = rec["w"][w]
                dist = float(np.linalg.norm(wclr - m))
                row.update({"window_gc": f"{wgc:.4f}", "gc_percentile": percentile(gcs, wgc),
                            "gc_robust_z": robust_z(gcs, wgc), "composition_distance": f"{dist:.3f}",
                            "composition_percentile": percentile(d, dist), "composition_robust_z": robust_z(d, dist),
                            "orf_fraction_window": f"{worf:.3f}", "orf_robust_z": robust_z(orfs, worf)})
            if WINDOWS[0] in rec["w"] and WINDOWS[0] in stats:
                gcs, m, d, orfs = stats[WINDOWS[0]]
                wgc, wclr, worf = rec["w"][WINDOWS[0]]
                row.update({"gc_robust_z_2kb": robust_z(gcs, wgc),
                            "composition_robust_z_2kb": robust_z(d, float(np.linalg.norm(wclr - m))),
                            "orf_robust_z_2kb": robust_z(orfs, worf)})
            if tiara and n not in tiara and rec["length"] < 3000: row["tiara_first_stage"] = "too short"
            rows.append(row)
        ls = sorted(lengths.values(), reverse=True)
        half, acc, n50 = sum(ls) / 2, 0, 0
        for L in ls:
            acc += L
            if acc >= half: n50 = L; break
        gr = {"genome": g, "assembly_length": sum(ls), "contigs": len(ls), "N50": n50}
        for w in WINDOWS:
            kb = f"{w // 1000}kb"
            gr[f"fraction_in_contigs_{kb}"] = f"{sum(L for L in ls if L >= w) / sum(ls):.3f}"
            gr[f"median_gc_{kb}"] = f"{np.median(stats[w][0]):.4f}" if w in stats else ""
            gr[f"median_orf_fraction_{kb}"] = f"{np.median(stats[w][3]):.3f}" if w in stats else ""
            gr[f"composition_distance_p99_{kb}"] = f"{np.percentile(stats[w][2], 99):.3f}" if w in stats else ""
        gr["tiara"] = "run" if tiara else "not run"
        if tiara:
            bases = collections.Counter()
            for n, (c1, _) in tiara.items(): bases[c1] += lengths.get(n, 0)
            gr["tiara_bases"] = "; ".join(f"{c} {b}" for c, b in bases.most_common())
        else:
            gr["tiara_bases"] = ""
        genome_rows.append(gr)
        print(f"{g}: {len(lengths)} contigs, {len(contigs)} with predictions", flush=True)
    for path, out in ((a.table, rows), (a.genome_table, genome_rows)):
        if not out: continue
        with open(path, "w") as o:
            wr = csv.DictWriter(o, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
            wr.writeheader(); wr.writerows(out)


if __name__ == "__main__":
    main()
