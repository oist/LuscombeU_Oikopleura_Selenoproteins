#!/usr/bin/env python3
"""Genes that end with a TGA followed closely by another stop codon, and SECIS-like hits after them.

A path to a C-terminal selenoprotein (question of 2026-10-27): a gene whose stop codon is a TGA followed, within one or
two codons, by another stop codon in the same frame would gain one to three residues, the first a selenocysteine, if a
SECIS element appeared downstream (as KKSU, whose TGA is followed directly by a stop codon). How often is this path
available, and are SECIS-like hits more frequent after such genes than after the others?

For each genome of --genomes with an annotation (stage1/samplesheet.csv), one coding sequence per gene (the longest;
genes: the gene attribute, else the ANISEED transcript name without its isoform suffix, else the transcript):
1. Stop codons. The primary stop codon (the last codon of the CDS rows, or the codon after them when the CDS rows end
   before it; genes with neither are left out) and the in-frame codons after it, read in the genome: the first stop
   codon among the next --look codons (1: directly after the primary stop codon). Background: the frequency of stop
   codons per in-frame codon at positions --bg-from to --bg-to after the primary stop codon. Per primary stop codon
   (TAA, TAG, TGA): genes with a stop codon at +1, +2, +3, and the enrichment at +1 over the background.
2. SECIS-like hits. In the --window nt after the primary stop codon on the gene strand (cut at the next annotated
   coding sequence on that strand), the best hit of the genome scan of step 4 (Oikopleura SECIS model, cmsearch, hits
   of 15 bits or more: work/stage2/uga_screen/scan/<genome>.genome.tbl) and, at the same coordinates, of the scan of
   the genome shuffled in windows of 1,000 nt (<genome>.shuffled.tbl; the shuffling keeps the windows in place). Genes
   in classes: TGA with a stop codon at +1 or +2 (the path), TGA with the next stop farther, TAA or TAG; known
   selenoprotein genes (overlapping a gold-standard model, or with an annotated selenocysteine) apart. Per class:
   windows with a hit of at least --bits (20, the threshold of 3' windows of step 3) and at least the genome threshold
   of step 4 (results/tables/stage2_secis_scan_summary.tsv), in the genome and in the shuffled copy, and a Fisher
   exact test of the path against TAA or TAG (genome scan, --bits).
3. Candidates: genes of the path class (and known selenoprotein genes, as positive controls) with a hit at the genome
   threshold or above.
Outputs: results/tables/stage2_tandem_stops_codons.tsv (1), _secis.tsv (2), _candidates.tsv (3); one row per gene in
work/stage2/tandem_stops/<genome>.genes.tsv.
"""
import argparse, bisect, collections, csv, gzip, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import read_fasta, rc, opener, read_gff_transcripts, read_tbl, attributes   # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
GENOMES = ["OKI", "OSKA", "OdB3", "Cint", "Crob", "Scla", "Tdem", "Bsch", "Pmam", "Bflo", "Blan", "Hsap", "Mmus", "Drer"]


def gene_key(tid, t):
    if t.get("gene"): return t["gene"]
    m = re.match(r"(.+)\.\d+\.t$", tid)              # ANISEED: Oidioi.CG.OSKA2016.S27.g15032.01.t
    return m.group(1) if m else tid


def products(gff, ids):
    """product attribute of the CDS rows of the transcripts ids."""
    out = {}
    for l in opener(gff):
        if "\tCDS\t" not in l: continue
        f = l.rstrip("\n").split("\t")
        at = attributes(f[8])
        p = at.get("Parent", "").split(",")[0]
        if p in ids and p not in out: out[p] = at.get("product", "").replace("%2C", ",")
    return out


def gold_spans(path):
    out = collections.defaultdict(list)
    if not os.path.exists(path): return out
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "mRNA": out[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    return out


class Hits:
    """Scan hits indexed by sequence and strand: best score in a window (0-based half-open; hits starting inside)."""

    def __init__(self, path, min_bits=15.0):
        self.ok = os.path.exists(path)
        idx = collections.defaultdict(list)
        for t, st, x, y, bits, ev in read_tbl(path, min_bits):
            idx[(t, st)].append((min(x, y) - 1, max(x, y), bits))
        self.idx = {k: sorted(v) for k, v in idx.items()}
        self.starts = {k: [h[0] for h in v] for k, v in self.idx.items()}

    def best(self, contig, strand, lo, hi):
        k = (contig, strand)
        if k not in self.idx: return 0.0, None
        v, s = self.idx[k], self.starts[k]
        i, j = bisect.bisect_left(s, lo), bisect.bisect_left(s, hi)
        if i >= j: return 0.0, None
        b = max(v[i:j], key=lambda h: h[2])
        return b[2], b


def fisher(a, b, c, d):
    try:
        from scipy.stats import fisher_exact
        return fisher_exact([[a, b], [c, d]], alternative="greater")[1]
    except Exception:
        return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", nargs="+", default=GENOMES)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--scan", default="work/stage2/uga_screen/scan")
    ap.add_argument("--scan-summary", default="results/tables/stage2_secis_scan_summary.tsv")
    ap.add_argument("--look", type=int, default=3)
    ap.add_argument("--bg-from", type=int, default=4)
    ap.add_argument("--bg-to", type=int, default=33)
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--bits", type=float, default=20.0)
    ap.add_argument("--work", default="work/stage2/tandem_stops")
    ap.add_argument("--prefix", default="results/tables/stage2_tandem_stops")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    thr = {}
    if os.path.exists(a.scan_summary):
        for r in csv.DictReader(open(a.scan_summary), delimiter="\t"):
            if r["set"] == "genome" and r.get("genome threshold (bits)", "-") not in ("", "-"):
                thr[r["genome"]] = float(r["genome threshold (bits)"])
    codon_rows, secis_rows, cand_rows = [], [], []
    for g in a.genomes:
        if g not in sheet or not sheet[g].get("annotation") or not os.path.exists(sheet[g]["annotation"]):
            print(f"{g}: no annotation, left out", file=sys.stderr); continue
        txs = read_gff_transcripts(sheet[g]["annotation"])
        best = {}
        for k, t in txs.items():
            if not t["cds"]: continue
            n = sum(y - x for x, y in t["cds"])
            gk = gene_key(k, t)
            if gk not in best or n > best[gk][1]: best[gk] = (k, n)
        chosen = {k for k, _ in best.values()}
        # coding starts and ends per sequence and strand, to cut the windows at the next coding sequence
        cds_lo, cds_hi = collections.defaultdict(list), collections.defaultdict(list)
        for k, t in txs.items():
            if t["cds"]:
                cds_lo[(t["contig"], t["strand"])].append(t["cds"][0][0]); cds_hi[(t["contig"], t["strand"])].append(t["cds"][-1][1])
        for v in list(cds_lo.values()) + list(cds_hi.values()): v.sort()
        gold = gold_spans(f"{a.gold}/{g}.gff3")
        real, shuf = Hits(f"{a.scan}/{g}.genome.tbl"), Hits(f"{a.scan}/{g}.shuffled.tbl")
        by_contig = collections.defaultdict(list)
        for k in chosen: by_contig[txs[k]["contig"]].append(k)
        genes = []
        for contig, seq in read_fasta(sheet[g]["fasta"]):
            if contig not in by_contig: continue
            seq = seq.upper()
            for k in by_contig[contig]:
                t = txs[k]
                st = t["strand"]
                lo, hi = t["cds"][0][0], t["cds"][-1][1]
                if st == "+":
                    if seq[hi - 3:hi] in STOPS: s0 = hi - 3
                    elif seq[hi:hi + 3] in STOPS: s0 = hi
                    else: continue
                    stop = seq[s0:s0 + 3]
                    after = [seq[s0 + 3 * i:s0 + 3 * i + 3] for i in range(1, a.bg_to + 1)]
                    w_lo = s0 + 3
                    nxt = cds_lo[(contig, st)]
                    j = bisect.bisect_right(nxt, w_lo)
                    w_hi = min(len(seq), w_lo + a.window, nxt[j] if j < len(nxt) else len(seq))
                else:
                    if rc(seq[lo:lo + 3]) in STOPS: s0 = lo
                    elif lo >= 3 and rc(seq[lo - 3:lo]) in STOPS: s0 = lo - 3
                    else: continue
                    stop = rc(seq[s0:s0 + 3])
                    after = [rc(seq[s0 - 3 * i:s0 - 3 * i + 3]) if s0 - 3 * i >= 0 else "" for i in range(1, a.bg_to + 1)]
                    w_hi = s0
                    prv = cds_hi[(contig, st)]
                    j = bisect.bisect_left(prv, w_hi) - 1
                    w_lo = max(0, w_hi - a.window, prv[j] if j >= 0 else 0)
                nxt_stop = next((i + 1 for i, c in enumerate(after[:a.look]) if c in STOPS), "")
                bg = after[a.bg_from - 1:a.bg_to]
                known = bool(t["sec"]) or any(x < hi and lo < y for x, y in gold.get((contig, st), []))
                rb, rh = real.best(contig, st, w_lo, w_hi)
                sb, _ = shuf.best(contig, st, w_lo, w_hi)
                dist = ""
                if rh: dist = (rh[0] - w_lo) if st == "+" else (w_hi - rh[1])
                cls = ("known selenoprotein" if known else
                       f"TGA, stop at +1 or +2" if stop == "TGA" and nxt_stop in (1, 2) else
                       "TGA, next stop farther" if stop == "TGA" else "TAA or TAG")
                genes.append({"genome": g, "transcript": k, "gene": gene_key(k, t), "protein": t.get("protein", ""),
                              "contig": contig, "strand": st, "stop": stop, "next_stop": nxt_stop,
                              "bg_stops": sum(c in STOPS for c in bg), "bg_codons": sum(1 for c in bg if len(c) == 3),
                              "window": w_hi - w_lo, "best_bits": rb, "hit_distance": dist, "shuffled_bits": sb,
                              "class": cls, "known": known})
        cols = ["genome", "transcript", "gene", "protein", "contig", "strand", "stop", "next_stop", "bg_stops",
                "bg_codons", "window", "best_bits", "hit_distance", "shuffled_bits", "class", "known"]
        with open(f"{a.work}/{g}.genes.tsv", "w") as o:
            o.write("\t".join(cols) + "\n")
            for r in genes: o.write("\t".join(str(r[c]) for c in cols) + "\n")
        # 1. stop codons
        for stop in ("TAA", "TAG", "TGA"):
            sel = [r for r in genes if r["stop"] == stop and not r["known"]]
            if not sel: continue
            n = len(sel)
            at = {i: sum(r["next_stop"] == i for r in sel) for i in (1, 2, 3)}
            bgf = sum(r["bg_stops"] for r in sel) / max(1, sum(r["bg_codons"] for r in sel))
            codon_rows.append([g, stop, n, f"{100 * n / sum(1 for r in genes if not r['known']):.1f}",
                               at[1], f"{100 * at[1] / n:.1f}", at[2], f"{100 * at[2] / n:.1f}", at[3], f"{100 * at[3] / n:.1f}",
                               f"{100 * bgf:.1f}", f"{(at[1] / n) / bgf:.2f}" if bgf else ""])
        # 2. SECIS-like hits
        if real.ok:
            t_g = thr.get(g)
            classes = ["TGA, stop at +1 or +2", "TGA, next stop farther", "TAA or TAG", "known selenoprotein"]
            counts = {}
            for c in classes:
                sel = [r for r in genes if r["class"] == c]
                n = len(sel)
                rb = sum(r["best_bits"] >= a.bits for r in sel); sb = sum(r["shuffled_bits"] >= a.bits for r in sel)
                rt = sum(r["best_bits"] >= t_g for r in sel) if t_g else ""; stt = sum(r["shuffled_bits"] >= t_g for r in sel) if t_g else ""
                counts[c] = (n, rb)
                secis_rows.append([g, c, n, rb, f"{100 * rb / n:.2f}" if n else "", sb, f"{100 * sb / n:.2f}" if n else "",
                                   t_g if t_g else "", rt, stt, ""])
            n1, k1 = counts["TGA, stop at +1 or +2"]; n2, k2 = counts["TAA or TAG"]
            if n1 and n2:
                secis_rows[-4][-1] = f"{fisher(k1, n1 - k1, k2, n2 - k2):.2g}"
            for r in genes:
                if t_g and r["best_bits"] >= t_g and r["class"] in ("TGA, stop at +1 or +2", "known selenoprotein") \
                        and r["stop"] == "TGA" and r["next_stop"] in (1, 2):
                    cand_rows.append(r)
        print(f"{g}: {len(genes)} genes; scan {'found' if real.ok else 'missing'}", flush=True)
    # products of the candidates
    for g in {r["genome"] for r in cand_rows}:
        prod = products(sheet[g]["annotation"], {r["transcript"] for r in cand_rows if r["genome"] == g})
        for r in cand_rows:
            if r["genome"] == g: r["product"] = prod.get(r["transcript"], "")
    with open(f"{a.prefix}_codons.tsv", "w") as o:
        o.write("genome\tstop_codon\tgenes\tpercent_of_genes\tstop_at_+1\tpercent\tstop_at_+2\tpercent\tstop_at_+3\t"
                f"percent\tbackground_stops_per_codon_percent (+{a.bg_from} to +{a.bg_to})\tenrichment_at_+1\n")
        for r in codon_rows: o.write("\t".join(str(x) for x in r) + "\n")
    with open(f"{a.prefix}_secis.tsv", "w") as o:
        o.write(f"genome\tclass\tgenes\twindows_with_hit_>={a.bits:g}_bits\tpercent\tshuffled_>={a.bits:g}\tpercent\t"
                "genome_threshold\twindows_at_threshold\tshuffled_at_threshold\tfisher_p_path_vs_TAA_TAG\n")
        for r in secis_rows: o.write("\t".join(str(x) for x in r) + "\n")
    cols = ["genome", "transcript", "gene", "protein", "product", "contig", "strand", "stop", "next_stop", "best_bits",
            "hit_distance", "shuffled_bits", "class"]
    with open(f"{a.prefix}_candidates.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sorted(cand_rows, key=lambda r: (r["genome"], -r["best_bits"])):
            o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"-> {a.prefix}_codons.tsv, _secis.tsv, _candidates.tsv ({len(cand_rows)} candidates)")


if __name__ == "__main__":
    main()
