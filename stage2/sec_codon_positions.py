#!/usr/bin/env python3
"""Position of the selenocysteine codon in the gold-standard gene models, against the other codons of the same gene.

For each gold-standard model (curation/gold_standard/<genome>.gff3; coding segments and the selenocysteine codon of
transl_except), on the coding sequence in transcript order:
- the coding exon that holds the selenocysteine codon (first, internal, last, single; split when an intron interrupts
  the codon), the distance (nt) from the codon to the nearest coding exon junction (0: the codon touches a junction;
  split codons apart), and to the last junction when the codon lies upstream of it (the rule of exon-junction
  dependent nonsense-mediated decay: a stop codon more than about 50 nt upstream of the last junction elicits decay);
- codons_to_stop: the position of the stop codon counted from the selenocysteine codon (1: the stop codon follows
  directly);
- for each measure, where the selenocysteine codon falls among all the sense codons of the same coding sequence
  (excluding the start and stop codons): the mid-rank fraction of codons with a smaller value (P(X < x) + P(X = x)/2),
  which is uniform between 0 and 1 when the codon position is unrelated to the measure.
Summary (--summary) per set of models: the mean of these fractions over the genes with at least one coding intron
(expected 0.5; small values: the selenocysteine codon is closer to a junction, or to the stop codon, than the other
codons of the gene), with a one-sided p-value from the sum of independent uniforms (Irwin-Hall, normal approximation
above 30 genes), for one assembly (--assembly, default OSKA2016, the only one with the 20 genes complete) and with one
gene per family (paralogs of a family share their gene structure; the gene with the median value), and the counts
of split codons, of codons within --near nt of a junction, of codons in the last coding exon and of codons among the
last --tail codons, against their expectations (sums over genes of the fraction of codons of each gene in that class;
for the last count, the probability of at least the observed number, Poisson-binomial).
Outputs: results/tables/stage2_sec_codon_positions.tsv (one row per model) and _summary.tsv.
"""
import argparse, collections, math, re


def models(path):
    """Gene -> (strand, coding segments (1-based, genome order), selenocysteine codon (lo, hi))."""
    out = collections.defaultdict(lambda: ["+", [], None])
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        gene = re.search(r"Parent=rna-([^;]+)", f[8]).group(1).split("_", 1)[1]
        m = out[gene]
        m[0] = f[6]; m[1].append((int(f[3]), int(f[4])))
        s = re.search(r"pos:(?:complement\()?(\d+)\.\.(\d+)", f[8])
        if s: m[2] = tuple(map(int, s.groups()))
    return {g: (s, sorted(c), sec) for g, (s, c, sec) in out.items() if sec}


def coding_positions(cds, strand):
    out = []
    for x, y in (cds if strand == "+" else cds[::-1]):
        out += list(range(x, y + 1)) if strand == "+" else list(range(y, x - 1, -1))
    return out


def measures(offset, junctions, n_codons, last_exon_start):
    """For a codon starting at CDS offset (0-based): split, distance to the nearest junction, distance upstream of the
    last junction (None in the last exon), exon class, codons to the stop codon."""
    s, e = offset, offset + 3                      # junction j: first offset of the next coding exon
    split = any(s < j < e for j in junctions)
    dist = min((j - e if j >= e else s - j) for j in junctions) if junctions and not split else (0 if split else None)
    last = junctions[-1] if junctions else None
    up_last = (last - e) if last is not None and e <= last else None
    k = offset // 3
    return {"split": split, "dist": dist, "up_last": up_last, "in_last": last is None or s >= last,
            "to_stop": n_codons - 1 - k}


def mid_rank(x, values):
    lt = sum(v < x for v in values); eq = sum(v == x for v in values)
    return (lt + eq / 2) / len(values) if values else None


def irwin_hall_p(fracs):
    """P(sum of n independent uniforms <= observed sum) (exact for n <= 30, normal approximation above)."""
    n, s = len(fracs), sum(fracs)
    if n == 0: return None
    if n > 30: return 0.5 * math.erfc(-(s - n / 2) / math.sqrt(n / 12) / math.sqrt(2))
    return sum((-1) ** k * math.comb(n, k) * max(s - k, 0) ** n for k in range(n + 1)) / math.factorial(n)


FAMILY = lambda g: re.sub(r"_(PAR|chr\w+|pair\d)$", "", g)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--assembly", default="OSKA")
    ap.add_argument("--near", type=int, default=30, help="nt between the codon and a junction counted as near")
    ap.add_argument("--tail", type=int, default=3, help="codons before the stop codon counted as C-terminal")
    ap.add_argument("--out", default="results/tables/stage2_sec_codon_positions.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_sec_codon_positions_summary.tsv")
    a = ap.parse_args()
    rows = []
    for g in a.genomes:
        for gene, (strand, cds, sec) in sorted(models(f"{a.gold}/{g}.gff3").items()):
            pos = coding_positions(cds, strand)
            n_codons = len(pos) // 3
            lens = [y - x + 1 for x, y in (cds if strand == "+" else cds[::-1])]
            junctions = [sum(lens[:i + 1]) for i in range(len(lens) - 1)]
            first = pos.index(sec[0] if strand == "+" else sec[1])
            if first % 3: continue                       # not in frame: should not happen
            obs = measures(first, junctions, n_codons, junctions[-1] if junctions else 0)
            sense = [measures(3 * k, junctions, n_codons, 0) for k in range(1, n_codons - 1)]   # without start and stop
            row = {"genome": g, "gene": gene, "codons": n_codons, "coding_exons": len(cds),
                   "sec_codon": first // 3 + 1, "exon": ("single" if not junctions else "split" if obs["split"] else
                                                         "last" if obs["in_last"] else
                                                         "first" if first < junctions[0] else "internal"),
                   "dist_junction": "" if obs["dist"] is None else obs["dist"],
                   "upstream_of_last_junction": "" if obs["up_last"] is None else obs["up_last"],
                   "codons_to_stop": obs["to_stop"]}
            if junctions:
                row["frac_dist"] = mid_rank(obs["dist"], [m["dist"] for m in sense])
                row["frac_split"] = sum(m["split"] for m in sense) / len(sense)
                row["frac_near"] = sum(m["dist"] <= a.near for m in sense) / len(sense)
                row["frac_last_exon"] = sum(m["in_last"] for m in sense) / len(sense)
            row["frac_to_stop"] = mid_rank(obs["to_stop"], [m["to_stop"] for m in sense])
            rows.append(row)
    cols = ["genome", "gene", "codons", "coding_exons", "sec_codon", "exon", "dist_junction",
            "upstream_of_last_junction", "codons_to_stop", "frac_dist", "frac_to_stop", "frac_split", "frac_near",
            "frac_last_exon"]
    fmt = lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")

    def summarise(sel, label):
        out = []
        intr = [r for r in sel if r["coding_exons"] > 1]
        for meas, name in (("frac_dist", "distance to the nearest junction"), ("frac_to_stop", "codons to the stop codon")):
            vals = [r[meas] for r in (intr if meas == "frac_dist" else sel)]
            out.append([label, f"rank of the Sec codon by {name} (mean mid-rank; expected 0.5)", len(vals),
                        f"{sum(vals) / len(vals):.3f}" if vals else "", f"{irwin_hall_p(vals):.2g}" if vals else ""])
        for cls, frac, name in (("split", "frac_split", "codon split by an intron"),
                                ("near", "frac_near", f"codon within {a.near} nt of a junction (split included)"),
                                ("last", "frac_last_exon", "codon in the last coding exon")):
            obs = sum(1 for r in intr if (r["exon"] == "split" if cls == "split" else
                                          (r["exon"] == "split" or (r["dist_junction"] != "" and r["dist_junction"] <= a.near)) if cls == "near" else
                                          r["exon"] in ("last", "split") and r["upstream_of_last_junction"] == "" if cls == "last" else False))
            exp = sum(r[frac] for r in intr)
            out.append([label, name, len(intr), f"observed {obs}", f"expected {exp:.1f}"])
        # C-terminal codons: Sec codon among the last --tail codons before the stop codon, against a uniform position
        # among the sense codons of each gene (Poisson-binomial tail)
        probs = [min(a.tail, r["codons"] - 2) / (r["codons"] - 2) for r in sel]
        obs = sum(r["codons_to_stop"] <= a.tail for r in sel)
        dist = [1.0]
        for q in probs:
            new = [0.0] * (len(dist) + 1)
            for k, v in enumerate(dist): new[k] += v * (1 - q); new[k + 1] += v * q
            dist = new
        out.append([label, f"codon among the last {a.tail} before the stop codon", len(sel), f"observed {obs}",
                    f"expected {sum(probs):.2f}; P(X >= observed) = {sum(dist[obs:]):.1e}"])
        return out

    summ = []
    one = [r for r in rows if r["genome"] == a.assembly]
    summ += summarise(one, f"{a.assembly}, all genes")
    fam = collections.defaultdict(list)
    for r in one: fam[FAMILY(r["gene"])].append(r)
    per_fam = []
    for f, rs in fam.items():
        rs = sorted(rs, key=lambda r: r.get("frac_dist", 1) if r["coding_exons"] > 1 else 1)
        per_fam.append(rs[len(rs) // 2])
    summ += summarise(per_fam, f"{a.assembly}, one gene per family ({len(per_fam)} families)")
    with open(a.summary, "w") as o:
        o.write("set\tmeasure\tgenes\tvalue\tp_or_expected\n")
        for r in summ: o.write("\t".join(str(x) for x in r) + "\n")
    for r in summ: print("\t".join(str(x) for x in r))


if __name__ == "__main__":
    main()
