#!/usr/bin/env python3
"""Stage 1 and stage 2 settings of the pipeline measured against the reference sets of stage2/truth_sets.py.

For each run (--runs NAME=DIR: a pipeline output directory), each genome of the run (compare/<genome>.sec_sites.tsv)
and each gene of the reference set of that genome (--truth; gold standard, lifted models whatever their status, and
validated stage 1 predictions), the Sec site of the run at the reference codon (same contig, strand and lowest
coordinate) is looked up:
- sec_selenoprofiles: the Selenoprofiles prediction holding the codon, and its label (a Sec site counts when the
  label is selenocysteine); sec_miniprot: miniprot aligned a selenocysteine of a profile sequence to the codon;
  found: either (the Sec sites that the pipeline reports and that its expectation check counts);
- family: the family of the site (Selenoprofiles, else miniprot), and whether it is the reference family;
- locus: when the codon is missed, the labels of the Selenoprofiles predictions that overlap the reference coding
  sequence on the same strand (locus found, Sec codon missed);
- structure, for references with a gene model and a Selenoprofiles prediction at the codon: introns of the
  prediction identical to those of the reference (exact), and the fractions of reference coding bases covered by the
  prediction (sensitivity) and of predicted coding bases inside the reference (precision); the terminal stop codon
  is ignored on both sides;
- secis: the best SECIS score in the 3' sequence of the prediction (reported by the pipeline only above the threshold
  of the model in that run).
The summary counts the gold-standard and stage 1 references and the lifted models read as complete (the gene table
lists every reference with a codon). Sec sites of a run (selenocysteine label or miniprot) at no reference codon
are listed as other sites, for every genome of the run, including those without reference (held-out genomes).
Outputs: --genes (one row per run, genome and reference gene), --summary (one row per run and genome), --other.
"""
import argparse, collections, csv, os, re, sys

RUNS = ["stage1=work/stage1/pipeline", "default=work/stage2/pipeline/default", "eval=work/stage2/pipeline/eval",
        "eval_miniprot=work/stage2/pipeline/eval_miniprot", "loo_OKI=work/stage2/pipeline/loo_OKI",
        "loo_OSKA=work/stage2/pipeline/loo_OSKA", "loo_OdB3=work/stage2/pipeline/loo_OdB3"]
GENE_ID = re.compile(r'gene_id "([^"]+)"')
GENOMES = ["OKI", "OSKA", "OdB3", "KUM", "AOM", "Bar2", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]


def usable(r):
    """References counted in the summary: gold standard, stage 1 references, and lifts read as complete."""
    return r["status"] in ("gold", "reference") or r["status"].startswith("complete")


def segments(text):
    return sorted(tuple(map(int, x.split("-"))) for x in text.split(";") if x) if text else []


def introns(cds):
    return {(a[1], b[0]) for a, b in zip(cds, cds[1:]) if b[0] - a[1] > 1}


def bases(cds):
    return {i for a, b in cds for i in range(a, b + 1)}


def strip_stop(cds, strand):
    """Coding segments without the last 3 bases (stop codon), on the strand of the gene."""
    if not cds: return cds
    cds = [list(x) for x in cds]
    left = 3
    order = range(len(cds) - 1, -1, -1) if strand == "+" else range(len(cds))
    for i in order:
        a, b = cds[i]
        n = min(left, b - a + 1)
        if strand == "+": cds[i][1] -= n
        else: cds[i][0] += n
        left -= n
        if not left: break
    return [tuple(x) for x in cds if x[1] >= x[0]]


def read_sites(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[(r["contig"], int(r["start"]), r["strand"])] = r
    return out


def read_gtf(path):
    cds = collections.defaultdict(list)
    where = {}
    if not os.path.exists(path): return None, None
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        m = GENE_ID.search(f[8])
        if not m: continue
        cds[m.group(1)].append((int(f[3]), int(f[4])))
        where[m.group(1)] = (f[0], f[6])
    return {k: sorted(v) for k, v in cds.items()}, where


def is_sec(site):
    return site.get("selenoprofiles_label") == "selenocysteine" or site.get("support") in ("miniprot", "both")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=RUNS)
    ap.add_argument("--truth", default="results/tables/stage2_truth_sets.tsv")
    ap.add_argument("--genomes", nargs="+", default=GENOMES, help="genomes reported (others of the runs are skipped)")
    ap.add_argument("--genes", default="results/tables/stage2_benchmark_genes.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_benchmark_summary.tsv")
    ap.add_argument("--other", default="results/tables/stage2_benchmark_other_sites.tsv")
    a = ap.parse_args()
    truth = collections.defaultdict(list)
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        if r["sec_start"]: truth[r["genome"]].append(r)
    gene_rows, summary, other = [], [], []
    for spec in a.runs:
        run, d = spec.split("=", 1)
        if not os.path.isdir(f"{d}/compare"):
            print(f"{run}: no {d}/compare, skipped", file=sys.stderr); continue
        present = {f[:-len(".sec_sites.tsv")] for f in os.listdir(f"{d}/compare") if f.endswith(".sec_sites.tsv")}
        genomes = [g for g in a.genomes if g in present]
        for g in genomes:
            sites = read_sites(f"{d}/compare/{g}.sec_sites.tsv")
            preds, where = read_gtf(f"{d}/selenoprofiles/{g}.selenoprofiles.gtf")
            refs = truth.get(g, [])
            ref_keys = {(r["contig"], int(r["sec_start"]), r["strand"]) for r in refs}
            n = collections.Counter()
            for r in refs:
                key = (r["contig"], int(r["sec_start"]), r["strand"])
                s = sites.get(key, {})
                sp_id = s.get("selenoprofiles_id", ".") if s else "."
                label = s.get("selenoprofiles_label", ".") if s else "."
                sp = label == "selenocysteine"
                mp = bool(s) and s.get("support") in ("miniprot", "both")
                fam = s.get("family", "") if s else ""
                row = {"run": run, "genome": g, "gene": r["gene"], "profile": r["profile"], "reference": r["source"],
                       "reference_status": r["status"], "found": "yes" if (sp or mp) else "no",
                       "sec_selenoprofiles": label if sp_id != "." else "", "sec_miniprot": "yes" if mp else "no",
                       "prediction": sp_id if sp_id != "." else "", "family": fam,
                       "same_family": ("yes" if fam == r["profile"] else "no") if fam else "",
                       "locus": "", "introns_exact": "", "cds_sensitivity": "", "cds_precision": "",
                       "secis": s.get("secis_score", "") if s and s.get("secis_score") not in (".", None) else ""}
                ref_cds = segments(r["cds"])
                if preds is not None and ref_cds:
                    if sp_id == ".":
                        lo, hi = ref_cds[0][0], ref_cds[-1][1]
                        labels = sorted({p.rsplit(".", 1)[1] for p, c in preds.items()
                                         if where[p] == (r["contig"], r["strand"]) and bases(c) & bases(ref_cds)})
                        row["locus"] = ",".join(labels) if labels else "none"
                    elif sp_id in preds:
                        pc = strip_stop(preds[sp_id], r["strand"])
                        rc_ = strip_stop(ref_cds, r["strand"])
                        bp, br = bases(pc), bases(rc_)
                        row["introns_exact"] = "yes" if introns(pc) == introns(rc_) else "no"
                        row["cds_sensitivity"] = f"{len(bp & br) / len(br):.3f}"
                        row["cds_precision"] = f"{len(bp & br) / len(bp):.3f}" if bp else ""
                gene_rows.append(row)
                if not usable(r): continue
                n["genes"] += 1
                n["found"] += row["found"] == "yes"
                n["selenoprofiles"] += sp
                n["miniprot"] += mp
                n["both"] += sp and mp
                n["same_family"] += row["same_family"] == "yes"
                n["locus_only"] += row["found"] == "no" and row["locus"] not in ("", "none")
                n["introns_exact"] += row["introns_exact"] == "yes"
                n["secis"] += bool(row["secis"])
            extra = [(k, s) for k, s in sites.items() if is_sec(s) and k not in ref_keys]
            for (contig, start, strand), s in sorted(extra):
                other.append({"run": run, "genome": g, "contig": contig, "start": start, "strand": strand,
                              "family": s["family"], "selenoprofiles_id": s["selenoprofiles_id"],
                              "selenoprofiles_label": s["selenoprofiles_label"], "support": s["support"],
                              "secis": s["secis_score"] if s["secis_score"] != "." else ""})
            fams = collections.Counter(s["family"] for _, s in extra)
            summary.append({"run": run, "genome": g, "reference_genes": n["genes"], "found": n["found"],
                            "found_selenoprofiles": n["selenoprofiles"], "found_miniprot": n["miniprot"],
                            "found_both": n["both"], "same_family": n["same_family"],
                            "locus_found_codon_missed": n["locus_only"],
                            "introns_exact": n["introns_exact"] if preds is not None else "",
                            "secis": n["secis"], "sec_sites": sum(1 for s in sites.values() if is_sec(s)),
                            "other_sec_sites": len(extra),
                            "other_families": ", ".join(f"{f} {c}" for f, c in sorted(fams.items()))})
    for path, rows in ((a.genes, gene_rows), (a.summary, summary), (a.other, other)):
        if not rows: continue
        with open(path, "w") as o:
            w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)
    for s in summary:
        print("\t".join(str(s[k]) for k in ("run", "genome", "reference_genes", "found", "found_selenoprofiles",
                                            "found_miniprot", "introns_exact", "secis", "other_sec_sites")))


if __name__ == "__main__":
    main()
