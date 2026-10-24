#!/usr/bin/env python3
"""Step 1 of the plan (docs/plan.md, 2026-10-24): SECIS search in the validated selenoprotein hits.

Hits: the stage 1 Selenoprofiles predictions with a selenocysteine codon (label "selenocysteine") in the six
Oikopleura assemblies (OKI2018, OSKA2016, OdB3, O. longicauda, O. albicans, O. vanhoeffeni); B. stygius,
M. erythrocephalus and F. borealis are held out (decision of 2026-10-24) and reported apart.
A hit is left out only on evidence against it (hit validation, round 1):
  - its best human hit (reciprocal best-hit tables) is in an unrelated family, that is, neither a human gene of its
    own family nor a related gene (RELATED: MIEN1 for SELENOW, SEPHS1 for SEPHS2); "no hit" is not evidence against;
  - in the tree of its family, its nearest relatives (the smallest side of a bipartition that holds it and other
    than appendicularian leaves) are all bacteria, archaea, fungi or other non-animal groups.
The curated O. dioica genes are always kept.
Scores (no new search): the best score in the 6,000 nt after each prediction with Rfam RF00031
(results/tables/stage1_secis_scores.tsv) and with the current appendicularian model
(results/tables/stage2_secis_scores.tsv), against the 200 random windows of each genome. The current model was
trained on the selenocysteine predictions of OKI2018, OdB3, O. longicauda, O. albicans and O. vanhoeffeni (and of the
two held-out oikopleurids), so that its scores there are overstated; OSKA2016 was held out of its training.
Cysteine homologs of the same genomes are given as a second negative control.
Outputs: --table (one row per hit), --summary (one row per genome and model) and --families (one row per
family, kept hits of the six Oikopleura assemblies).
"""
import argparse, collections, csv, glob, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import machinery_lba as lba
from machinery_rbh import MACHINERY_GENES, SELENOPROTEIN_GENES
from selenoprotein_orthology import baseline_names, sides

TRAINING = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"]
HELD_OUT = ["Bsty", "Mery", "Fbor"]
CURRENT_MODEL_TRAINING = {"OKI", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery"}
RELATED = {"SELENOW": ["MIEN1"], "SEPHS2": ["SEPHS1"]}
NONANIMAL = {"Bacteria", "Archaea", "Fungi", "Green algae", "Red algae", "Stramenopiles", "Alveolates", "Amoebozoans",
             "Kinetoplastids", "Cryptophytes", "Heteroloboseans", "Apusozoans"}


def scores(path):
    """Best score and length of the 3' sequence of each prediction (Selenoprofiles writes up to 6,000 nt, less at
    the end of a contig), and the best scores of the random windows of each genome (0 when no hit)."""
    out, length, background = {}, {}, collections.defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        b = float(r["best_bits"]) if r["best_bits"] not in ("", ".") else None
        m = re.match(r"region_of_(\d+)_bp_downstream_of:(.+)", r["sequence"])
        if r["group"] == "random window": background[r["genome"]].append(b or 0.0)
        elif m:
            out[(r["genome"], m.group(2))] = b; length[(r["genome"], m.group(2))] = int(m.group(1))
    return out, length, background


def nonanimal_neighbours(trees, taxa):
    """(genome, prediction) of appendicularian leaves whose nearest relatives are all non-animal."""
    flagged = set()
    for tf in glob.glob(f"{trees}/*.treefile"):
        fam = os.path.basename(tf)[:-len(".treefile")]
        groups = dict(l.rstrip("\n").split("\t") for l in open(f"{trees}/{fam}.groups.tsv"))
        names, sd = sides(lba.read_tree(tf))
        tax = {n: lba.taxon(n, groups, taxa) for n in names}
        for n in names:
            if tax[n] != "Appendicularians": continue
            near = min((s for s in sd if n in s and any(tax[x] != "Appendicularians" for x in s)), key=len)
            others = {tax[x] for x in near if tax[x] != "Appendicularians"}
            if others and others <= NONANIMAL:
                g, _, p = n.partition("_"); flagged.add((g, p))
    return flagged


def pct(values, q):
    v = sorted(values)
    return v[min(len(v) - 1, int(round(q * (len(v) - 1))))] if v else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rfam", default="results/tables/stage1_secis_scores.tsv")
    ap.add_argument("--current", default="results/tables/stage2_secis_scores.tsv")
    ap.add_argument("--rbh", nargs="+", default=["results/tables/stage2_selenoprotein_rbh.tsv",
                                                 "results/tables/stage1_machinery_rbh.tsv"])
    ap.add_argument("--trees", default="results/trees/selenoproteins")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--baseline", default="curation/baseline_sec_genes.tsv")
    ap.add_argument("--annotation", default="results/tables/stage1_annotation.tsv")
    ap.add_argument("--table", default="results/tables/stage2_secis_hits.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_secis_hits_summary.tsv")
    ap.add_argument("--families", default="results/tables/stage2_secis_hits_families.tsv")
    a = ap.parse_args()
    rfam, length, bg_rfam = scores(a.rfam)
    cur_model, _, bg_cur = scores(a.current)
    rbh = {}
    for f in a.rbh:
        for r in csv.DictReader(open(f), delimiter="\t"): rbh.setdefault((r["genome"], r["prediction"]), r)
    curated = baseline_names(a.baseline, a.annotation)
    flagged = nonanimal_neighbours(a.trees, lba.load_taxa(a.taxa))
    genes = {**SELENOPROTEIN_GENES, **MACHINERY_GENES}
    rows = []
    for (g, pred), b in sorted(rfam.items()):
        if g not in TRAINING + HELD_OUT: continue
        label = pred.rsplit(".", 1)[-1]
        if label not in ("selenocysteine", "cysteine"): continue
        fam = pred.split(".")[0]
        r = rbh.get((g, pred))
        best = r["best_human_gene"] if r else "not in the tables"
        status = r["status"] if r else "-"
        reasons = []
        if g in HELD_OUT: reasons.append("held-out genome")
        if r and best not in ("-", "") and best not in genes.get(fam, []) + RELATED.get(fam, []):
            reasons.append(f"best human hit in another family ({best})")
        if (g, pred) in flagged: reasons.append("nearest relatives non-animal")
        name = curated.get((g, pred), "")
        if name: reasons = [x for x in reasons if x == "held-out genome"]   # curated genes kept
        rows.append({"genome": g, "prediction": pred, "family": fam, "label": label, "curated gene": name or "-",
                     "best human hit": best, "status": status,
                     "decision": "kept" if not reasons else "left out", "reason": "; ".join(reasons) or "-",
                     "3' sequence (nt)": length.get((g, pred), ""),
                     "RF00031 bits": "" if b is None else f"{b:.1f}",
                     "current model bits": "" if cur_model.get((g, pred)) is None else f"{cur_model[(g, pred)]:.1f}",
                     "in the training of the current model": "yes" if g in CURRENT_MODEL_TRAINING else "no"})
    cols = list(rows[0])
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")

    summ = [("genome", "model", "selenocysteine hits kept", "left out", "kept with 3' sequence under 500 nt",
             "with a hit", "at least 20 bits",
             "at least 25 bits", "above the best random window", "cysteine homologs", "cysteine homologs at least 20 bits",
             "best random window", "99th percentile of random windows", "random windows at least 20 bits")]
    for g in TRAINING + HELD_OUT:
        for model, col, bg in (("RF00031", "RF00031 bits", bg_rfam), ("current", "current model bits", bg_cur)):
            sec = [r for r in rows if r["genome"] == g and r["label"] == "selenocysteine"]
            kept = [r for r in sec if r["decision"] == "kept"]
            cys = [r for r in rows if r["genome"] == g and r["label"] == "cysteine"]
            val = lambda r: float(r[col]) if r[col] else None
            bmax = max(bg[g]) if bg[g] else 0.0
            kv = [val(r) for r in kept]
            summ.append((g, model + (" (training genome)" if model == "current" and g in CURRENT_MODEL_TRAINING else ""),
                         len(kept), len(sec) - len(kept), sum(int(r["3' sequence (nt)"] or 0) < 500 for r in kept),
                         sum(v is not None for v in kv),
                         sum(v is not None and v >= 20 for v in kv), sum(v is not None and v >= 25 for v in kv),
                         sum(v is not None and v > bmax for v in kv), len(cys),
                         sum(val(r) is not None and val(r) >= 20 for r in cys), f"{bmax:.1f}", f"{pct(bg[g], 0.99):.1f}",
                         sum(x >= 20 for x in bg[g])))
    with open(a.summary, "w") as o:
        for r in summ: o.write("\t".join(map(str, r)) + "\n")
    # per family, kept selenocysteine hits of the six Oikopleura assemblies
    fam_rows = [("family", "kept hits", "3' sequence under 500 nt", "RF00031 at least 20 bits",
                 "current model at least 25 bits (OSKA2016, held out)", "current model at least 25 bits (training genomes)")]
    kept = [r for r in rows if r["decision"] == "kept" and r["label"] == "selenocysteine" and r["genome"] in TRAINING]
    f = lambda r, c: float(r[c]) if r[c] else -1
    for fam in sorted({r["family"] for r in kept}):
        k = [r for r in kept if r["family"] == fam]
        osk = [r for r in k if r["genome"] == "OSKA"]; trn = [r for r in k if r["genome"] != "OSKA"]
        fam_rows.append((fam, len(k), sum(int(r["3' sequence (nt)"] or 0) < 500 for r in k),
                         sum(f(r, "RF00031 bits") >= 20 for r in k),
                         f"{sum(f(r, 'current model bits') >= 25 for r in osk)} of {len(osk)}",
                         f"{sum(f(r, 'current model bits') >= 25 for r in trn)} of {len(trn)}"))
    with open(a.families, "w") as o:
        for r in fam_rows: o.write("\t".join(map(str, r)) + "\n")
    print(open(a.families).read())
    for r in rows:
        if r["decision"] == "left out" and r["genome"] in TRAINING and r["label"] == "selenocysteine":
            print("left out:", r["genome"], r["prediction"], r["reason"])
    print(open(a.summary).read())


if __name__ == "__main__":
    main()
