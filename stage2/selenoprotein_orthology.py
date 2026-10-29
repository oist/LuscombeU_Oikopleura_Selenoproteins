#!/usr/bin/env python3
"""First-pass orthology of the appendicularian selenoproteins, from the trees and the reciprocal best hits.

For each family tree (stage2/machinery_phylogeny.py, make selenoprotein-phylogeny; and the LG+C60+F+G PMSF tree of
make selenoprotein-pmsf when present), the appendicularian sequences are split into appendicularian clades: the
largest groups of appendicularian leaves that form a clade (a side of a bipartition) on their own; a single leaf is
a clade of one. For each clade, read on the unrooted tree:
  - its support (ultrafast bootstrap);
  - its sister group: the smallest side of a bipartition that contains the clade and other leaves (the sister group
    of a tree rooted on the larger remaining side), described by taxonomic group (config/taxon_groups.tsv), with its
    support;
  - the nearest vertebrate genes: the human genes found in the smallest such side that contains vertebrate
    sequences with a known human gene (human RefSeq proteins of the family, named Human_<gene>, and vertebrate
    predictions with their best human hit), with the support of that side and its numbers of other tunicate and
    vertebrate sequences;
  - the best human hits of its members (results/tables/stage2_selenoprotein_rbh.tsv), as gene and status;
  - the Selenoprofiles labels of its members (selenocysteine, cysteine, homologue...; "gold standard" for the
    gold-standard proteins, leaves <genome>_gold_<gene>), and the curated O. dioica genes among them: the genes of the
    gold-standard proteins, and, with --annotation (trees of stage 1 predictions), the genes of
    curation/baseline_sec_genes.tsv matched to the stage 1 predictions by their selenocysteine codon in
    results/tables/stage1_annotation.tsv, so that only predictions with a selenocysteine are named.
The best human hits of a gold-standard protein are those of the predictions it replaced (<family>.removed.tsv).
This is a summary for review: the reading is root-free but depends on the tree, and a side that holds most of the
tree, or genes of several paralogs, does not assign the clade to one gene.
Outputs: --table (one row per family, tree and appendicularian clade) and --markdown (the same, readable).
"""
import argparse, collections, csv, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import machinery_lba as lba   # taxon assignment and tree reading
from machinery_phylogeny import APPENDICULARIANS

FAMILIES = ["DIO", "GPX", "MsrA", "SELENOF", "SELENOK", "SELENOL", "SELENOO", "SELENOS", "SELENOT", "SELENOU",
            "SELENOW", "TXNRD"]
VERTEBRATE_GROUPS = {"Mammals", "Birds", "Reptiles", "Amphibians", "Ray-finned fishes", "Lobe-finned fishes",
                     "Jawless fishes", "Cartilaginous fishes"}


def split_name(leaf):
    """Genome code and prediction of a leaf named <genome>_<prediction> by machinery_phylogeny.py."""
    g, _, pred = leaf.partition("_")
    return g, pred


def gold_gene(leaf):
    """Gene of a gold-standard leaf (<genome>_gold_<gene>), else None."""
    m = re.match(r"[A-Za-z0-9]+_gold_(.+)$", leaf)
    return m.group(1) if m else None


def replaced(path):
    """gold-standard leaf -> predictions it replaced, from <family>.removed.tsv of machinery_phylogeny.py."""
    out = collections.defaultdict(list)
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"):
            m = re.fullmatch(r"gold-standard gene (\S+) used", r["reason"])
            if m: out[f"{split_name(r['sequence'])[0]}_gold_{m.group(1)}"].append(r["sequence"])
    return out


def human_gene(leaf, groups, rbh):
    """Human gene of a vertebrate leaf: Human_<gene>_<protein>, or the best human hit of a vertebrate prediction."""
    m = re.match(r"Human_([A-Za-z0-9]+)_", leaf)
    if m: return m.group(1)
    if groups.get(leaf) == "vertebrate":
        r = rbh.get(split_name(leaf))
        if r and r["best_human_gene"] not in ("-", ""): return r["best_human_gene"]
    return None


def sides(tree):
    """All sides of the bipartitions of an unrooted tree (both sides of each internal branch), with supports."""
    names = frozenset(x.name for x in tree.get_terminals())
    out = {}
    for c in tree.get_nonterminals():
        s = frozenset(x.name for x in c.get_terminals())
        for side in (s, names - s):
            if 0 < len(side) < len(names): out.setdefault(side, c.confidence)
    for x in names: out.setdefault(frozenset([x]), None)
    return names, out


def describe(leaves, tax):
    c = collections.Counter(tax[n] for n in leaves)
    return ", ".join(f"{k} {v}" for k, v in sorted(c.items(), key=lambda kv: (-kv[1], kv[0])))


def fmt(x): return "NA" if x is None else f"{x:.0f}"


def baseline_names(baseline, annotation):
    """(genome, prediction) -> curated gene, for stage 1 predictions whose selenocysteine codon is a curated one."""
    cur = {}
    for r in csv.DictReader(open(baseline), delimiter="\t"):
        cur[(r["genome"], r["contig"], int(r["sec_codon"]))] = r["gene"]
    out = {}
    for r in csv.DictReader(open(annotation), delimiter="\t"):
        if r["sec_start"] in (".", ""): continue
        for d in (0, -2, 2):
            k = (r["genome"], r["contig"], int(r["sec_start"]) + d)
            if k in cur: out[(r["genome"], r["prediction"])] = cur[k]; break
    return out


def analyse(fam, label, treefile, groups, tax, rbh, curated=None, gold=None):
    curated, gold = curated or {}, gold or {}
    tree = lba.read_tree(treefile)
    names, sd = sides(tree)
    app = {n for n in names if tax[n] == "Appendicularians"}
    app_sides = [s for s in sd if s <= app]
    clades = [s for s in app_sides if not any(s < t for t in app_sides)]
    rows = []
    for k, T in enumerate(sorted(clades, key=lambda s: (-len(s), sorted(s))), 1):
        containing = sorted((s for s in sd if T < s), key=len)
        sister = containing[0] if containing else None
        vert = next((s for s in containing if any(human_gene(n, groups, rbh) for n in s - T)), None)
        genes = collections.Counter(g for n in (vert - T if vert else ()) for g in [human_gene(n, groups, rbh)] if g)
        members = sorted(T)
        hits = collections.Counter()
        for n in members:
            for p in (gold.get(n) or ["(none)"] if gold_gene(n) else [n]):
                r = rbh.get(split_name(p))
                hits[f"{r['best_human_gene']} ({r['status']})" if r else "not in the best-hit table"] += 1
        labels = collections.Counter("gold standard" if gold_gene(n) else n.rsplit(".", 1)[-1] for n in members)
        named = sorted([f"{split_name(n)[0]}:{curated[split_name(n)]}" for n in members if split_name(n) in curated] +
                       [f"{split_name(n)[0]}:{gold_gene(n)}" for n in members if gold_gene(n)])
        rows.append({
            "family": fam, "tree": label, "clade": k, "sequences": len(T),
            "labels": ", ".join(f"{l} {v}" for l, v in sorted(labels.items(), key=lambda kv: (-kv[1], kv[0]))),
            "curated O. dioica genes": " ".join(named) or "-",
            "genomes": ",".join(sorted({split_name(n)[0] for n in T}, key=lambda g: sorted(APPENDICULARIANS).index(g)
                                       if g in APPENDICULARIANS else 99)),
            "members": " ".join(members), "clade support": fmt(sd.get(T)) if len(T) > 1 else "-",
            "sister group": describe(sister - T, tax) if sister else "-",
            "sister support": fmt(sd.get(sister)) if sister else "-",
            "nearest vertebrate genes": ", ".join(f"{g} ({v})" for g, v in sorted(genes.items())) or "-",
            "their side": f"{len(vert)} of {len(names)} leaves" if vert else "-",
            "side support": fmt(sd.get(vert)) if vert else "-",
            "other tunicates in side": sum(tax[n] == "Other tunicates" for n in vert) if vert else 0,
            "vertebrates in side": sum(tax[n] in VERTEBRATE_GROUPS for n in vert) if vert else 0,
            "best human hits": "; ".join(f"{h} x{v}" for h, v in hits.most_common()),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=FAMILIES)
    ap.add_argument("--trees", default="results/trees/selenoproteins")
    ap.add_argument("--pmsf", default="results/trees/selenoproteins/lba")
    ap.add_argument("--rbh", default="results/tables/stage2_selenoprotein_rbh.tsv")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--baseline", default="curation/baseline_sec_genes.tsv")
    ap.add_argument("--annotation", default="", help="stage 1 annotation table, to name the curated genes among stage "
                    "1 predictions (trees of the first pass: results/tables/stage1_annotation.tsv)")
    ap.add_argument("--table", default="results/tables/stage2_selenoprotein_orthology.tsv")
    ap.add_argument("--markdown", default="results/tables/stage2_selenoprotein_orthology.md")
    a = ap.parse_args()
    rbh = {}
    if os.path.exists(a.rbh):
        for r in csv.DictReader(open(a.rbh), delimiter="\t"): rbh[(r["genome"], r["prediction"])] = r
    else:
        print(f"WARNING: {a.rbh} not found; best human hits left empty", flush=True)
    taxa = lba.load_taxa(a.taxa)
    curated = baseline_names(a.baseline, a.annotation) if a.annotation and os.path.exists(a.baseline) and \
        os.path.exists(a.annotation) else {}
    rows = []
    for fam in a.families:
        gfile = f"{a.trees}/{fam}.groups.tsv"
        if not os.path.exists(f"{a.trees}/{fam}.treefile") or not os.path.exists(gfile):
            print(f"{fam}: no tree; skipped", flush=True); continue
        groups = dict(l.rstrip("\n").split("\t") for l in open(gfile))
        for label, tf in (("ML", f"{a.trees}/{fam}.treefile"), ("PMSF", f"{a.pmsf}/{fam}.pmsf.treefile")):
            if not os.path.exists(tf): continue
            leaves = {x.name for x in lba.read_tree(tf).get_terminals()}
            tax = {n: lba.taxon(n, groups, taxa) for n in leaves}
            rows += analyse(fam, label, tf, groups, tax, rbh, curated, replaced(f"{a.trees}/{fam}.removed.tsv"))
    if not rows: raise SystemExit("no tree found")
    cols = list(rows[0])
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(a.markdown, "w") as o:
        o.write("# Appendicularian selenoprotein clades: first-pass orthology\n\n"
                "Generated by stage2/selenoprotein_orthology.py; see its docstring for the reading of the trees.\n")
        fam = None
        for r in rows:
            if r["family"] != fam:
                fam = r["family"]; o.write(f"\n## {fam}\n\n")
            o.write(f"- {r['tree']} clade {r['clade']} ({r['sequences']} sequences: {r['labels']}; {r['genomes']}; "
                    f"curated: {r['curated O. dioica genes']}; support {r['clade support']}): sister group {r['sister group']} ({r['sister support']}); nearest "
                    f"vertebrate genes {r['nearest vertebrate genes']} (side of {r['their side']}, support "
                    f"{r['side support']}, {r['other tunicates in side']} other tunicate sequences); best human hits: "
                    f"{r['best human hits']}\n")
    print(open(a.markdown).read())


if __name__ == "__main__":
    main()
