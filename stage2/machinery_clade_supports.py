#!/usr/bin/env python3
"""Supports of the clades that the text cites, read from every machinery tree.

The long-branch summary (stage2/machinery_lba.py) records three clades per family; the text also cites other clades
(paralogous forms, the position of the appendicularian clade, vertebrate paralogs). This script reads them from all
trees of a family: maximum likelihood (results/trees/<family>.treefile), C20 (results/trees/C20/), PMSF and site
removal (results/trees/lba/), the tree without the outgroup long branches (lba/<family>.taxa), the PSTK tree with one
leaf per appendicularian species (results/trees/pstk_one_leaf/), the PSTK tree without its outgroup
(results/trees/lba_no_outgroup/, clades only, as it has no root), and the trees of the two parts of the alignment
(results/trees/domains/, selected model and PMSF). Trees missing on disk are skipped, and so are entries whose required leaves (for example the fungi, in the tree
without them) are absent.
Three kinds of entries:
  clade     support of the bipartition that separates exactly these leaves (as present in the tree), or "no";
  mrca      smallest clade of the rooted tree that holds these leaves: its size, its support and its composition;
  position  for leaves that form a clade in the rooted tree: support of its parent node and composition of its
            sister group (the parent minus the clade).
Trees are rooted as in figures/phylogeny_collapsed.py (KTI12 for PSTK, prokaryotic proteins for SEPHS, non-animal
proteins for the SECIS-binding proteins); supports are those of the bipartitions, so rerooting does not move them.
Compositions count leaves by taxonomic group (as in stage2/machinery_lba.py).
Output (--out, TSV): family, tree, entry, kind, leaves, support, composition.
"""
import argparse, collections, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from machinery_lba import read_tree, taxon   # noqa: E402
from taxonomy import load_taxa, ANIMALS      # noqa: E402

APP_SEC = r"^(OdB3|OSKA|OKI)_gold_SEPHS2$|^Oikopleura_dioica\.SPS2|^(Oalb|Ovan|Olon)_SEPHS2\.\d+\.selenocysteine"
FUNGI = r"Smittium|Zancudomyces"


def T(*groups):
    return lambda n, tax: tax[n] in groups


def R(rx):
    return lambda n, tax: re.search(rx, n) is not None


def OR(*fs): return lambda n, tax: any(f(n, tax) for f in fs)


def AND(*fs): return lambda n, tax: all(f(n, tax) for f in fs)


def NOT(f): return lambda n, tax: not f(n, tax)


ENTRIES = {
    "pstk": [
        ("appendicularians", "clade", T("Appendicularians")),
        ("other tunicates", "clade", T("Other tunicates")),
        ("appendicularians and other tunicates", "clade", T("Appendicularians", "Other tunicates")),
        ("PSTK other than appendicularians", "clade", AND(NOT(T("outgroup")), NOT(T("Appendicularians"))), [T("outgroup")]),
        ("Smittium and Zancudomyces", "clade", R(FUNGI)),
        ("appendicularians with Smittium and Zancudomyces", "clade", OR(T("Appendicularians"), R(FUNGI)), [R(FUNGI)]),
        ("appendicularians", "position", T("Appendicularians")),
        ("Smittium and Zancudomyces", "position", R(FUNGI)),
    ],
    "SEPHS2": [
        ("appendicularian selenocysteine SEPHS2", "clade", R(APP_SEC)),
        ("appendicularian selenocysteine SEPHS2 and F. borealis", "clade", OR(R(APP_SEC), R("^Fbor_")), [R("^Fbor_")]),
        ("appendicularian second copy", "clade", R(r"gold_SEPHS_b|^Olon_SEPHS2\.5\.")),
        ("appendicularian second copy and tunicate glycine forms", "clade",
         R(r"gold_SEPHS_b|^Olon_SEPHS2\.5\.|glycine|Botryllus_schlosseri\.SPS1")),
        ("tunicate SEPHS except the appendicularian selenocysteine clade", "clade",
         AND(T("Other tunicates", "Appendicularians"), NOT(R(APP_SEC)), NOT(R("^Fbor_")))),
        ("appendicularian second copy", "position", R(r"gold_SEPHS_b|^Olon_SEPHS2\.5\.")),
        ("appendicularian selenocysteine SEPHS2", "position", R(APP_SEC)),
        ("appendicularian selenocysteine SEPHS2 and F. borealis", "position", OR(R(APP_SEC), R("^Fbor_")), [R("^Fbor_")]),
    ],
    "SBP2": [
        ("oikopleurid proteins", "clade", R("_gold_SBP2")),
        ("oikopleurid proteins and other tunicates", "clade", OR(R("_gold_SBP2"), T("Other tunicates"))),
        ("appendicularians (with F. borealis) and other tunicates", "clade", T("Appendicularians", "Other tunicates"),
         [R("^Fbor_")]),
        ("tunicates and human SECISBP2L", "mrca", OR(T("Other tunicates"), R("^Human_SECISBP2L"))),
        ("tunicates and amphioxus", "mrca", OR(T("Other tunicates"), R("^(Bflo|Blan)_"))),
        ("human SECISBP2 and SECISBP2L", "mrca", R("^Human_SECISBP2")),
        ("human SECISBP2L, amphioxus and echinoderm", "mrca", OR(R("^Human_SECISBP2L"), R("^(Bflo|Blan)_"), T("Echinoderms"))),
        ("F. borealis", "position", R("^Fbor_")),
    ],
}


def trees(fam, root):
    out = [("maximum likelihood", f"{root}/{fam}.treefile"), ("C20", f"{root}/C20/{fam}.treefile"),
           ("PMSF", f"{root}/lba/{fam}.pmsf.treefile")]
    out += [(f"fastest {k} % of sites removed", f"{root}/lba/{fam}.slowfast{k}.treefile") for k in (10, 20, 30, 40)]
    out += [("without Smittium and Zancudomyces (machinery_lba.py --exclude)", f"{root}/lba/{fam}.taxa.treefile"),
            ("without the outgroup (make pstk-no-outgroup)", f"{root}/lba_no_outgroup/{fam}.taxa.treefile")]
    out += [("one leaf per species", f"{root}/{fam}_one_leaf/{fam}.treefile"),
            ("one leaf per species, C20", f"{root}/{fam}_one_leaf/C20/{fam}.treefile")]
    for p in sorted(os.listdir(f"{root}/domains")) if os.path.isdir(f"{root}/domains") else []:
        m = re.fullmatch(rf"{re.escape(fam)}\.(\w+?)(\.pmsf)?\.treefile", p)
        if m and m.group(1) not in ("full",) and not re.search(r"\.(au|constraint\d+)\.", p):
            out.append((f"domain {m.group(1).replace('_', ' ')}" + (", PMSF" if m.group(2) else ""), f"{root}/domains/{p}"))
    return [(k, v) for k, v in out if os.path.exists(v)]


def bip(t, side):
    names = {x.name for x in t.get_terminals()}
    if not side or side == names: return "NA"
    for c in t.get_nonterminals():
        L = {x.name for x in c.get_terminals()}
        if L == side or L == names - side:
            return f"{c.confidence:.0f}" if c.confidence is not None else "root"
    return "no"


def composition(leaves, tax):
    c = collections.Counter(tax[n] for n in leaves)
    return ", ".join(f"{k} {v}" for k, v in sorted(c.items(), key=lambda x: (-x[1], x[0])))


def rooted(p, fam, tax):
    t = read_tree(p)
    leaves = t.get_terminals()
    if fam == "pstk": out = [x for x in leaves if tax[x.name] == "outgroup"]
    elif fam == "SEPHS2": out = [x for x in leaves if tax[x.name] in ("Bacteria", "Archaea")]
    else: out = [x for x in leaves if tax[x.name] not in ANIMALS | {"Appendicularians", "Other tunicates", "unclassified",
                                                                       "Choanoflagellates", "Other holozoans"}]
    if out: t.root_with_outgroup(out if t.is_monophyletic(out) else out[0])
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=["pstk", "SEPHS2", "SBP2"])
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--out", default="results/tables/stage2_machinery_clade_supports.tsv")
    a = ap.parse_args()
    taxa = load_taxa(a.taxa)
    rows = []
    for fam in a.families:
        groups = {}
        for l in open(f"{a.trees}/{fam}.groups.tsv"):
            f = l.rstrip("\n").split("\t"); groups[f[0]] = f[1]
        for label, path in trees(fam, a.trees):
            t = read_tree(path)
            names = {x.name for x in t.get_terminals()}
            tax = {n: taxon(n, groups, taxa) for n in names}
            rt = rooted(path, fam, tax)
            unrooted = fam == "pstk" and not any(tax[n] == "outgroup" for n in names)   # no root: no position
            for entry, kind, sel, *req in ENTRIES.get(fam, []):
                if unrooted and kind != "clade": continue
                S = {n for n in names if sel(n, tax)}
                if not S or any(not any(f(n, tax) for n in names) for f in (req[0] if req else [])): continue
                row = dict(family=fam, tree=label, entry=entry, kind=kind, leaves=len(S), support="", composition="")
                if kind == "clade":
                    row["support"] = bip(t, S)
                else:
                    m = rt.common_ancestor([x for x in rt.get_terminals() if x.name in S])
                    M = {x.name for x in m.get_terminals()}
                    if kind == "mrca":
                        row.update(leaves=len(M), support=bip(t, M), composition=composition(M, tax))
                    elif M != S:
                        row.update(support="not a clade")
                    else:
                        path_ = rt.get_path(m)
                        parent = path_[-2] if len(path_) > 1 else rt.root
                        P = {x.name for x in parent.get_terminals()}
                        row.update(support=bip(t, P), composition="sister: " + composition(P - S, tax))
                rows.append(row)
    cols = ["family", "tree", "entry", "kind", "leaves", "support", "composition"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows)} rows written to {a.out}")


if __name__ == "__main__":
    main()
