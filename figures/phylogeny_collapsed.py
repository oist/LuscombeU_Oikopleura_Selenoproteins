#!/usr/bin/env python3
"""Display versions of the machinery trees: taxonomic groups collapsed, weak nodes as polytomies.

Two operations, both common in published trees:
  1. nodes with ultrafast bootstrap support below --min-support (default 50) are collapsed into polytomies;
  2. clades whose leaves all belong to one taxonomic group (config/taxon_groups.tsv, assigned from the genus or the
     genome code) are drawn as one tip, "group (number of sequences)", whatever their internal supports. The
     appendicularian and other tunicate sequences are never collapsed, so that their placement stays visible.
Unclassified leaves (accessions without species names) are kept as they are. Supports of at least 70 are printed.
The trees and alignments themselves are not changed. Inputs: results/trees/<family>.treefile and .groups.tsv
(or another directory with --trees). Outputs: <figures>/phylogeny_<family>_collapsed.{png,pdf}.
"""
import argparse, copy, re

KEEP = {"Appendicularians", "Other tunicates"}
COLOURS = {"Appendicularians": "#c0392b", "Other tunicates": "#e67e22", "Mammals": "#2471a3", "Birds": "#2471a3",
           "Reptiles": "#2471a3", "Amphibians": "#2471a3", "Ray-finned fishes": "#2471a3", "outgroup": "#7f8c8d"}


def leaf_group(name, groups, taxa):
    g = groups.get(name, "profile")
    if g == "outgroup": return "outgroup"
    if g == "curated": return "Appendicularians"
    first = name.split("|")[0]
    first = re.sub(r"^(sps1|sps2|sps)_", "", first)
    if first.startswith("Human_"): return "Mammals"
    return taxa.get(first.split("_")[0], "unclassified")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=["pstk", "SEPHS2", "SBP2"])
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--groups-dir", default="results/trees")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--min-support", type=float, default=50)
    ap.add_argument("--figures", default="figures")
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    taxa = {}
    for l in open(a.taxa):
        if l.startswith("#") or l.startswith("genus\t"): continue
        g, grp = l.rstrip("\n").split("\t"); taxa[g] = grp
    for fam in a.families:
        draw(fam, a, taxa)


def draw(fam, a, taxa):
    from Bio import Phylo
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = dict(l.rstrip("\n").split("\t") for l in open(f"{a.groups_dir}/{fam}.groups.tsv"))
    tree = Phylo.read(f"{a.trees}/{fam}.treefile", "newick")
    for c in tree.get_nonterminals():
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    grp = {c.name: leaf_group(c.name, groups, taxa) for c in tree.get_terminals()}
    # rooting, as in stage2/machinery_phylogeny.py and the analyses: outgroup for PSTK, non-animal sequences for the
    # SECIS-binding proteins, midpoint for SEPHS
    out = [c for c in tree.get_terminals() if grp[c.name] == "outgroup"]
    nonanimal = [c for c in tree.get_terminals() if grp[c.name] in
                 ("Green algae", "Red algae", "Stramenopiles", "Cryptophytes", "Amoebozoans", "Heteroloboseans")]
    try:
        if out: tree.root_with_outgroup(out if tree.is_monophyletic(out) else out[0])
        elif fam == "SBP2" and nonanimal: tree.root_with_outgroup(nonanimal if tree.is_monophyletic(nonanimal) else nonanimal[0])
        else: tree.root_at_midpoint()
    except Exception:
        tree.root_at_midpoint()
    # 1. weak nodes become polytomies
    for c in list(tree.find_clades(order="postorder")):
        if c is tree.root or c.is_terminal(): continue
        if c.confidence is not None and c.confidence < a.min_support:
            try: tree.collapse(c)
            except ValueError: pass
    # 2. single-group clades become one tip
    def collapse_groups(clade):
        if clade.is_terminal(): return
        names = {grp[t.name] for t in clade.get_terminals()}
        if len(names) == 1 and len(clade.get_terminals()) > 1:
            g = names.pop()
            if g not in KEEP and g != "unclassified":
                n = len(clade.get_terminals())
                depth = max(clade.distance(t) for t in clade.get_terminals())
                clade.clades = []
                clade.name = f"{g} ({n})" if g != "outgroup" else f"KTI12, outgroup ({n})"
                clade.branch_length = (clade.branch_length or 0) + depth / 2
                grp[clade.name] = g
                return
        for ch in clade.clades: collapse_groups(ch)
    collapse_groups(tree.root)
    tree.ladderize()

    def short(name):
        g = grp.get(name, "")
        if g in KEEP: return re.sub(r"\.(homologue|selenocysteine|glycine|alanine|serine|unaligned|pseudo)$", r" (\1)",
                                    name.replace("|", " ").replace("_curated", " curated"))
        return name if "(" in name else name.split("|")[0]
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(8, max(4, 0.22 * n)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: short(c.name) if c.is_terminal() else None,
               label_colors=lambda lab: COLOURS.get(next((v for k, v in grp.items() if short(k) == lab), ""), "#2c3e50"),
               branch_labels=lambda c: f"{c.confidence:.0f}" if (c.confidence is not None and c.confidence >= 70
                                                                 and not c.is_terminal()) else "")
    for t in ax.texts: t.set_fontsize(6)
    ax.set_yticks([]); ax.set_ylabel("")
    ax.set_title(f"{fam}: groups of one taxon collapsed; nodes with ultrafast bootstrap < {a.min_support:.0f} collapsed; "
                 f"supports >= 70 shown", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(f"{a.figures}/phylogeny_{fam}_collapsed{a.suffix}.{ext}", dpi=200)
    plt.close(fig)
    print(fam, n, "tips after collapsing")


if __name__ == "__main__":
    main()
