#!/usr/bin/env python3
"""Display versions of the machinery trees: taxonomic groups collapsed, weak nodes as polytomies.

Two operations, both common in published trees:
  1. nodes with ultrafast bootstrap support below --min-support (default 50) are collapsed into polytomies;
  2. clades whose leaves all belong to one taxonomic group (config/taxon_groups.tsv, assigned from the genus or the
     genome code) are drawn as one tip, "group (number of sequences)", whatever their internal supports. The
     appendicularian and other tunicate sequences are never collapsed, so that their placement stays visible.
Unclassified leaves (accessions without species names) are kept as they are. Supports of at least 70 are printed.
A column gives the length of each protein (unaligned), or for a collapsed tip the median and the range of its
proteins, to show which sequences are partial. The trees and alignments themselves are not changed. Inputs: results/trees/<family>.treefile and .groups.tsv
(or another directory with --trees). Outputs: <figures>/phylogeny_<family>_collapsed.{png,pdf}.
"""
import argparse, copy, re

KEEP = {"Appendicularians", "Other tunicates"}
# Highlighted references: the related family used as outgroup, the paralogs that the appendicularian proteins could be
# confused with, and the human gene that names each family.
REFERENCE = {"pstk": {"outgroup": r"^KTI12_|^Human_KTI12", "named": r"^Human_PSTK_"},
             "SEPHS2": {"paralog": r"^Human_SEPHS1_|^sps1_|threonine", "named": r"^Human_SEPHS2_"},
             "SBP2": {"paralog": r"^Human_SECISBP2L_", "named": r"^Human_SECISBP2_"}}
BANDS = {"outgroup": ("#d5d8dc", "related family (outgroup)"), "paralog": ("#e8daef", "paralog of the named gene"),
         "named": ("#d4efdf", "clade of the human gene that names the family")}
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
    lengths = sequence_lengths(f"{a.trees}/{fam}.treefile") or sequence_lengths(f"{a.groups_dir}/{fam}.treefile")
    members = {}
    refs = REFERENCE.get(fam, {})
    cat = {c.name: next((k for k, pat in refs.items() if re.search(pat, c.name)), None) for c in tree.get_terminals()}
    # rooting, as in stage2/machinery_phylogeny.py and the analyses: outgroup for PSTK, non-animal sequences for the
    # SECIS-binding proteins, midpoint for SEPHS
    out = [c for c in tree.get_terminals() if grp[c.name] == "outgroup"]
    nonanimal = [c for c in tree.get_terminals() if grp[c.name] in
                 ("Green algae", "Red algae", "Stramenopiles", "Cryptophytes", "Amoebozoans", "Heteroloboseans")]
    def reroot(t):
        try:
            if out: t.root_with_outgroup(out if t.is_monophyletic(out) else out[0])
            elif fam == "SBP2" and nonanimal: t.root_with_outgroup(nonanimal if t.is_monophyletic(nonanimal) else nonanimal[0])
            else: t.root_at_midpoint()
        except Exception:
            t.root_at_midpoint()
    reroot_keeping_supports(tree, reroot)
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
                member_names = [t.name for t in clade.get_terminals()]
                cats = {cat.get(t.name) for t in clade.get_terminals()} - {None}
                human = sorted({re.sub(r"^Human_([A-Z0-9]+)_.*", r"\1", t.name) for t in clade.get_terminals()
                                if t.name.startswith("Human_")})
                clade.clades = []
                clade.name = (f"{g} ({n})" if g != "outgroup" else f"KTI12, outgroup ({n})") + \
                             (f", with human {', '.join(human)}" if human else "")
                cat[clade.name] = "outgroup" if "outgroup" in cats else next(iter(cats), None)
                members[clade.name] = member_names
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
        return name if "(" in name else name.split("|")[0][:40]
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(9, max(4, 0.22 * n)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: short(c.name) if c.is_terminal() else None,
               label_colors=lambda lab: COLOURS.get(next((v for k, v in grp.items() if short(k) == lab), ""), "#2c3e50"),
               branch_labels=lambda c: f"{c.confidence:.0f}" if (c.confidence is not None and c.confidence >= 70
                                                                 and not c.is_terminal()) else "")
    for t in ax.texts: t.set_fontsize(6)
    length_column(ax, tree, lengths, members)
    # background bands: outgroup, paralogs, and the clade of the human gene that names the family; the bands of the
    # named gene extend to the smallest clade that holds its human sequence and other reference tips
    from matplotlib.patches import Patch
    tips = tree.get_terminals()
    # extend the paralog and named-gene bands to the largest clade around the human tip that holds only vertebrate
    # tips (or unclassified accessions), so that the band covers the vertebrate clade of that gene
    vert = {"Mammals", "Birds", "Reptiles", "Amphibians", "Ray-finned fishes", "unclassified"}
    for tip in list(tips):
        k = cat.get(tip.name)
        if k not in ("paralog", "named") or not tip.name.startswith(grp.get(tip.name, "") or "?"):
            if k not in ("paralog", "named"): continue
        path = tree.get_path(tip); best = None
        for anc in reversed(path[:-1]):
            if all(grp.get(t.name) in vert and cat.get(t.name) in (None, k) for t in anc.get_terminals()): best = anc
            else: break
        if best is not None:
            for t in best.get_terminals(): cat[t.name] = k
    for i, tip in enumerate(tips):
        k = cat.get(tip.name)
        if k in BANDS: ax.axhspan(i + 0.5, i + 1.5, color=BANDS[k][0], zorder=0, lw=0)
    used = [k for k in BANDS if any(cat.get(t.name) == k for t in tips)]
    if used: ax.legend(handles=[Patch(color=BANDS[k][0], label=BANDS[k][1]) for k in used], loc="upper left",
                       bbox_to_anchor=(0, -0.04), ncol=len(used), fontsize=6, frameon=False)
    ax.set_yticks([]); ax.set_ylabel("")
    ax.set_title(f"{fam}: groups of one taxon collapsed; nodes with ultrafast bootstrap < {a.min_support:.0f} collapsed; "
                 f"supports >= 70 shown", fontsize=8)
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    for ext in ("png", "pdf"): fig.savefig(f"{a.figures}/phylogeny_{fam}_collapsed{a.suffix}.{ext}", dpi=200)
    plt.close(fig)
    print(fam, n, "tips after collapsing")


def reroot_keeping_supports(tree, reroot):
    """Reroot with reroot(tree), then put each support back on its bipartition: Bio.Phylo keeps supports on node
    objects, so that rerooting moves them by one branch along the path between the old and the new root."""
    names = frozenset(x.name for x in tree.get_terminals())

    def key(c):
        leaves = frozenset(x.name for x in c.get_terminals())
        return min(leaves, names - leaves, key=lambda s: (len(s), sorted(s)))
    before = {key(c): c.confidence for c in tree.get_nonterminals() if c is not tree.root}
    reroot(tree)
    for c in tree.get_nonterminals():
        c.confidence = None if c is tree.root else before.get(key(c))


def sequence_lengths(treefile):
    """Unaligned lengths of the sequences of a tree, from <family>.faa next to the tree or in the parent directory."""
    import os
    base = treefile[:-len(".treefile")]
    for faa in (f"{base}.faa", os.path.join(os.path.dirname(os.path.dirname(base)), os.path.basename(base) + ".faa")):
        if os.path.exists(faa):
            lengths, name = {}, None
            for l in open(faa):
                if l.startswith(">"): name = l[1:].strip(); lengths[name] = 0
                else: lengths[name] += len(l.strip().replace("-", ""))
            return lengths
    return {}


def length_column(ax, tree, lengths, members=None, fontsize=6):
    """Protein length to the right of each tip; for a collapsed tip, the median and the range of its sequences."""
    import statistics
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.01, 0.2, "length (aa)", transform=tr, fontsize=fontsize, fontweight="bold", va="center")
    for i, tip in enumerate(tree.get_terminals()):
        names = (members or {}).get(tip.name, [tip.name])
        ls = [lengths[n] for n in names if n in lengths]
        if not ls: continue
        text = f"{ls[0]}" if len(ls) == 1 else f"{statistics.median(ls):.0f} ({min(ls)}-{max(ls)})"
        ax.text(1.01, i + 1, text, transform=tr, fontsize=fontsize, va="center")


if __name__ == "__main__":
    main()
