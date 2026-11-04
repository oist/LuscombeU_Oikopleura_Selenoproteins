#!/usr/bin/env python3
"""Simplified view of a machinery tree, for a figure panel next to the gene models.

The tree is rooted as in figures/phylogeny_collapsed.py (supports kept on their bipartitions), nodes supported by less
than --min-support become polytomies, and each leaf gets a display group: the appendicularian species (O. dioica, with
its assemblies counted), the other tunicates, the vertebrates, the outgroup (KTI12 for PSTK), the archaea, and "other" for
all other sequences (other animals, fungi, stramenopiles, alveolates, green algae); a tip of one or two sequences of
the groups "other" or "archaea" is named by their genera. Clades whose leaves all share a display group are drawn as one tip
with the number of sequences; tips of the group "other" that are siblings in the same polytomy are drawn as one tip
("unresolved"), since their branching order is not supported. The tree is drawn as a phylogram: a collapsed tip is a
wedge from its ancestor to the deepest of its sequences. Supports of at least --show-support are printed.
Standalone use writes <out>.{png,pdf}; figures/reannotation_page.py --tree-panel draws it as panel D.
"""
import argparse, collections, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "stage2"))
from phylogeny_collapsed import leaf_group, reroot_keeping_supports   # noqa: E402
from taxonomy import load_taxa                                       # noqa: E402

VERTEBRATES = {"Mammals", "Birds", "Reptiles", "Amphibians", "Ray-finned fishes", "Lobe-finned fishes", "Jawless fishes"}
SPECIES = [(r"^(OKI|OSKA|OdB3)_", "O. dioica"), (r"^Olon_", "O. longicauda"), (r"^Oalb_", "O. albicans"),
           (r"^Ovan_", "O. vanhoeffeni"), (r"^Bsty_", "B. stygius"), (r"^Mery_", "M. erythrocephalus"),
           (r"^Fbor_", "F. borealis")]
COLOURS = {"appendicularian": "#c0392b", "other tunicates": "#e67e22", "vertebrates": "#2471a3", "archaea": "#5d6d7e",
           "outgroup": "#7f8c8d", "other": "#2c3e50"}


def display_group(name, fine):
    """(display group, label) of a leaf."""
    if fine == "outgroup": return "outgroup", None
    if fine == "Archaea": return "archaea", None
    if fine == "Appendicularians":
        sp = next((s for pat, s in SPECIES if re.search(pat, name)), "appendicularian")
        return "appendicularian", sp
    if fine == "Other tunicates": return "other tunicates", None
    if fine in VERTEBRATES: return "vertebrates", None
    return "other", None


ASSEMBLY = {"OKI": "OKI2018", "OSKA": "OSKA2016", "OdB3": "OdB3"}
SINGULAR = {"fungi": "fungus", "amphioxus": "amphioxus"}
NON_ANIMALS = {"Fungi": "fungi", "Stramenopiles": "stramenopiles", "Alveolates": "alveolates",
               "Green algae": "green algae", "Archaea": "archaea"}
PLURAL = {"Flies": "insects", "Mosquitoes": "insects", "Other insects": "insects", "Stramenopiles": "stramenopiles",
          "Cephalochordates": "amphioxus", "Brachiopods": "brachiopods", "Flatworms": "flatworms"}


class Tip:
    """A drawn tip: label as (italic part, plain part), display group, leaves, support of its clade, depth of the
    deepest leaf."""
    def __init__(self, label, group, leaves, support, depth1):
        self.label, self.group, self.leaves, self.support, self.depth1 = label, group, leaves, support, depth1


def simplify(treefile, groups_file, taxa_file, min_support=50, outgroup_label="KTI12, outgroup"):
    """Nested lists (inner nodes: [support, children...]) and Tip objects, with depths from the root."""
    from Bio import Phylo
    taxa = load_taxa(taxa_file)
    groups = dict(l.rstrip("\n").split("\t")[:2] for l in open(groups_file))
    tree = Phylo.read(treefile, "newick")
    for c in tree.get_nonterminals():
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    fine = {c.name: leaf_group(c.name, groups, taxa) for c in tree.get_terminals()}
    disp = {n: display_group(n, f) for n, f in fine.items()}
    out = [c for c in tree.get_terminals() if fine[c.name] == "outgroup"]

    def reroot(t):
        if out:   # on an ingroup leaf first, so that a monophyletic outgroup is a clade
            t.root_with_outgroup(next(x for x in t.get_terminals() if x not in out))
            m = t.is_monophyletic(out)   # a clade: root on the middle of its stem (Bio.Phylo would make its base a trifurcation)
            if m: t.root_with_outgroup(m, outgroup_branch_length=(m.branch_length or 0) / 2)
            else: t.root_with_outgroup(out[0])
        else: t.root_at_midpoint()
    reroot_keeping_supports(tree, reroot)
    for c in list(tree.find_clades(order="postorder")):
        if c is tree.root or c.is_terminal(): continue
        if c.confidence is not None and c.confidence < min_support:
            try: tree.collapse(c)
            except ValueError: pass
    depth = tree.depths()

    def key(n):
        g, sp = disp[n]
        return (g, sp if g == "appendicularian" else None)

    def label(k, leaves):
        g, sp = k
        n = len(leaves)
        if g == "appendicularian":
            if sp == "O. dioica":
                asm = [ASSEMBLY.get(a_, a_) for a_ in ASSEMBLY if any(x.startswith(a_ + "_") for x in leaves)]
                return sp, f" ({', '.join(asm)})" if len(asm) > 1 else ""
            return sp, f" ({n} loci)" if n > 1 else ""
        if g == "outgroup": return "", f"{outgroup_label} ({n})"
        if g == "vertebrates":
            human = sorted({re.sub(r"^Human_([A-Za-z0-9]+)_.*", r"\1", x) for x in leaves if x.startswith("Human_")})
            return "", f"vertebrates ({n})" + (f", with human {', '.join(human)}" if human else "")
        if g in ("other", "archaea"):
            count = collections.Counter(PLURAL.get(fine[x], fine[x].lower()) for x in leaves)
            names = [SINGULAR.get(k_, k_[:-1] if k_.endswith("s") else k_) if count[k_] == 1 else k_ for k_ in count]
            if len(names) > 3:
                kinds = ["other animals" if any(fine[x] not in NON_ANIMALS for x in leaves) else None] + \
                        [lab for grp, lab in NON_ANIMALS.items() if any(fine[x] == grp for x in leaves)]
                names = [k_ for k_ in kinds if k_]
            text = " and ".join([", ".join(names[:-1]), names[-1]]) if len(names) > 1 else names[0]
            if n <= 2:   # one or two sequences: their genera, then their groups in the same order
                by = sorted({(re.split(r"[_.]", x)[0], PLURAL.get(fine[x], fine[x].lower())) for x in leaves})
                kinds = [SINGULAR.get(k_, k_[:-1] if k_.endswith("s") else k_) for _, k_ in by]
                if len(by) > 1 and len({k_ for _, k_ in by}) == 1: kinds = [by[0][1]]   # two of one group: plural
                return ", ".join(g_ for g_, _ in by), f" ({', '.join(kinds)})"
            return "", f"{text} ({n})"
        return "", f"{g} ({n})"

    def build(c):
        leaves = [t.name for t in c.get_terminals()]
        ks = {key(n) for n in leaves}
        if len(ks) == 1:
            k = ks.pop()
            d1 = max(depth[t] for t in c.get_terminals())
            return Tip(label(k, leaves), k[0], leaves, None if c is tree.root or c.is_terminal() else c.confidence, d1)
        kids = [build(ch) for ch in c.clades]
        # sibling "other" tips of a polytomy: one tip, their order being unsupported
        others = [k for k in kids if isinstance(k, Tip) and k.group == "other"]
        if len(c.clades) > 2 and len(others) > 1:
            kids = [k for k in kids if k not in others]
            leaves_o = [x for k in others for x in k.leaves]
            it, plain = label(("other", None), leaves_o)
            kids.append(Tip((it, plain + ", unresolved"), "other", leaves_o, None, max(k.depth1 for k in others)))
        return [None if c is tree.root else c.confidence, depth.get(c, 0.0)] + kids

    return build(tree.root)


def tips_of(node):
    if isinstance(node, Tip): return [node]
    return [t for ch in node[2:] for t in tips_of(ch)]


def order(node):
    """Ladderize: smaller subtrees first, the outgroup on top."""
    if isinstance(node, Tip): return node
    kids = [order(ch) for ch in node[2:]]
    kids.sort(key=lambda ch: (0 if isinstance(ch, Tip) and ch.group == "outgroup" else 1, len(tips_of(ch))))
    return node[:2] + kids


def draw(ax, root, x0, y0, width, row, fontsize=4.8, show_support=70, label_w=None):
    """Draw on ax (mm coordinates, y downwards negative) from (x0, y0); returns the y of the last row."""
    from matplotlib.patches import Polygon
    from matplotlib.textpath import TextToPath
    root = order(root)
    tips = tips_of(root)
    tp = TextToPath()

    def text_w(s, italic=False):
        from matplotlib.font_manager import FontProperties
        w, _, _ = tp.get_text_width_height_descent(s, FontProperties(size=fontsize, style="italic" if italic else "normal"),
                                                   ismath=False)
        return w / 72 * 25.4
    lab_w = label_w if label_w is not None else max(text_w(t.label[0], True) + text_w(t.label[1]) for t in tips) + 1.5
    dmax = max(t.depth1 for t in tips)
    scale = (width - lab_w - 2) / dmax
    ys = {}
    for i, t in enumerate(tips): ys[id(t)] = y0 - (i + 0.5) * row

    def place(node):
        if isinstance(node, Tip): return ys[id(node)]
        yy = [place(ch) for ch in node[2:]]
        ys[id(node)] = (min(yy) + max(yy)) / 2
        return ys[id(node)]
    place(root)
    lw = 0.5

    def x_of(d): return x0 + d * scale

    def walk(node, d_parent):
        if isinstance(node, Tip):
            y, col = ys[id(node)], COLOURS.get(node.group, "black")
            xa, xb = x_of(d_parent), x_of(node.depth1)
            if len(node.leaves) > 1:
                ax.add_patch(Polygon([(xa, y), (xb, y + row * 0.32), (xb, y - row * 0.32)], closed=True,
                                     facecolor=col, edgecolor="none", alpha=0.35))
            ax.plot([xa, xb], [y, y], color=col, lw=lw, solid_capstyle="butt")
            it, plain = node.label
            xt = xb + 0.8
            if it:
                ax.text(xt, y, it, fontsize=fontsize, va="center", color=col, style="italic")
                xt += text_w(it, True)
            if plain: ax.text(xt, y, plain, fontsize=fontsize, va="center", color=col)
            if node.support is not None and node.support >= show_support:
                ax.text(xa + 0.4, y + 0.25, f"{node.support:.0f}", fontsize=fontsize - 0.8, ha="left", va="bottom",
                        color="0.25")
            return
        sup, d = node[0], node[1]
        y, x = ys[id(node)], x_of(d)
        ax.plot([x_of(d_parent), x], [y, y], color="black", lw=lw, solid_capstyle="butt")
        yy = [ys[id(ch)] for ch in node[2:]]
        ax.plot([x, x], [min(yy), max(yy)], color="black", lw=lw)
        if sup is not None and sup >= show_support:
            ax.text(x - 0.4, y + 0.25, f"{sup:.0f}", fontsize=fontsize - 0.8, ha="right", va="bottom", color="0.25")
        for ch in node[2:]: walk(ch, d)
    walk(root, root[1])
    # scale bar: 0.5 substitutions per site
    yb = y0 - len(tips) * row - 1.5
    ax.plot([x0, x0 + 0.5 * scale], [yb, yb], color="black", lw=lw)
    ax.text(x0 + 0.5 * scale + 1, yb, "0.5 substitutions per site", fontsize=fontsize - 0.8, va="center", color="0.25")
    return yb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", default="results/trees/lba/pstk.pmsf.treefile")
    ap.add_argument("--groups", default="results/trees/pstk.groups.tsv")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--min-support", type=float, default=50)
    ap.add_argument("--out", default="figures/simplified_tree_pstk")
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    root = simplify(a.tree, a.groups, a.taxa, a.min_support)
    n = len(tips_of(root))
    W, H = 80.0, n * 2.6 + 8
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(-H, 0); ax.axis("off")
    draw(ax, root, 2, -2, W - 4, 2.6)
    for ext in ("png", "pdf"): fig.savefig(f"{a.out}.{ext}", dpi=300)
    print(n, "tips;", f"{a.out}.png")


if __name__ == "__main__":
    main()
