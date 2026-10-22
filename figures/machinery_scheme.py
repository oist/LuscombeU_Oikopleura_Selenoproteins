#!/usr/bin/env python3
"""Schematic of the selenocysteine machinery: synthesis of selenocysteine on tRNA-Sec and decoding of UGA.

Drawn from scratch with matplotlib (no template). Gene names as in Table MACHINERY. Glyph shapes follow the
Systems Biology Graphical Notation (SBGN, Process Description): metabolites as ellipses, macromolecules (proteins) as
rounded rectangles, nucleic acids as rectangles with rounded bottom corners, complexes as rectangles with cut corners;
catalysis as a line ending in an open circle. Proteins of Table MACHINERY are coloured, other proteins grey, the
associated factor TRNAU1AP is outlined with a dashed line.

Biochemistry drawn (references of manuscript/references.bib): tRNA-Sec is charged with serine by seryl-tRNA
synthetase; PSTK phosphorylates the seryl group; SEPSECS converts O-phosphoseryl-tRNA-Sec to selenocysteinyl-tRNA-Sec
with selenophosphate, made from selenide and ATP by SEPHS2 [Xu2006; Labunskyy2014]; EEFSEC carries
selenocysteinyl-tRNA-Sec to the ribosome at UGA, recruited by SECISBP2 bound to the SECIS element in the 3' UTR
[Copeland2003; Labunskyy2014], which it also binds; TRNAU1AP (SECp43) binds selenocysteinyl-tRNA-Sec and is
found in complexes with SEPSECS [SmallHoward2006].
Outputs: <out>.pdf and <out>.svg (editable, text kept as text), and <out>.png with --png.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Polygon, FancyArrowPatch, Circle, PathPatch
from matplotlib.path import Path

FS = 7.5   # font size (points)
STYLE = {   # fill, edge
    "metabolite": ("#fbe7a6", "#9a7400"),
    "rna": ("#c8e6c0", "#3d7a33"),
    "protein": ("#b9d5ee", "#2c5f8d"),
    "other": ("#e3e3e3", "#6b6b6b"),
    "complex": ("#efefef", "#6b6b6b"),
}
LINE = "#333333"


def width(text, pad=3.0):
    return 0.82 * max(len(s) for s in text.split("\n")) + pad


def label(ax, x, y, text, size=FS, **kw):
    ax.text(x, y, text, ha="center", va="center", fontsize=size, zorder=5, **kw)


def metabolite(ax, x, y, text, w=None, h=4.4):
    w = w or width(text, 4)
    f, e = STYLE["metabolite"]
    ax.add_patch(Ellipse((x, y), w, h, fc=f, ec=e, lw=0.8, zorder=3))
    label(ax, x, y, text)
    return dict(x=x, y=y, w=w, h=h)


def protein(ax, x, y, text, kind="protein", w=None, h=4.6, dashed=False, bold=True):
    w = w or width(text) * (1.1 if bold else 1)   # bold text is wider
    f, e = STYLE[kind]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=1.4",
                                fc=f, ec=e, lw=0.9, ls=(0, (3, 2)) if dashed else "-", zorder=3))
    label(ax, x, y, text, fontweight="bold" if bold else "normal")
    return dict(x=x, y=y, w=w, h=h)


def rna_path(x0, y0, w, h, r):
    """Rectangle with square top corners and rounded bottom corners (SBGN nucleic acid feature)."""
    k = 0.5523 * r   # cubic Bezier approximation of a quarter circle
    verts = [(x0, y0 + h), (x0 + w, y0 + h), (x0 + w, y0 + r),
             (x0 + w, y0 + r - k), (x0 + w - r + k, y0), (x0 + w - r, y0),
             (x0 + r, y0),
             (x0 + r - k, y0), (x0, y0 + r - k), (x0, y0 + r),
             (x0, y0 + h)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return Path(verts, codes)


def rna(ax, x, y, text, w=None, h=4.6, text_kw=None):
    w = w or width(text)
    f, e = STYLE["rna"]
    ax.add_patch(PathPatch(rna_path(x - w / 2, y - h / 2, w, h, min(1.6, h / 2)), fc=f, ec=e, lw=0.9, zorder=3))
    if text: label(ax, x, y, text, **(text_kw or {}))
    return dict(x=x, y=y, w=w, h=h)


def complex_glyph(ax, x, y, w, h, text, c=2.0):
    f, e = STYLE["complex"]
    x0, x1, y0, y1 = x - w / 2, x + w / 2, y - h / 2, y + h / 2
    pts = [(x0 + c, y0), (x1 - c, y0), (x1, y0 + c), (x1, y1 - c), (x1 - c, y1), (x0 + c, y1), (x0, y1 - c), (x0, y0 + c)]
    ax.add_patch(Polygon(pts, closed=True, fc=f, ec=e, lw=0.9, zorder=2))
    if text: label(ax, x, y1 - 2.4, text)
    return dict(x=x, y=y, w=w, h=h)


def edge(g, side):
    """Point on the bounding box of a glyph: 'top', 'bottom', 'left', 'right'."""
    return {"top": (g["x"], g["y"] + g["h"] / 2), "bottom": (g["x"], g["y"] - g["h"] / 2),
            "left": (g["x"] - g["w"] / 2, g["y"]), "right": (g["x"] + g["w"] / 2, g["y"])}[side]


def arrow(ax, p, q, rad=0.0, lw=1.0, ls="-", color=LINE, z=1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>,head_length=3.2,head_width=1.8", mutation_scale=1,
                                 connectionstyle=f"arc3,rad={rad}", lw=lw, ls=ls, color=color, zorder=z,
                                 shrinkA=0, shrinkB=0))


def catalysis(ax, p, q, r=0.75):
    """Line from the enzyme to the process, ending in an open circle (SBGN catalysis)."""
    import math
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = math.hypot(dx, dy)
    end = (q[0] - dx / d * r, q[1] - dy / d * r)
    ax.plot([p[0], end[0]], [p[1], end[1]], color=LINE, lw=0.9, zorder=1, solid_capstyle="butt")
    ax.add_patch(Circle(q, r, fc="white", ec=LINE, lw=0.9, zorder=4))


def binding(ax, p, q):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=LINE, lw=0.9, ls=(0, (2.5, 1.8)), zorder=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="figures/machinery_scheme")
    ap.add_argument("--png", action="store_true", help="also write a PNG (for checking; not committed)")
    a = ap.parse_args()
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,   # text kept as text in the SVG
                         "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"]})

    X0, X1, Y0, Y1 = 0, 100, -17, 74
    fig = plt.figure(figsize=(7.2, 7.2 * (Y1 - Y0) / (X1 - X0)))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal"); ax.axis("off")

    # panels
    for y0, y1, title in ((28.5, 73, "Synthesis of selenocysteine on its tRNA"),
                          (-0.5, 27.5, "Decoding of UGA as selenocysteine")):
        ax.add_patch(FancyBboxPatch((1, y0), 98, y1 - y0, boxstyle="round,pad=0,rounding_size=1.5",
                                    fc="#f7f7f7", ec="none", zorder=0))
        ax.text(2.5, y1 - 2.2, title, fontsize=FS + 0.5, fontweight="bold", va="center", zorder=5)

    # selenium branch
    yse = 60
    selenide = metabolite(ax, 44, yse, "selenide")
    sephos = metabolite(ax, 69, yse, "selenophosphate")
    sephs2 = protein(ax, 56.5, 67, "SEPHS2")
    arrow(ax, edge(selenide, "right"), edge(sephos, "left"))
    catalysis(ax, edge(sephs2, "bottom"), (56.5, yse))
    label(ax, 56.5, yse - 2.4, "ATP → AMP + Pi", size=FS - 1.5, color="#555555")

    # tRNA states
    yt = 46
    t0 = rna(ax, 9, yt, "tRNA-Sec")
    t1 = rna(ax, 31, yt, "Ser-tRNA-Sec")
    t2 = rna(ax, 54, yt, "Sep-tRNA-Sec")
    t3 = rna(ax, 83, yt, "Sec-tRNA-Sec")
    ye = 36.5
    steps = [(t0, t1, protein(ax, 20, ye, "SerRS", kind="other", bold=False)),
             (t1, t2, protein(ax, 42.5, ye, "PSTK")),
             (t2, t3, protein(ax, 68.5, ye, "SEPSECS"))]
    for s, t, enz in steps:
        p, q = edge(s, "right"), edge(t, "left")
        arrow(ax, p, q)
        catalysis(ax, edge(enz, "top"), (enz["x"], yt))
    serine = metabolite(ax, 20, 55, "serine")
    arrow(ax, edge(serine, "bottom"), (20, yt + 0.4))
    ax.text(21, 50.2, "ATP → AMP + PPi", fontsize=FS - 1.5, color="#555555", va="center", zorder=5)
    label(ax, 42.5, yt + 2.6, "ATP → ADP", size=FS - 1.5, color="#555555")
    arrow(ax, edge(sephos, "bottom"), (68.5, yt + 0.4))
    secp43 = protein(ax, 85, 36.5, "TRNAU1AP", dashed=True)
    ax.text(85, 32.6, "associated factor", ha="center", va="center", fontsize=FS - 1.5, color="#555555", zorder=5)
    binding(ax, edge(secp43, "left"), edge(steps[2][2], "right"))
    binding(ax, edge(secp43, "top"), (85, yt - 2.3))

    # decoding
    ym = 5.5
    mrna = rna(ax, 50, ym, "", w=92, h=2.6)
    label(ax, 2.6, ym, "5′", size=FS - 1)
    label(ax, 97.6, ym, "3′", size=FS - 1)
    ax.text(8, ym - 3.3, "selenoprotein mRNA", fontsize=FS - 1, va="center", zorder=5)
    # SECIS: stem-loop in the 3' UTR
    xs = 84
    stem = Path([(xs - 1.2, ym + 1.3), (xs - 1.2, ym + 9.5), (xs - 1.2, ym + 12.5), (xs + 1.2, ym + 12.5),
                 (xs + 1.2, ym + 9.5), (xs + 1.2, ym + 1.3)],
                [Path.MOVETO, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO])
    ax.add_patch(PathPatch(stem, fc="none", ec=STYLE["rna"][1], lw=1.4, zorder=3))
    for yy in (ym + 3, ym + 5, ym + 7, ym + 9):
        ax.plot([xs - 1.2, xs + 1.2], [yy, yy], color=STYLE["rna"][1], lw=0.6, zorder=3)
    label(ax, xs + 4.8, ym + 7, "SECIS", size=FS - 0.5, color=STYLE["rna"][1], fontweight="bold")
    ax.text(xs, ym - 3.3, "3′ UTR", fontsize=FS - 1, va="center", ha="center", zorder=5)
    sbp = protein(ax, xs, 22.5, "SECISBP2/SECISBP2L")
    ax.plot([xs, xs], [ym + 12.6, edge(sbp, "bottom")[1]], color=LINE, lw=0.9, ls=(0, (2.5, 1.8)), zorder=1)
    # ribosome on UGA
    xr = 40
    ribo = complex_glyph(ax, xr, 10.6, 17, 14.5, "ribosome")
    ax.add_patch(FancyBboxPatch((xr - 2.2, ym - 1.3), 4.4, 2.6, boxstyle="square,pad=0", fc="#fdfdfd",
                                ec=STYLE["rna"][1], lw=0.9, zorder=4))
    label(ax, xr, ym, "UGA", size=FS - 1.5, fontweight="bold")
    ax.text(xr, ym - 3.3, "coding sequence", fontsize=FS - 1, va="center", ha="center", zorder=5)
    # EEFSEC carrying Sec-tRNA-Sec
    eef = protein(ax, 62, 21, "EEFSEC")
    ax.text(62, 17.3, "+ Sec-tRNA-Sec", ha="center", va="center", fontsize=FS - 1.5, color="#555555", zorder=5)
    arrow(ax, edge(t3, "bottom"), (64.5, edge(eef, "top")[1]), rad=-0.15)
    arrow(ax, edge(eef, "left"), (xr + 8.5, 14), rad=0.2)
    binding(ax, edge(sbp, "left"), edge(eef, "right"))
    binding(ax, (xs - 7, 20.2), (xr + 8.5, 9))
    # product
    sel = protein(ax, 15, 18, "selenoprotein", kind="other", bold=False)
    ax.text(15, 14.3, "Sec at the UGA position", ha="center", va="center", fontsize=FS - 1.5, color="#555555",
            zorder=5)
    arrow(ax, (xr - 8.5, 15), edge(sel, "right"))

    # legend
    yl = -6.5
    def key(x, y, text): ax.text(x, y, text, fontsize=FS - 1, va="center")
    metabolite(ax, 6, yl, "", w=7, h=3.6); key(10.5, yl, "metabolite")
    rna(ax, 27, yl, "", w=7, h=3.6); key(31.5, yl, "RNA")
    protein(ax, 43, yl, "", w=7, h=3.6); key(47.5, yl, "protein in Table MACHINERY")
    protein(ax, 75, yl, "", kind="other", w=7, h=3.6); key(79.5, yl, "other protein")
    yl2 = -12.5
    protein(ax, 6, yl2, "", w=7, h=3.6, dashed=True); key(10.5, yl2, "associated factor")
    complex_glyph(ax, 30, yl2, 7, 3.6, "", c=1.0); key(34.5, yl2, "complex")
    arrow(ax, (45, yl2), (51, yl2)); key(52.5, yl2, "conversion, transfer")
    catalysis(ax, (69, yl2), (74.5, yl2)); key(76, yl2, "catalysis")
    binding(ax, (86, yl2), (90.5, yl2)); key(91.5, yl2, "binding")

    for ext in ("pdf", "svg") + (("png",) if a.png else ()):
        fig.savefig(f"{a.out}.{ext}", dpi=200)
    print("written", a.out)


if __name__ == "__main__":
    main()
