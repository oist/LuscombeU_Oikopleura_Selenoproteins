#!/usr/bin/env python3
"""Figure LOOP: the loop from the pipeline run to the gold standard and the final settings (schematic; no data).

Steps drawn clockwise from the top; the exit arrow leads to the final run, whose results the manuscript reports.
Outputs: figures/iteration_loop.{svg,png,pdf}.
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

STEPS = [
    ("Pipeline run", "nine appendicularian genomes"),
    ("Comparison with the baseline", "curated O. dioica genes"),
    ("Targeted searches", "extended profiles, tblastn absence tests"),
    ("Gene models from transcripts", "gold standard, annotation fixes"),
    ("SECIS model", "trained on gold-standard 3′ regions"),
    ("UGA/SECIS screens", "transcripts, genomes, comparative"),
    ("Candidates", "rules fixed in advance; baseline update"),
    ("Appendicularian profiles", "tested on held-out assemblies"),
]
INK, ACCENT, MUTED, FILL = "#1E2A33", "#2B5F9E", "#56616B", "#F3F1EC"


def main():
    fig, ax = plt.subplots(figsize=(9.2, 7.6))
    ax.set_xlim(-6.2, 6.2); ax.set_ylim(-5.0, 5.6); ax.set_aspect("equal"); ax.axis("off")
    n, rx, ry = len(STEPS), 4.2, 3.9
    pos = [(rx * math.sin(2 * math.pi * k / n), ry * math.cos(2 * math.pi * k / n)) for k in range(n)]
    w, h = 3.0, 1.0
    boxes = []
    for (x, y), (title, sub) in zip(pos, STEPS):
        main_box = title == "Pipeline run"
        b = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.18", fc=FILL,
                           ec=ACCENT if main_box else MUTED, lw=1.6 if main_box else 1)
        ax.add_patch(b); boxes.append(b)
        ax.text(x, y + 0.17, title, ha="center", va="center", fontsize=9.5, fontweight="bold", color=INK)
        ax.text(x, y - 0.21, sub, ha="center", va="center", fontsize=7.2, color=MUTED)
    for k in range(n):
        ax.add_patch(FancyArrowPatch(pos[k], pos[(k + 1) % n], patchA=boxes[k], patchB=boxes[(k + 1) % n],
                                     arrowstyle="-|>", mutation_scale=12, color=MUTED, lw=1, shrinkA=3, shrinkB=3,
                                     connectionstyle="arc3,rad=-0.15"))
    ax.text(0, 0.35, "until a round adds\nno O. dioica gene", ha="center", va="center", fontsize=10, color=INK)
    ax.text(0, -0.6, "each round recorded in\nthe Git repository", ha="center", va="center", fontsize=8.5, color=MUTED)
    fb = FancyBboxPatch((3.3, 4.55), 2.6, 0.9, boxstyle="round,pad=0.02,rounding_size=0.18", fc="#E6ECF0", ec=ACCENT,
                        lw=1.6)
    ax.add_patch(fb)
    ax.text(4.6, 5.17, "Final run", ha="center", va="center", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(4.6, 4.8, "census, orthology, operons", ha="center", va="center", fontsize=7.2, color=MUTED)
    ax.add_patch(FancyArrowPatch(pos[0], (4.6, 5.0), patchA=boxes[0], patchB=fb, arrowstyle="-|>", mutation_scale=12,
                                 color=ACCENT, lw=1.4, shrinkA=3, shrinkB=3))
    fig.tight_layout()
    for ext in ("svg", "png", "pdf"):
        fig.savefig(f"figures/iteration_loop.{ext}", dpi=200)


if __name__ == "__main__":
    main()
