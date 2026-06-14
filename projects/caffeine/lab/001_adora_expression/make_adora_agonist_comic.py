#!/usr/bin/env python
"""Comic-panel illustration: what happens when adenosine (the agonist) binds each
ADORA receptor.

Four panels (A1, A2A, A2B, A3). Each shows the receptor in the cell membrane,
adenosine docking into the extracellular pocket, the G-protein it couples to,
the effect on adenylyl cyclase, and the resulting change in cAMP. Inhibitory
receptors (Gi -> down cAMP) use a cool palette; stimulatory (Gs -> up cAMP) warm.

Biology summary (standard pharmacology; see wiki/concepts/camp-signaling.md):
    A1   Gi/o  -> inhibits AC -> down cAMP   (brain "brake"; the one caffeine blocks for alertness)
    A2A  Gs/Golf-> activates AC -> up cAMP    (striatal medium spiny neurons; classic caffeine target)
    A2B  Gs    -> activates AC -> up cAMP     (low-affinity; epithelium/vasculature, high-adenosine)
    A3   Gi/o  -> inhibits AC -> down cAMP    (immune: mast cells, macrophages, microglia)

Output: figures/adora_agonist_comic.{png,svg}
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

LAB_DIR = Path(__file__).resolve().parent
FIG_DIR = LAB_DIR / "figures"

PANELS = [
    {
        "name": "ADORA1  (A1)",
        "g": "Gαi/o",
        "direction": "down",
        "accent": "#2b6cb0",      # cool blue = inhibitory
        "soft": "#e6f0fb",
        "caption": "BRAKES ON. Gi inhibits adenylyl cyclase, so cAMP drops. "
                   "Calms neuron firing across the brain. This is the receptor "
                   "caffeine blocks to keep you alert.",
    },
    {
        "name": "ADORA2A  (A2A)",
        "g": "Gαs/olf",
        "direction": "up",
        "accent": "#dd6b20",      # warm orange = stimulatory
        "soft": "#fdebd8",
        "caption": "GAS PEDAL. Gs activates adenylyl cyclase, so cAMP rises. "
                   "Sits on striatal medium spiny neurons — the canonical "
                   "caffeine target (top A2A cell type in our data).",
    },
    {
        "name": "ADORA2B  (A2B)",
        "g": "Gαs",
        "direction": "up",
        "accent": "#b7791f",      # amber = stimulatory, low-affinity
        "soft": "#fbf2da",
        "caption": "LOW-AFFINITY BACKUP. Also Gs -> up cAMP, but only when "
                   "adenosine runs high (stress, inflammation, hypoxia). "
                   "Barrier epithelium and vasculature.",
    },
    {
        "name": "ADORA3  (A3)",
        "g": "Gαi/o",
        "direction": "down",
        "accent": "#2c7a7b",      # teal = inhibitory
        "soft": "#ddf0f0",
        "caption": "BRAKES, ON IMMUNE CELLS. Gi -> down cAMP. Concentrated in "
                   "mast cells, macrophages and microglia (top ADORA3 cells "
                   "in our data).",
    },
]


def draw_panel(ax, cfg: dict) -> None:
    accent, soft = cfg["accent"], cfg["soft"]
    up = cfg["direction"] == "up"
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Panel frame
    ax.add_patch(FancyBboxPatch((0.15, 0.15), 9.7, 9.7,
                 boxstyle="round,pad=0.02,rounding_size=0.35",
                 fill=True, facecolor="white", edgecolor=accent, linewidth=3))

    # Title banner
    ax.add_patch(FancyBboxPatch((0.5, 9.0), 9.0, 0.7,
                 boxstyle="round,pad=0.02,rounding_size=0.2",
                 facecolor=accent, edgecolor="none"))
    ax.text(5, 9.35, cfg["name"], ha="center", va="center",
            fontsize=15, fontweight="bold", color="white")

    # Extracellular / intracellular labels
    ax.text(0.75, 8.35, "outside cell", fontsize=7, style="italic", color="#666")
    ax.text(0.75, 5.05, "inside cell", fontsize=7, style="italic", color="#666")

    # Membrane bilayer band
    mem_y, mem_h = 5.7, 1.5
    ax.add_patch(mpatches.Rectangle((0.5, mem_y), 9.0, mem_h,
                 facecolor="#f4e8c1", edgecolor="#e0cf94", linewidth=1))
    # phospholipid heads (two rows of dots)
    for x in [0.8 + 0.5 * i for i in range(18)]:
        ax.add_patch(mpatches.Circle((x, mem_y + mem_h - 0.18), 0.12,
                     facecolor="#d9b94f", edgecolor="none"))
        ax.add_patch(mpatches.Circle((x, mem_y + 0.18), 0.12,
                     facecolor="#d9b94f", edgecolor="none"))
    ax.text(9.25, mem_y + mem_h / 2, "membrane", fontsize=6.5, ha="right",
            va="center", color="#8a7320", rotation=0)

    # The 7-transmembrane receptor (barrel straddling the membrane)
    rx = 4.35
    ax.add_patch(FancyBboxPatch((rx, mem_y - 0.35), 1.3, mem_h + 0.7,
                 boxstyle="round,pad=0.02,rounding_size=0.25",
                 facecolor=accent, edgecolor="black", linewidth=1.5, alpha=0.92))
    ax.text(rx + 0.65, mem_y + mem_h / 2, "7TM", ha="center", va="center",
            fontsize=8, color="white", fontweight="bold")

    # Adenosine ligand (the agonist) docking from above
    ax.add_patch(mpatches.Circle((5.0, 8.15), 0.38, facecolor="#e53e3e",
                 edgecolor="black", linewidth=1.3))
    ax.text(5.0, 8.15, "A", ha="center", va="center", fontsize=12,
            color="white", fontweight="bold")
    ax.text(4.45, 8.15, "Adenosine\n(agonist)", ha="right", va="center",
            fontsize=7.5, color="#9b2c2c", fontweight="bold")
    ax.add_patch(FancyArrowPatch((5.0, 7.72), (5.0, 7.5),
                 arrowstyle="-|>", mutation_scale=16, color="#9b2c2c", linewidth=2))

    # G-protein blob (intracellular, under receptor)
    gx, gy = 3.0, 4.6
    ax.add_patch(mpatches.Ellipse((gx, gy), 1.7, 0.95, facecolor="#9f7aea",
                 edgecolor="black", linewidth=1.3))
    ax.text(gx, gy, cfg["g"], ha="center", va="center", fontsize=9.5,
            color="white", fontweight="bold")
    # receptor -> G-protein activation arrow
    ax.add_patch(FancyArrowPatch((rx + 0.2, mem_y - 0.3), (gx + 0.55, gy + 0.45),
                 arrowstyle="-|>", mutation_scale=14, color="black", linewidth=1.6))

    # Adenylyl cyclase box
    acx, acy = 6.9, 4.6
    ax.add_patch(FancyBboxPatch((acx - 0.95, acy - 0.5), 1.9, 1.0,
                 boxstyle="round,pad=0.02,rounding_size=0.15",
                 facecolor="#edf2f7", edgecolor="#4a5568", linewidth=1.4))
    ax.text(acx, acy, "adenylyl\ncyclase", ha="center", va="center",
            fontsize=8, color="#2d3748", fontweight="bold")
    # G-protein -> AC arrow: activation (->) or inhibition (-| )
    if up:
        ax.add_patch(FancyArrowPatch((gx + 0.9, gy), (acx - 1.0, acy),
                     arrowstyle="-|>", mutation_scale=15, color="#2f855a", linewidth=2.2))
        ax.text((gx + acx) / 2, gy + 0.45, "activates", ha="center",
                fontsize=7.5, color="#2f855a", fontweight="bold")
    else:
        ax.add_patch(FancyArrowPatch((gx + 0.9, gy), (acx - 1.0, acy),
                     arrowstyle="-[", mutation_scale=14, color="#c53030", linewidth=2.2))
        ax.text((gx + acx) / 2, gy + 0.45, "inhibits", ha="center",
                fontsize=7.5, color="#c53030", fontweight="bold")

    # cAMP outcome
    out_color = "#2f855a" if up else "#c53030"
    arrow_char = "↑" if up else "↓"
    ax.text(acx, 3.05, f"cAMP {arrow_char}", ha="center", va="center",
            fontsize=20, color=out_color, fontweight="bold")
    ax.add_patch(FancyArrowPatch((acx, acy - 0.55), (acx, 3.45),
                 arrowstyle="-|>", mutation_scale=13, color=out_color, linewidth=1.8))

    # Caption box
    ax.add_patch(FancyBboxPatch((0.5, 0.55), 9.0, 1.9,
                 boxstyle="round,pad=0.03,rounding_size=0.2",
                 facecolor=soft, edgecolor=accent, linewidth=1.2))
    wrapped = textwrap.fill(cfg["caption"], width=58)
    ax.text(5, 1.5, wrapped, ha="center", va="center", fontsize=8.3,
            color="#1a202c")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 13))
    for ax, cfg in zip(axes.flat, PANELS):
        draw_panel(ax, cfg)

    fig.suptitle("When adenosine binds: the four ADORA receptors",
                 fontsize=20, fontweight="bold", y=0.985)
    fig.text(0.5, 0.012,
             "Adenosine is the natural agonist. Caffeine is the antagonist — it "
             "occupies these pockets without flipping the switch, blocking "
             "adenosine (most relevant at A1 and A2A).",
             ha="center", fontsize=9, style="italic", color="#444")
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    for ext in ("png", "svg"):
        out = FIG_DIR / f"adora_agonist_comic.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"  -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
