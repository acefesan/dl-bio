#!/usr/bin/env python3
"""Subtree HOG coloring analysis.

For each of the top N root HOGs (by protein count in the sample), generates a UMAP
plot where proteins within that root HOG are colored by their level-1 sub-HOG
descendant, and all other proteins are shown as gray background.

This reveals whether proteins in the same evolutionary family (HOG sub-tree) cluster
together in the embedding space, at a finer resolution than root HOGs.

Usage:
    python 04_subtree_hog_analysis.py \\
        --umap chapters/chapter2/runs/.../clustering/umap_coordinates.csv \\
        --data chapters/chapter2/runs/.../dataset/cafa3_annotations.feather \\
        --output-dir chapters/chapter2/lab/003_subtree_hog_coloring/figures \\
        --top-hogs 4 \\
        --min-subhog-size 5
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent


def parse_hog_level(hog_id: str, level: int = 1) -> str | None:
    """Extract HOG path truncated to the given depth level.

    E.g. 'HOG:E0801468.10osu.6761a' at level=1 → 'HOG:E0801468.10osu'
    """
    if not isinstance(hog_id, str) or not hog_id.startswith("HOG:"):
        return None
    parts = hog_id.split(".")
    return ".".join(parts[: level + 1]) if len(parts) > level else hog_id


def load_umap_coords(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    print(f"Loaded {len(df):,} UMAP coordinates from {path.name}")
    return df


def load_hog_ids(data_path: Path) -> pd.DataFrame:
    """Load EntryID → hog_id mapping from annotations feather."""
    print(f"Loading HOG IDs from {data_path.name}...")
    df = pd.read_feather(data_path, columns=["EntryID", "hog_id"])
    df = df.drop_duplicates("EntryID")
    print(f"  {len(df):,} unique proteins")
    return df[["EntryID", "hog_id"]]


def plot_subtree_hog(
    df: pd.DataFrame,
    roothog_id: float,
    output_path: Path,
    min_subhog_size: int = 5,
    hog_level: int = 1,
):
    """UMAP plot: gray background, proteins in roothog_id colored by sub-HOG."""
    in_hog = df["roothog_id"] == roothog_id
    hog_df = df[in_hog].copy()
    bg_df = df[~in_hog]

    if len(hog_df) == 0:
        print(f"  No proteins found for roothog_id={roothog_id}, skipping.")
        return

    # Parse sub-HOG at given level
    hog_df["sub_hog"] = hog_df["hog_id"].apply(lambda x: parse_hog_level(x, hog_level))

    # Group small sub-HOGs into "Other"
    subhog_counts = hog_df["sub_hog"].value_counts()
    small = subhog_counts[subhog_counts < min_subhog_size].index
    hog_df["sub_hog_label"] = hog_df["sub_hog"].where(
        ~hog_df["sub_hog"].isin(small), other="Other (small sub-HOGs)"
    )

    n_subhogs = hog_df["sub_hog_label"].nunique()
    root_id_str = str(int(roothog_id))

    fig, ax = plt.subplots(figsize=(12, 9))

    # Background
    ax.scatter(
        bg_df["umap_x"], bg_df["umap_y"],
        c="lightgray", alpha=0.2, s=8, label=f"Other proteins (n={len(bg_df):,})",
        zorder=1,
    )

    # Colored sub-HOGs
    colors = plt.cm.tab20.colors + plt.cm.tab20b.colors
    labels_order = (
        hog_df["sub_hog_label"]
        .value_counts()
        .index.tolist()
    )
    # Put "Other" last
    if "Other (small sub-HOGs)" in labels_order:
        labels_order.remove("Other (small sub-HOGs)")
        labels_order.append("Other (small sub-HOGs)")

    for i, label in enumerate(labels_order):
        group = hog_df[hog_df["sub_hog_label"] == label]
        color = "silver" if label == "Other (small sub-HOGs)" else colors[i % len(colors)]
        alpha = 0.4 if label == "Other (small sub-HOGs)" else 0.85
        ax.scatter(
            group["umap_x"], group["umap_y"],
            c=[color], alpha=alpha, s=25,
            label=f"{label.split('.')[-1] if '.' in label else label} (n={len(group)})",
            zorder=2,
        )

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(
        f"Root HOG {root_id_str}: Level-{hog_level} sub-HOG coloring\n"
        f"({len(hog_df):,} proteins in HOG, {n_subhogs} sub-HOG groups shown)",
        fontsize=13,
    )
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7, title=f"Sub-HOG (level {hog_level})")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Subtree HOG coloring on existing UMAP")
    parser.add_argument("--umap", required=True,
                        help="Path to umap_coordinates.csv from a clustering run")
    parser.add_argument("--data", required=True,
                        help="Path to cafa3_annotations.feather (for hog_id column)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for plots")
    parser.add_argument("--top-hogs", type=int, default=4,
                        help="Number of top root HOGs to plot (default: 4)")
    parser.add_argument("--min-subhog-size", type=int, default=5,
                        help="Min proteins for a sub-HOG to get its own color (default: 5)")
    parser.add_argument("--hog-level", type=int, default=1,
                        help="Sub-HOG depth level to color by (default: 1)")
    parser.add_argument("--hog-ids", type=str, default=None,
                        help="Comma-separated explicit roothog_ids to plot (overrides --top-hogs)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    umap_df = load_umap_coords(Path(args.umap))
    hog_map = load_hog_ids(Path(args.data))

    # Join hog_id onto UMAP coordinates
    df = umap_df.merge(hog_map, on="EntryID", how="left")
    print(f"Joined: {len(df):,} proteins, {df['hog_id'].notna().sum():,} with hog_id")

    # Select root HOGs to plot
    if args.hog_ids:
        selected_hogs = [float(x.strip()) for x in args.hog_ids.split(",")]
    else:
        hog_counts = df[df["roothog_id"] != 0]["roothog_id"].value_counts()
        selected_hogs = hog_counts.head(args.top_hogs).index.tolist()

    print(f"\nSelected root HOGs: {[int(h) for h in selected_hogs]}")

    # Summary stats
    summary = {}
    for roothog_id in selected_hogs:
        n = (df["roothog_id"] == roothog_id).sum()
        pct = 100 * n / len(df)
        print(f"  HOG {int(roothog_id)}: {n:,} proteins ({pct:.1f}% of sample)")
        summary[int(roothog_id)] = {"n_proteins": int(n), "pct_of_sample": round(float(pct), 2)}

    # Generate one plot per root HOG
    print(f"\nGenerating subtree HOG plots (level={args.hog_level})...")
    for roothog_id in selected_hogs:
        out_file = output_dir / f"subtree_hog_{int(roothog_id)}_level{args.hog_level}.png"
        plot_subtree_hog(df, roothog_id, out_file,
                         min_subhog_size=args.min_subhog_size,
                         hog_level=args.hog_level)

    # Also generate a combined overview: all top HOGs on one plot with one color per HOG
    print("\nGenerating combined overview plot...")
    fig, ax = plt.subplots(figsize=(14, 10))

    # Background: proteins not in any of the selected HOGs
    mask_any = df["roothog_id"].isin(selected_hogs)
    bg_df = df[~mask_any]
    ax.scatter(bg_df["umap_x"], bg_df["umap_y"],
               c="lightgray", alpha=0.15, s=8, label=f"Other (n={len(bg_df):,})", zorder=1)

    colors = plt.cm.tab10.colors
    for i, roothog_id in enumerate(selected_hogs):
        group = df[df["roothog_id"] == roothog_id]
        ax.scatter(group["umap_x"], group["umap_y"],
                   c=[colors[i % len(colors)]], alpha=0.75, s=20,
                   label=f"HOG {int(roothog_id)} (n={len(group):,})",
                   zorder=2)

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(f"Top {len(selected_hogs)} root HOGs highlighted in UMAP", fontsize=13)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=9)
    plt.tight_layout()
    overview_path = output_dir / "subtree_overview.png"
    plt.savefig(overview_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {overview_path}")

    # Save summary
    with open(output_dir / "subtree_summary.json", "w") as f:
        json.dump({"hog_level": args.hog_level, "selected_hogs": summary}, f, indent=2)

    print(f"\nDone. Output: {output_dir}")


if __name__ == "__main__":
    main()
