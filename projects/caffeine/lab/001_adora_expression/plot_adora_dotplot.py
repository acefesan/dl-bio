#!/usr/bin/env python
"""Dotplot ADORA receptor expression by Tabula Sapiens cell type."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
DEFAULT_OUT = LAB_DIR / "figures" / "tabula_sapiens_adora_dotplot_cell_type.png"
DEFAULT_TABLE = LAB_DIR / "figures" / "tabula_sapiens_adora_dotplot_cell_type.csv"
GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_categorical(h5ad_path: Path, obs_column: str) -> tuple[np.ndarray, list[str]]:
    with h5py.File(h5ad_path, "r") as f:
        node = f["obs"][obs_column]
        if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
            raise ValueError(f"obs/{obs_column!r} is not an AnnData categorical column")
        return node["codes"][:], decode_values(node["categories"][:])


def summarize_by_group(
    expression: np.ndarray,
    codes: np.ndarray,
    categories: list[str],
    min_cells: int,
) -> pd.DataFrame:
    rows = []
    for code, group in enumerate(categories):
        mask = codes == code
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        group_expr = expression[mask]
        pct_expr = (group_expr > 0).mean(axis=0) * 100.0
        mean_expr = group_expr.mean(axis=0)
        for i, gene in enumerate(GENES):
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "n_cells": n_cells,
                    "pct_expressing": float(pct_expr[i]),
                    "mean_expression": float(mean_expr[i]),
                }
            )
    return pd.DataFrame(rows)


def select_groups(summary: pd.DataFrame, top_n: int) -> list[str]:
    scores = (
        summary.assign(score=summary["pct_expressing"] * np.log1p(summary["mean_expression"]))
        .groupby("group", as_index=True)
        .agg(max_pct_expressing=("pct_expressing", "max"), max_mean_expression=("mean_expression", "max"), max_score=("score", "max"))
        .sort_values(["max_score", "max_pct_expressing", "max_mean_expression"], ascending=False)
    )
    return scores.head(top_n).index[::-1].tolist()


def display_label(value: str) -> str:
    return value.replace("_", " ")


def plot_dotplot(summary: pd.DataFrame, groups: list[str], output_path: Path, title: str, ylabel: str, dpi: int) -> None:
    plot_df = summary[summary["group"].isin(groups)].copy()
    group_to_y = {group: i for i, group in enumerate(groups)}
    gene_to_x = {gene: i for i, gene in enumerate(GENES)}
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["group"].map(group_to_y)

    fig_height = max(6.0, 0.28 * len(groups) + 2.0)
    fig, ax = plt.subplots(figsize=(9.8, fig_height))
    sizes = 8 + plot_df["pct_expressing"].to_numpy() * 5.0
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=sizes,
        c=plot_df["mean_expression"],
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )

    ax.set_xticks(range(len(GENES)), GENES)
    ax.tick_params(axis="x", labelrotation=25)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.set_yticks(range(len(groups)), [display_label(group) for group in groups])
    ax.set_xlim(-0.6, len(GENES) - 0.4)
    ax.set_ylim(-0.8, len(groups) - 0.2)
    ax.grid(axis="both", color="#e6e6e6", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Adenosine receptor gene")
    ax.set_ylabel(ylabel)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.03)
    cbar.set_label("Mean expression")

    legend_pcts = [1, 5, 15, 30]
    handles = [
        ax.scatter([], [], s=8 + pct * 5.0, facecolors="none", edgecolors="#333333", linewidths=0.5)
        for pct in legend_pcts
    ]
    ax.legend(
        handles,
        [f"{pct}% expressing" for pct in legend_pcts],
        title="Dot size",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.22, 0.14),
    )

    fig.tight_layout(rect=(0, 0, 0.86, 1))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--expression-cache", type=Path, default=DEFAULT_EXPR_CACHE)
    parser.add_argument("--group-by", default="cell_type")
    parser.add_argument("--min-cells", type=int, default=50)
    parser.add_argument("--top-n", type=int, default=35)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cached = np.load(args.expression_cache, allow_pickle=False)
    expression = cached["expression"]
    codes, categories = read_categorical(args.h5ad, args.group_by)
    if expression.shape[0] != codes.shape[0]:
        raise ValueError(f"Expression rows ({expression.shape[0]}) do not match obs rows ({codes.shape[0]})")

    summary = summarize_by_group(expression, codes, categories, args.min_cells)
    groups = select_groups(summary, args.top_n)

    args.table.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_values(["group", "gene"]).to_csv(args.table, index=False)
    print(f"Wrote {args.table}")

    group_label = display_label(args.group_by)
    title = f"ADORA receptor expression by {group_label} (top {len(groups)} ADORA-enriched groups)"
    plot_dotplot(summary, groups, args.out, title, group_label.title(), args.dpi)


if __name__ == "__main__":
    main()
