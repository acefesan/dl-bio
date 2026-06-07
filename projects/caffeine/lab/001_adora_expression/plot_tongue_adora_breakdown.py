#!/usr/bin/env python
"""Break down ADORA expression within Tabula Sapiens tongue cells."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
OUT = LAB_DIR / "figures" / "tabula_sapiens_tongue_adora_cell_type_breakdown.png"
TABLE = LAB_DIR / "figures" / "tabula_sapiens_tongue_adora_cell_type_breakdown.csv"
GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")


def decode(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_categorical(f: h5py.File, column: str) -> np.ndarray:
    node = f[f"obs/{column}"]
    categories = np.array(decode(node["categories"][:]), dtype=object)
    codes = node["codes"][:]
    values = np.full(codes.shape, "missing", dtype=object)
    ok = codes >= 0
    values[ok] = categories[codes[ok]]
    return values


def summarize(expression: np.ndarray, tissue: np.ndarray, cell_type: np.ndarray) -> pd.DataFrame:
    tongue = tissue == "Tongue"
    rows = []
    for group in sorted(set(cell_type[tongue])):
        mask = tongue & (cell_type == group)
        n_cells = int(mask.sum())
        if n_cells < 20:
            continue
        group_expr = expression[mask]
        pct_expr = (group_expr > 0).mean(axis=0) * 100.0
        mean_expr = group_expr.mean(axis=0)
        for i, gene in enumerate(GENES):
            rows.append(
                {
                    "cell_type": group,
                    "gene": gene,
                    "n_cells": n_cells,
                    "pct_expressing": float(pct_expr[i]),
                    "mean_expression": float(mean_expr[i]),
                }
            )
    return pd.DataFrame(rows)


def plot(summary: pd.DataFrame) -> None:
    adora2b = summary[summary["gene"] == "ADORA2B"].sort_values(
        ["pct_expressing", "mean_expression"], ascending=False
    )
    groups = adora2b.head(16)["cell_type"].iloc[::-1].tolist()
    plot_df = summary[summary["cell_type"].isin(groups)].copy()
    group_to_y = {group: i for i, group in enumerate(groups)}
    gene_to_x = {gene: i for i, gene in enumerate(GENES)}
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["cell_type"].map(group_to_y)

    fig, ax = plt.subplots(figsize=(8.8, 6.6))
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=10 + plot_df["pct_expressing"].to_numpy() * 7.0,
        c=plot_df["mean_expression"],
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )
    ax.set_xticks(range(len(GENES)), GENES)
    ax.set_yticks(range(len(groups)), groups)
    ax.set_xlabel("Adenosine receptor gene")
    ax.set_ylabel("Tongue cell type")
    ax.set_title("Tabula Sapiens Tongue ADORA expression by cell type")
    ax.grid(axis="both", color="#e6e6e6", linewidth=0.7)
    ax.set_axisbelow(True)
    cbar = fig.colorbar(scatter, ax=ax, pad=0.03)
    cbar.set_label("Mean expression")
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")


def main() -> None:
    expression = np.load(EXPR_CACHE, allow_pickle=False)["expression"]
    with h5py.File(H5AD, "r") as f:
        tissue = read_categorical(f, "tissue_in_publication")
        cell_type = read_categorical(f, "cell_type")
    summary = summarize(expression, tissue, cell_type)
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_values(["gene", "pct_expressing"], ascending=[True, False]).to_csv(TABLE, index=False)
    print(f"Wrote {TABLE}")
    plot(summary)


if __name__ == "__main__":
    main()
