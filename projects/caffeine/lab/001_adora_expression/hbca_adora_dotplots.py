#!/usr/bin/env python
"""Extract and dotplot ADORA receptor expression in HBCA non-neuronal cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
HBCA_CACHE_DIR = LAB_DIR / "cache" / "human_brain_cell_atlas"
DEFAULT_H5AD = HBCA_CACHE_DIR / "hbca_all_non_neuronal_b165f033.h5ad"
DEFAULT_EXPR_CACHE = HBCA_CACHE_DIR / "hbca_adora_expression.npz"
DEFAULT_FIGURES_DIR = LAB_DIR / "figures"
GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_encoded_column(group: h5py.Group) -> np.ndarray:
    if isinstance(group, h5py.Group):
        if "categories" in group and "codes" in group:
            categories = decode_values(group["categories"][:])
            codes = group["codes"][:]
            return np.array([categories[code] if code >= 0 else "" for code in codes], dtype=object)
        raise ValueError(f"Unsupported group column with keys: {list(group.keys())}")
    return np.array(decode_values(group[:]), dtype=object)


def read_var_names(f: h5py.File) -> np.ndarray:
    if "feature_name" in f["var"]:
        return read_encoded_column(f["var/feature_name"])
    if "_index" in f["var"]:
        return read_encoded_column(f["var/_index"])
    raise ValueError("Could not find var/feature_name or var/_index")


def gene_indices(f: h5py.File, genes: tuple[str, ...]) -> np.ndarray:
    names = read_var_names(f)
    indices = []
    for gene in genes:
        matches = np.flatnonzero(names == gene)
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one var/feature_name match for {gene}, found {len(matches)}")
        indices.append(int(matches[0]))
    return np.array(indices, dtype=np.int64)


def extract_gene_expression(
    h5ad_path: Path,
    matrix_group: str,
    genes: tuple[str, ...],
    row_chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(h5ad_path, "r") as f:
        matrix = f[matrix_group]
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{matrix_group} is not a CSR matrix")

        gene_idx = gene_indices(f, genes)
        n_obs, _ = matrix.attrs["shape"]
        indptr = matrix["indptr"][:]
        indices_ds = matrix["indices"]
        data_ds = matrix["data"]
        expression = np.zeros((int(n_obs), len(genes)), dtype=np.float32)

        order = np.argsort(gene_idx)
        sorted_gene_idx = gene_idx[order]
        n_nnz = int(indices_ds.shape[0])
        for row_start in range(0, int(n_obs), row_chunk_size):
            row_stop = min(row_start + row_chunk_size, int(n_obs))
            nnz_start = int(indptr[row_start])
            nnz_stop = int(indptr[row_stop])
            if nnz_start == nnz_stop:
                continue

            indices = indices_ds[nnz_start:nnz_stop]
            positions = np.searchsorted(sorted_gene_idx, indices)
            in_range = positions < len(sorted_gene_idx)
            mask = np.zeros(indices.shape, dtype=bool)
            mask[in_range] = sorted_gene_idx[positions[in_range]] == indices[in_range]
            if np.any(mask):
                offsets = np.flatnonzero(mask)
                local_indptr = indptr[row_start : row_stop + 1] - nnz_start
                rows = np.searchsorted(local_indptr, offsets, side="right") - 1 + row_start
                cols = order[positions[mask]]
                expression[rows, cols] = data_ds[nnz_start:nnz_stop][mask]

            if row_stop == int(n_obs) or row_stop % (row_chunk_size * 10) == 0:
                print(
                    f"processed {row_stop:,}/{int(n_obs):,} cells; "
                    f"{int(np.count_nonzero(expression)):,} ADORA values found "
                    f"after {nnz_stop:,}/{n_nnz:,} nnz"
                )

    return expression, gene_idx


def load_or_extract_expression(
    h5ad_path: Path,
    cache_path: Path,
    matrix_group: str,
    genes: tuple[str, ...],
    row_chunk_size: int,
    force: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if cache_path.exists() and not force:
        cached = np.load(cache_path, allow_pickle=False)
        cached_genes = tuple(str(x) for x in cached["genes"])
        cached_matrix = str(cached["matrix_group"])
        if cached_genes == genes and cached_matrix == matrix_group:
            print(f"Using expression cache {cache_path}")
            return cached["expression"], cached["gene_indices"]

    expression, indices = extract_gene_expression(h5ad_path, matrix_group, genes, row_chunk_size)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        expression=expression,
        gene_indices=indices,
        genes=np.array(genes),
        matrix_group=np.array(matrix_group),
    )
    print(f"Wrote expression cache {cache_path}")
    return expression, indices


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
    group_by: str,
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
        positive_mean = np.divide(
            group_expr.sum(axis=0),
            (group_expr > 0).sum(axis=0),
            out=np.zeros(len(GENES), dtype=np.float32),
            where=(group_expr > 0).sum(axis=0) > 0,
        )
        for i, gene in enumerate(GENES):
            rows.append(
                {
                    "group_by": group_by,
                    "group": group,
                    "gene": gene,
                    "n_cells": n_cells,
                    "n_expressing": int(np.count_nonzero(group_expr[:, i] > 0)),
                    "pct_expressing": float(pct_expr[i]),
                    "mean_expression": float(mean_expr[i]),
                    "mean_expression_positive_cells": float(positive_mean[i]),
                }
            )
    return pd.DataFrame(rows)


def select_groups(summary: pd.DataFrame, top_n: int) -> list[str]:
    scores = (
        summary.assign(score=summary["pct_expressing"] * np.log1p(summary["mean_expression"]))
        .groupby("group", as_index=True)
        .agg(
            max_pct_expressing=("pct_expressing", "max"),
            max_mean_expression=("mean_expression", "max"),
            max_score=("score", "max"),
        )
        .sort_values(["max_score", "max_pct_expressing", "max_mean_expression"], ascending=False)
    )
    return scores.head(top_n).index[::-1].tolist()


def display_label(value: str) -> str:
    return value.replace("_", " ")


def plot_dotplot(
    summary: pd.DataFrame,
    groups: list[str],
    output_path: Path,
    title: str,
    ylabel: str,
    dpi: int,
) -> None:
    plot_df = summary[summary["group"].isin(groups)].copy()
    group_to_y = {group: i for i, group in enumerate(groups)}
    gene_to_x = {gene: i for i, gene in enumerate(GENES)}
    plot_df["x"] = plot_df["gene"].map(gene_to_x)
    plot_df["y"] = plot_df["group"].map(group_to_y)

    fig_height = max(5.6, 0.3 * len(groups) + 2.0)
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


def expression_thresholds(expression: np.ndarray, quantile: float) -> np.ndarray:
    thresholds = np.zeros(len(GENES), dtype=np.float32)
    for i, gene in enumerate(GENES):
        positive = expression[:, i][expression[:, i] > 0]
        if len(positive) == 0:
            raise ValueError(f"{gene} has no positive expression values")
        thresholds[i] = np.quantile(positive, quantile)
    return thresholds


def write_gene_summary(
    output_json: Path,
    output_csv: Path,
    expression: np.ndarray,
    gene_idx: np.ndarray,
    quantile: float,
    matrix_group: str,
) -> None:
    thresholds = expression_thresholds(expression, quantile)
    high = expression >= thresholds
    rows = []
    for i, gene in enumerate(GENES):
        rows.append(
            {
                "gene": gene,
                "var_index": int(gene_idx[i]),
                "n_cells_total": int(expression.shape[0]),
                "n_cells_positive": int(np.count_nonzero(expression[:, i] > 0)),
                "pct_cells_positive": float(np.count_nonzero(expression[:, i] > 0) / expression.shape[0] * 100.0),
                "mean_expression_all_cells": float(expression[:, i].mean()),
                "mean_expression_positive_cells": float(expression[expression[:, i] > 0, i].mean()),
                "high_threshold": float(thresholds[i]),
                "n_cells_high": int(np.count_nonzero(high[:, i])),
                "pct_cells_high": float(np.count_nonzero(high[:, i]) / expression.shape[0] * 100.0),
            }
        )
    summary = {
        "dataset": "Human Brain Cell Atlas v1.0, all non-neuronal cells",
        "matrix_group": matrix_group,
        "threshold_rule": f"per-gene q{quantile:g} among cells with expression > 0",
        "n_cells_total": int(expression.shape[0]),
        "n_cells_any_positive": int(np.count_nonzero((expression > 0).any(axis=1))),
        "n_cells_any_high": int(np.count_nonzero(high.any(axis=1))),
        "genes": {row["gene"]: row for row in rows},
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame(rows).to_csv(output_csv, index=False)
    print(f"Wrote {output_json}")
    print(f"Wrote {output_csv}")


def run_dotplot(
    h5ad_path: Path,
    expression: np.ndarray,
    group_by: str,
    min_cells: int,
    top_n: int,
    figures_dir: Path,
    dpi: int,
) -> None:
    codes, categories = read_categorical(h5ad_path, group_by)
    if expression.shape[0] != codes.shape[0]:
        raise ValueError(f"Expression rows ({expression.shape[0]}) do not match obs rows ({codes.shape[0]})")

    summary = summarize_by_group(expression, codes, categories, group_by, min_cells)
    groups = select_groups(summary, top_n)
    table_path = figures_dir / f"hbca_adora_dotplot_{group_by}.csv"
    figure_path = figures_dir / f"hbca_adora_dotplot_{group_by}.png"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_values(["group", "gene"]).to_csv(table_path, index=False)
    print(f"Wrote {table_path}")

    title = f"HBCA non-neuronal ADORA expression by {display_label(group_by)}"
    plot_dotplot(summary, groups, figure_path, title, display_label(group_by).title(), dpi)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--matrix-group", default="X")
    parser.add_argument("--expression-cache", type=Path, default=DEFAULT_EXPR_CACHE)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--group-by", nargs="+", default=["cell_type", "supercluster_term"])
    parser.add_argument("--row-chunk-size", type=int, default=20_000)
    parser.add_argument("--min-cells", type=int, default=50)
    parser.add_argument("--top-n", type=int, default=35)
    parser.add_argument("--threshold-quantile", type=float, default=0.75)
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expression, gene_idx = load_or_extract_expression(
        args.h5ad,
        args.expression_cache,
        args.matrix_group,
        GENES,
        args.row_chunk_size,
        args.force,
    )
    write_gene_summary(
        args.figures_dir / "hbca_adora_expression_summary.json",
        args.figures_dir / "hbca_adora_expression_summary.csv",
        expression,
        gene_idx,
        args.threshold_quantile,
        args.matrix_group,
    )
    for group_by in args.group_by:
        run_dotplot(args.h5ad, expression, group_by, args.min_cells, args.top_n, args.figures_dir, args.dpi)


if __name__ == "__main__":
    main()
