#!/usr/bin/env python
"""Plot Tabula Sapiens UMAP cells with high ADORA receptor expression."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
DEFAULT_FIGURE = LAB_DIR / "figures" / "tabula_sapiens_adora_high_umap.png"
DEFAULT_SUMMARY = LAB_DIR / "figures" / "tabula_sapiens_adora_high_umap_summary.json"
GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")
COLORS = {
    "ADORA1": "#1b9e77",
    "ADORA2A": "#d95f02",
    "ADORA2B": "#7570b3",
    "ADORA3": "#e7298a",
}


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_var_names(f: h5py.File) -> np.ndarray:
    node = f["var/feature_name"]
    if isinstance(node, h5py.Group):
        categories = decode_values(node["categories"][:])
        codes = node["codes"][:]
        return np.array([categories[code] if code >= 0 else "" for code in codes], dtype=object)
    return np.array(decode_values(node[:]), dtype=object)


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
    chunk_nnz: int,
) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(h5ad_path, "r") as f:
        gene_idx = gene_indices(f, genes)
        matrix = f[matrix_group]
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{matrix_group} is not a CSR matrix")

        n_obs, _ = matrix.attrs["shape"]
        indptr = matrix["indptr"][:]
        indices_ds = matrix["indices"]
        data_ds = matrix["data"]
        expr = np.zeros((int(n_obs), len(genes)), dtype=np.float32)

        order = np.argsort(gene_idx)
        sorted_gene_idx = gene_idx[order]
        n_nnz = int(indices_ds.shape[0])
        for start in range(0, n_nnz, chunk_nnz):
            stop = min(start + chunk_nnz, n_nnz)
            indices = indices_ds[start:stop]
            positions = np.searchsorted(sorted_gene_idx, indices)
            in_range = positions < len(sorted_gene_idx)
            mask = np.zeros(indices.shape, dtype=bool)
            mask[in_range] = sorted_gene_idx[positions[in_range]] == indices[in_range]
            if not np.any(mask):
                continue

            offsets = np.flatnonzero(mask)
            rows = np.searchsorted(indptr, start + offsets, side="right") - 1
            cols = order[positions[mask]]
            values = data_ds[start:stop][mask]
            expr[rows, cols] = values
            print(f"extracted {int(np.count_nonzero(expr)):,} receptor values after {stop:,}/{n_nnz:,} nnz")

    return expr, gene_idx


def load_or_extract_expression(
    h5ad_path: Path,
    cache_path: Path,
    matrix_group: str,
    genes: tuple[str, ...],
    chunk_nnz: int,
    force: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if cache_path.exists() and not force:
        cached = np.load(cache_path, allow_pickle=False)
        cached_genes = tuple(str(x) for x in cached["genes"])
        cached_matrix = str(cached["matrix_group"])
        if cached_genes == genes and cached_matrix == matrix_group:
            return cached["expression"], cached["gene_indices"]

    expression, indices = extract_gene_expression(h5ad_path, matrix_group, genes, chunk_nnz)
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


def expression_thresholds(expression: np.ndarray, genes: tuple[str, ...], quantile: float) -> np.ndarray:
    thresholds = np.zeros(len(genes), dtype=np.float32)
    for i, gene in enumerate(genes):
        positive = expression[:, i][expression[:, i] > 0]
        if len(positive) == 0:
            raise ValueError(f"{gene} has no positive expression values")
        thresholds[i] = np.quantile(positive, quantile)
    return thresholds


def assign_high_cells(expression: np.ndarray, thresholds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    high = expression >= thresholds
    score = np.divide(
        expression,
        thresholds,
        out=np.zeros_like(expression, dtype=np.float32),
        where=thresholds[None, :] > 0,
    )
    assigned = np.full(expression.shape[0], -1, dtype=np.int8)
    any_high = high.any(axis=1)
    assigned[any_high] = np.argmax(score[any_high], axis=1).astype(np.int8)
    return high, assigned


def plot_umap(
    h5ad_path: Path,
    output_path: Path,
    expression: np.ndarray,
    thresholds: np.ndarray,
    assigned: np.ndarray,
    genes: tuple[str, ...],
    point_size_background: float,
    point_size_high: float,
    alpha_background: float,
    alpha_high: float,
    hide_background: bool,
    dpi: int,
) -> None:
    with h5py.File(h5ad_path, "r") as f:
        umap = f["obsm/X_umap"][:]

    fig, ax = plt.subplots(figsize=(8.0, 7.2))
    if not hide_background:
        ax.scatter(
            umap[:, 0],
            umap[:, 1],
            c="#d4d4d4",
            s=point_size_background,
            alpha=alpha_background,
            linewidths=0,
            rasterized=True,
        )

    for i, gene in enumerate(genes):
        mask = assigned == i
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=COLORS[gene],
            s=point_size_high,
            alpha=alpha_high,
            linewidths=0,
            rasterized=True,
            label=f"{gene} high",
        )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=7,
            markerfacecolor=COLORS[gene],
            markeredgecolor="none",
            label=f"{gene} high (>= {thresholds[i]:.3g})",
        )
        for i, gene in enumerate(genes)
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=9)
    ax.set_title("Tabula Sapiens cells with high adenosine receptor expression", fontsize=13)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#b8b8b8")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


def write_summary(
    summary_path: Path,
    expression: np.ndarray,
    thresholds: np.ndarray,
    high: np.ndarray,
    assigned: np.ndarray,
    genes: tuple[str, ...],
    gene_idx: np.ndarray,
    quantile: float,
    matrix_group: str,
) -> None:
    summary = {
        "matrix_group": matrix_group,
        "threshold_rule": f"per-gene q{quantile:g} among cells with expression > 0",
        "n_cells_total": int(expression.shape[0]),
        "n_cells_any_high": int(np.count_nonzero(assigned >= 0)),
        "n_cells_high_by_gene_independent": {
            gene: int(np.count_nonzero(high[:, i])) for i, gene in enumerate(genes)
        },
        "n_cells_assigned_to_gene": {
            gene: int(np.count_nonzero(assigned == i)) for i, gene in enumerate(genes)
        },
        "genes": {
            gene: {
                "var_index": int(gene_idx[i]),
                "n_cells_positive": int(np.count_nonzero(expression[:, i] > 0)),
                "threshold": float(thresholds[i]),
            }
            for i, gene in enumerate(genes)
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--matrix-group", default="X")
    parser.add_argument("--expression-cache", type=Path, default=DEFAULT_EXPR_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--threshold-quantile", type=float, default=0.75)
    parser.add_argument("--chunk-nnz", type=int, default=25_000_000)
    parser.add_argument("--force-extract", action="store_true")
    parser.add_argument("--point-size-background", type=float, default=0.035)
    parser.add_argument("--point-size-high", type=float, default=0.22)
    parser.add_argument("--alpha-background", type=float, default=0.18)
    parser.add_argument("--alpha-high", type=float, default=0.86)
    parser.add_argument("--hide-background", action="store_true", help="Plot only cells assigned to a receptor.")
    parser.add_argument("--dpi", type=int, default=260)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expression, gene_idx = load_or_extract_expression(
        h5ad_path=args.h5ad,
        cache_path=args.expression_cache,
        matrix_group=args.matrix_group,
        genes=GENES,
        chunk_nnz=args.chunk_nnz,
        force=args.force_extract,
    )
    thresholds = expression_thresholds(expression, GENES, args.threshold_quantile)
    high, assigned = assign_high_cells(expression, thresholds)
    plot_umap(
        h5ad_path=args.h5ad,
        output_path=args.out,
        expression=expression,
        thresholds=thresholds,
        assigned=assigned,
        genes=GENES,
        point_size_background=args.point_size_background,
        point_size_high=args.point_size_high,
        alpha_background=args.alpha_background,
        alpha_high=args.alpha_high,
        hide_background=args.hide_background,
        dpi=args.dpi,
    )
    write_summary(
        summary_path=args.summary,
        expression=expression,
        thresholds=thresholds,
        high=high,
        assigned=assigned,
        genes=GENES,
        gene_idx=gene_idx,
        quantile=args.threshold_quantile,
        matrix_group=args.matrix_group,
    )


if __name__ == "__main__":
    main()
