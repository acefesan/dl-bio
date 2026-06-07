#!/usr/bin/env python
"""Dotplot of cytoskeleton gene expression by Tabula Sapiens cell type.

Genes are shown on the x-axis, grouped by structural family with vertical
separators. Cell types are on the y-axis, ranked so that types with strong
specialized cytoskeleton (muscle actins, neurofilaments, vimentin, keratins)
float to the top.

Dot size  = % of cells in that type expressing the gene (> 0).
Dot color = mean expression across all cells of that type.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_cyto_pergene_expression.npz"
DEFAULT_OUT = LAB_DIR / "figures" / "tabula_sapiens_cyto_dotplot_cell_type.png"
DEFAULT_TABLE = LAB_DIR / "figures" / "tabula_sapiens_cyto_dotplot_cell_type.csv"

# Ordered gene list, grouped by family.  Order within each family is preserved
# on the x-axis so the plot reads left-to-right from ubiquitous → specialized.
GENE_FAMILIES: dict[str, list[str]] = {
    "Actin\n(ubiquitous)": ["ACTB", "ACTG1"],
    "Actin\n(muscle)": ["ACTA1", "ACTA2", "ACTC1", "ACTG2"],
    "Tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBB", "TUBB2A", "TUBB4B", "TUBG1"],
    "Neuro-\nfilaments": ["NEFL", "NEFM", "NEFH"],
    "Vimentin": ["VIM"],
    "Keratins": ["KRT1", "KRT5", "KRT8", "KRT14", "KRT18", "KRT19"],
}

FAMILY_COLORS: dict[str, str] = {
    "Actin\n(ubiquitous)": "#2196F3",
    "Actin\n(muscle)": "#F44336",
    "Tubulin": "#FF9800",
    "Neuro-\nfilaments": "#9C27B0",
    "Vimentin": "#009688",
    "Keratins": "#E91E63",
}

ALL_GENES_ORDERED = [g for genes in GENE_FAMILIES.values() for g in genes]

# Genes considered "specialized" for ranking cell types
SPECIALIZED = {"ACTA1", "ACTA2", "ACTC1", "ACTG2", "NEFL", "NEFM", "NEFH", "VIM",
               "KRT1", "KRT5", "KRT8", "KRT14", "KRT18", "KRT19"}


# ---------------------------------------------------------------------------
# H5AD helpers
# ---------------------------------------------------------------------------

def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_var_names(f: h5py.File) -> np.ndarray:
    node = f["var/feature_name"]
    if isinstance(node, h5py.Group):
        cats = decode_values(node["categories"][:])
        codes = node["codes"][:]
        return np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
    return np.array(decode_values(node[:]), dtype=object)


def read_categorical(h5ad_path: Path, obs_column: str) -> tuple[np.ndarray, list[str]]:
    with h5py.File(h5ad_path, "r") as f:
        node = f["obs"][obs_column]
        if not isinstance(node, h5py.Group) or "codes" not in node:
            raise ValueError(f"obs/{obs_column!r} is not categorical")
        return node["codes"][:], decode_values(node["categories"][:])


# ---------------------------------------------------------------------------
# Per-gene expression extraction (tolerant of missing genes)
# ---------------------------------------------------------------------------

def find_gene_indices(f: h5py.File, genes: list[str]) -> dict[str, int]:
    names = read_var_names(f)
    found: dict[str, int] = {}
    for gene in genes:
        matches = np.flatnonzero(names == gene)
        if len(matches) == 1:
            found[gene] = int(matches[0])
        elif len(matches) > 1:
            print(f"  Warning: {gene} has {len(matches)} matches — skipping")
        else:
            print(f"  Warning: {gene} not found — skipping")
    return found


def extract_per_gene_expression(
    h5ad_path: Path,
    matrix_group: str,
    genes: list[str],
    chunk_nnz: int,
) -> tuple[np.ndarray, list[str]]:
    """Return (n_obs, n_found) dense float32 array and the list of found gene names."""
    with h5py.File(h5ad_path, "r") as f:
        gene_idx_map = find_gene_indices(f, genes)
        found_genes = [g for g in genes if g in gene_idx_map]  # preserve order
        target_var = np.array([gene_idx_map[g] for g in found_genes], dtype=np.int64)

        matrix = f[matrix_group]
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{matrix_group} is not a CSR matrix")

        n_obs, _ = matrix.attrs["shape"]
        indptr = matrix["indptr"][:]
        indices_ds = matrix["indices"]
        data_ds = matrix["data"]
        expr = np.zeros((int(n_obs), len(found_genes)), dtype=np.float32)

        order = np.argsort(target_var)
        sorted_var = target_var[order]

        n_nnz = int(indices_ds.shape[0])
        for start in range(0, n_nnz, chunk_nnz):
            stop = min(start + chunk_nnz, n_nnz)
            col_indices = indices_ds[start:stop]
            positions = np.searchsorted(sorted_var, col_indices)
            in_range = positions < len(sorted_var)
            mask = np.zeros(col_indices.shape, dtype=bool)
            mask[in_range] = sorted_var[positions[in_range]] == col_indices[in_range]
            if not np.any(mask):
                continue

            offsets = np.flatnonzero(mask)
            rows = np.searchsorted(indptr, start + offsets, side="right") - 1
            cols = order[positions[mask]]
            expr[rows, cols] = data_ds[start:stop][mask].astype(np.float32)
            print(f"  extracted {int(np.count_nonzero(expr)):,} values after {stop:,}/{n_nnz:,} nnz")

    print(f"  Found {len(found_genes)}/{len(genes)} requested genes")
    return expr, found_genes


def load_or_extract(
    h5ad_path: Path,
    cache_path: Path,
    matrix_group: str,
    genes: list[str],
    chunk_nnz: int,
    force: bool,
) -> tuple[np.ndarray, list[str]]:
    if cache_path.exists() and not force:
        cached = np.load(cache_path, allow_pickle=True)
        found_genes = [str(g) for g in cached["genes_found"].tolist()]
        if (
            set(found_genes).issubset(set(genes))
            and str(cached["matrix_group"]) == matrix_group
        ):
            return cached["expression"], found_genes

    expr, found_genes = extract_per_gene_expression(h5ad_path, matrix_group, genes, chunk_nnz)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        expression=expr,
        genes_found=np.array(found_genes),
        matrix_group=np.array(matrix_group),
    )
    print(f"Wrote {cache_path}")
    return expr, found_genes


# ---------------------------------------------------------------------------
# Summarize expression per cell type
# ---------------------------------------------------------------------------

def summarize_by_cell_type(
    expression: np.ndarray,
    found_genes: list[str],
    codes: np.ndarray,
    categories: list[str],
    min_cells: int,
) -> pd.DataFrame:
    rows = []
    for code, ct in enumerate(categories):
        mask = codes == code
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        sub = expression[mask]
        pct = (sub > 0).mean(axis=0) * 100.0
        mean = sub.mean(axis=0)
        for j, gene in enumerate(found_genes):
            rows.append({"cell_type": ct, "gene": gene, "n_cells": n_cells,
                         "pct_expressing": float(pct[j]), "mean_expression": float(mean[j])})
    return pd.DataFrame(rows)


def rank_cell_types(summary: pd.DataFrame, top_n: int) -> list[str]:
    """Rank by peak specialized gene expression, then fall back to ubiquitous score."""
    spec_score = (
        summary[summary["gene"].isin(SPECIALIZED)]
        .assign(score=lambda d: d["pct_expressing"] * np.log1p(d["mean_expression"]))
        .groupby("cell_type")["score"]
        .max()
        .rename("spec_score")
    )
    all_score = (
        summary
        .assign(score=lambda d: d["pct_expressing"] * np.log1p(d["mean_expression"]))
        .groupby("cell_type")["score"]
        .max()
        .rename("all_score")
    )
    combined = pd.concat([spec_score, all_score], axis=1).fillna(0)
    combined = combined.sort_values(["spec_score", "all_score"], ascending=False)
    return combined.head(top_n).index[::-1].tolist()


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_dotplot(
    summary: pd.DataFrame,
    cell_types: list[str],
    found_genes: list[str],
    output_path: Path,
    title: str,
    dpi: int,
) -> None:
    # Only keep genes actually found, in the canonical family order
    genes_plot = [g for g in ALL_GENES_ORDERED if g in found_genes]

    ct_to_y = {ct: i for i, ct in enumerate(cell_types)}
    gene_to_x = {g: i for i, g in enumerate(genes_plot)}

    df = summary[summary["cell_type"].isin(cell_types) & summary["gene"].isin(genes_plot)].copy()
    df["x"] = df["gene"].map(gene_to_x)
    df["y"] = df["cell_type"].map(ct_to_y)

    fig_h = max(8.0, 0.26 * len(cell_types) + 3.0)
    fig_w = max(10.0, 0.52 * len(genes_plot) + 3.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    max_expr = df["mean_expression"].max()
    sizes = 4 + df["pct_expressing"].to_numpy() * 4.2
    scatter = ax.scatter(
        df["x"], df["y"],
        s=sizes,
        c=df["mean_expression"],
        cmap="viridis",
        vmin=0, vmax=max_expr,
        edgecolors="#444444",
        linewidths=0.2,
        zorder=3,
    )

    # Family separator lines and top-of-axes family labels
    x_cursor = 0
    family_label_y = len(cell_types) - 0.5
    for family, fgenes in GENE_FAMILIES.items():
        present = [g for g in fgenes if g in gene_to_x]
        if not present:
            continue
        xs = [gene_to_x[g] for g in present]
        color = FAMILY_COLORS[family]
        ax.axvline(min(xs) - 0.5, color="#bbbbbb", linewidth=0.7, zorder=1)
        mid = (min(xs) + max(xs)) / 2
        ax.text(mid, family_label_y + 0.8, family.replace("\n", " "),
                ha="center", va="bottom", fontsize=7.5, color=color, fontweight="bold")
        ax.hlines(family_label_y + 0.1, min(xs) - 0.4, max(xs) + 0.4,
                  color=color, linewidth=1.5)

    ax.set_xticks(range(len(genes_plot)), genes_plot, fontsize=7.5)
    ax.tick_params(axis="x", labelrotation=40)
    for label, gene in zip(ax.get_xticklabels(), genes_plot):
        label.set_horizontalalignment("right")
        for family, fgenes in GENE_FAMILIES.items():
            if gene in fgenes:
                label.set_color(FAMILY_COLORS[family])
                break

    ax.set_yticks(range(len(cell_types)), cell_types, fontsize=7.5)
    ax.set_xlim(-0.7, len(genes_plot) - 0.3)
    ax.set_ylim(-1.0, len(cell_types) + 1.8)
    ax.grid(axis="both", color="#ececec", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=11, pad=28)
    ax.set_xlabel("Cytoskeleton gene", labelpad=8)
    ax.set_ylabel("Cell type")

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, shrink=0.55)
    cbar.set_label("Mean expression", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    legend_pcts = [5, 25, 50, 80]
    handles = [
        ax.scatter([], [], s=4 + p * 4.2, facecolors="none",
                   edgecolors="#444444", linewidths=0.5)
        for p in legend_pcts
    ]
    ax.legend(handles, [f"{p}%" for p in legend_pcts],
              title="% expressing", title_fontsize=7.5,
              frameon=False, fontsize=7.5,
              loc="lower right", bbox_to_anchor=(1.18, 0.0))

    fig.tight_layout(rect=(0, 0, 0.88, 1))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    p.add_argument("--matrix-group", default="X")
    p.add_argument("--expression-cache", type=Path, default=DEFAULT_EXPR_CACHE)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    p.add_argument("--group-by", default="cell_type")
    p.add_argument("--min-cells", type=int, default=50)
    p.add_argument("--top-n", type=int, default=40)
    p.add_argument("--chunk-nnz", type=int, default=25_000_000)
    p.add_argument("--force-extract", action="store_true")
    p.add_argument("--dpi", type=int, default=240)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    print("Loading / extracting expression …")
    expr, found_genes = load_or_extract(
        h5ad_path=args.h5ad,
        cache_path=args.expression_cache,
        matrix_group=args.matrix_group,
        genes=ALL_GENES_ORDERED,
        chunk_nnz=args.chunk_nnz,
        force=args.force_extract,
    )

    print("Reading cell type labels …")
    codes, categories = read_categorical(args.h5ad, args.group_by)
    if expr.shape[0] != codes.shape[0]:
        raise ValueError(f"Expression rows {expr.shape[0]} != obs rows {codes.shape[0]}")

    print("Summarizing …")
    summary = summarize_by_cell_type(expr, found_genes, codes, categories, args.min_cells)
    cell_types = rank_cell_types(summary, args.top_n)

    args.table.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_values(["cell_type", "gene"]).to_csv(args.table, index=False)
    print(f"Wrote {args.table}")

    title = (
        f"Cytoskeleton gene expression by {args.group_by} "
        f"(top {len(cell_types)} types ranked by specialized cytoskeleton)"
    )
    plot_dotplot(summary, cell_types, found_genes, args.out, title, args.dpi)


if __name__ == "__main__":
    main()
