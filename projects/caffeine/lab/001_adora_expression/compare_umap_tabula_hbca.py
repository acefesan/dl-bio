#!/usr/bin/env python
"""Compare independent Tabula Sapiens and HBCA UMAP embeddings side by side."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_TABULA_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_HBCA_H5AD = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_all_non_neuronal_b165f033.h5ad"
DEFAULT_TABULA_EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
DEFAULT_HBCA_EXPR_CACHE = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_adora_expression.npz"
DEFAULT_CELL_TYPE_OUT = LAB_DIR / "figures" / "compare_umap_tabula_hbca_cell_type.png"
DEFAULT_ADORA_OUT = LAB_DIR / "figures" / "compare_umap_tabula_hbca_adora_high.png"
DEFAULT_SUMMARY = LAB_DIR / "figures" / "compare_umap_tabula_hbca_summary.json"
GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")
ADORA_COLORS = {
    "ADORA1": "#1b9e77",
    "ADORA2A": "#d95f02",
    "ADORA2B": "#7570b3",
    "ADORA3": "#e7298a",
}


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_categorical(f: h5py.File, obs_column: str) -> tuple[np.ndarray, list[str]]:
    node = f["obs"][obs_column]
    if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
        raise ValueError(f"obs/{obs_column!r} is not an AnnData categorical column")
    return node["codes"][:], decode_values(node["categories"][:])


def maybe_subsample(n_obs: int, max_points: int | None, seed: int) -> np.ndarray | slice:
    if max_points is None or max_points >= n_obs:
        return slice(None)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_points, replace=False))


def subset_array(values: np.ndarray, idx: np.ndarray | slice) -> np.ndarray:
    if isinstance(idx, slice):
        return values
    return values[idx]


def top_category_codes(codes: np.ndarray, top_n: int) -> np.ndarray:
    valid = codes[codes >= 0]
    if valid.size == 0:
        return np.array([], dtype=np.int64)
    counts = np.bincount(valid)
    ranked = np.argsort(counts)[::-1]
    return ranked[:top_n]


def dense_codes_for_top(codes: np.ndarray, top_codes: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    code_to_dense = {int(code): i for i, code in enumerate(top_codes)}
    dense = np.full(codes.shape[0], -1, dtype=np.int16)
    for code, dense_code in code_to_dense.items():
        dense[codes == code] = dense_code
    return dense, code_to_dense


def make_palette(n: int) -> list[tuple[float, float, float, float]]:
    colors = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        colors.extend(cmap(i) for i in range(cmap.N))
    if n > len(colors):
        colors.extend(plt.get_cmap("hsv")(np.linspace(0, 1, n - len(colors), endpoint=False)))
    return colors[:n]


def style_embedding_axis(ax: plt.Axes, xlabel: str, ylabel: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#b8b8b8")


def draw_local_cell_type_panel(
    ax: plt.Axes,
    coords: np.ndarray,
    codes: np.ndarray,
    categories: list[str],
    title: str,
    top_n: int,
    seed: int,
    point_size: float,
    alpha: float,
) -> list[Line2D]:
    top_codes = top_category_codes(codes, top_n)
    dense, _ = dense_codes_for_top(codes, top_codes)
    palette = make_palette(len(top_codes))
    order = np.random.default_rng(seed).permutation(coords.shape[0])

    other = dense[order] < 0
    if np.any(other):
        other_order = order[other]
        ax.scatter(
            coords[other_order, 0],
            coords[other_order, 1],
            c="#d6d6d6",
            s=point_size,
            alpha=0.18,
            linewidths=0,
            rasterized=True,
        )

    for dense_code, color in enumerate(palette):
        mask = dense[order] == dense_code
        if not np.any(mask):
            continue
        draw_idx = order[mask]
        ax.scatter(
            coords[draw_idx, 0],
            coords[draw_idx, 1],
            c=[color],
            s=point_size,
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

    ax.set_title(title, fontsize=12)
    style_embedding_axis(ax, "UMAP 1", "UMAP 2")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=4.5,
            markerfacecolor=palette[i],
            markeredgecolor="none",
            label=categories[int(code)],
        )
        for i, code in enumerate(top_codes)
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=4.5,
            markerfacecolor="#d6d6d6",
            markeredgecolor="none",
            label="other cell types",
        )
    )
    return handles


def plot_cell_type_comparison(
    tabula_h5ad: Path,
    hbca_h5ad: Path,
    output_path: Path,
    max_points: int,
    top_n_cell_types: int,
    seed: int,
    point_size: float,
    alpha: float,
    dpi: int,
) -> dict[str, int]:
    with h5py.File(tabula_h5ad, "r") as tabula, h5py.File(hbca_h5ad, "r") as hbca:
        tabula_codes, tabula_categories = read_categorical(tabula, "cell_type")
        hbca_codes, hbca_categories = read_categorical(hbca, "cell_type")
        tabula_idx = maybe_subsample(tabula_codes.shape[0], max_points, seed)
        hbca_idx = maybe_subsample(hbca_codes.shape[0], max_points, seed + 1)
        tabula_coords = tabula["obsm/X_umap"][tabula_idx, :2]
        hbca_coords = hbca["obsm/X_UMAP"][hbca_idx, :2]
        tabula_plot_codes = subset_array(tabula_codes, tabula_idx)
        hbca_plot_codes = subset_array(hbca_codes, hbca_idx)

    fig, axes = plt.subplots(1, 2, figsize=(16.8, 7.2))
    tabula_handles = draw_local_cell_type_panel(
        axes[0],
        tabula_coords,
        tabula_plot_codes,
        tabula_categories,
        f"Tabula Sapiens X_umap\ncell_type labels ({len(tabula_plot_codes):,}/{len(tabula_codes):,} cells)",
        top_n_cell_types,
        seed,
        point_size,
        alpha,
    )
    hbca_handles = draw_local_cell_type_panel(
        axes[1],
        hbca_coords,
        hbca_plot_codes,
        hbca_categories,
        f"Human Brain Cell Atlas X_UMAP\ncell_type labels ({len(hbca_plot_codes):,}/{len(hbca_codes):,} cells)",
        top_n_cell_types,
        seed + 1,
        point_size,
        alpha,
    )

    axes[0].legend(
        handles=tabula_handles,
        title=f"Tabula top {top_n_cell_types}",
        loc="center left",
        bbox_to_anchor=(-0.02, -0.17),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        ncol=2,
    )
    axes[1].legend(
        handles=hbca_handles,
        title=f"HBCA top {top_n_cell_types}",
        loc="center left",
        bbox_to_anchor=(-0.02, -0.17),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        ncol=2,
    )
    fig.suptitle(
        "Independent UMAPs: coordinates are dataset-local and are not a merged atlas",
        fontsize=14,
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        "Each panel uses its own UMAP fit and its own local cell_type labels; colors are not shared across panels.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")
    return {
        "tabula_cells_total": int(tabula_codes.shape[0]),
        "tabula_cells_plotted": int(len(tabula_plot_codes)),
        "hbca_cells_total": int(hbca_codes.shape[0]),
        "hbca_cells_plotted": int(len(hbca_plot_codes)),
        "top_n_cell_types_per_panel": int(top_n_cell_types),
    }


def load_expression_cache(path: Path, expected_n_cells: int) -> tuple[np.ndarray, tuple[str, ...]]:
    cached = np.load(path, allow_pickle=False)
    expression = cached["expression"]
    genes = tuple(str(x) for x in cached["genes"]) if "genes" in cached.files else GENES
    if expression.shape[0] != expected_n_cells:
        raise ValueError(f"{path} has {expression.shape[0]:,} rows, expected {expected_n_cells:,}")
    missing = [gene for gene in GENES if gene not in genes]
    if missing:
        raise ValueError(f"{path} is missing genes: {', '.join(missing)}")
    order = [genes.index(gene) for gene in GENES]
    return expression[:, order], GENES


def expression_thresholds(expression: np.ndarray, quantile: float) -> np.ndarray:
    thresholds = np.zeros(expression.shape[1], dtype=np.float32)
    for i, gene in enumerate(GENES):
        positive = expression[:, i][expression[:, i] > 0]
        if positive.size == 0:
            raise ValueError(f"{gene} has no positive expression values")
        thresholds[i] = np.quantile(positive, quantile)
    return thresholds


def assign_high_cells(expression: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
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
    return assigned


def draw_adora_high_panel(
    ax: plt.Axes,
    coords: np.ndarray,
    assigned: np.ndarray,
    title: str,
    seed: int,
    point_size_background: float,
    point_size_high: float,
) -> None:
    order = np.random.default_rng(seed).permutation(coords.shape[0])
    ax.scatter(
        coords[order, 0],
        coords[order, 1],
        c="#dddddd",
        s=point_size_background,
        alpha=0.16,
        linewidths=0,
        rasterized=True,
    )
    for i, gene in enumerate(GENES):
        mask = assigned[order] == i
        if not np.any(mask):
            continue
        draw_idx = order[mask]
        ax.scatter(
            coords[draw_idx, 0],
            coords[draw_idx, 1],
            c=ADORA_COLORS[gene],
            s=point_size_high,
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
    ax.set_title(title, fontsize=12)
    style_embedding_axis(ax, "UMAP 1", "UMAP 2")


def plot_adora_comparison(
    tabula_h5ad: Path,
    hbca_h5ad: Path,
    tabula_expression_cache: Path,
    hbca_expression_cache: Path,
    output_path: Path,
    quantile: float,
    dpi: int,
) -> dict[str, object] | None:
    if not tabula_expression_cache.exists() or not hbca_expression_cache.exists():
        print(
            "Skipping ADORA-high comparison; missing cache(s): "
            f"{tabula_expression_cache if not tabula_expression_cache.exists() else ''} "
            f"{hbca_expression_cache if not hbca_expression_cache.exists() else ''}"
        )
        return None

    with h5py.File(tabula_h5ad, "r") as tabula, h5py.File(hbca_h5ad, "r") as hbca:
        tabula_coords = tabula["obsm/X_umap"][:, :2]
        hbca_coords = hbca["obsm/X_UMAP"][:, :2]

    tabula_expr, _ = load_expression_cache(tabula_expression_cache, tabula_coords.shape[0])
    hbca_expr, _ = load_expression_cache(hbca_expression_cache, hbca_coords.shape[0])
    tabula_thresholds = expression_thresholds(tabula_expr, quantile)
    hbca_thresholds = expression_thresholds(hbca_expr, quantile)
    tabula_assigned = assign_high_cells(tabula_expr, tabula_thresholds)
    hbca_assigned = assign_high_cells(hbca_expr, hbca_thresholds)

    fig, axes = plt.subplots(1, 2, figsize=(16.4, 7.0))
    draw_adora_high_panel(
        axes[0],
        tabula_coords,
        tabula_assigned,
        f"Tabula Sapiens X_umap\nADORA-high cells ({np.count_nonzero(tabula_assigned >= 0):,}/{tabula_assigned.size:,})",
        seed=17,
        point_size_background=0.035,
        point_size_high=0.12,
    )
    draw_adora_high_panel(
        axes[1],
        hbca_coords,
        hbca_assigned,
        f"Human Brain Cell Atlas X_UMAP\nADORA-high cells ({np.count_nonzero(hbca_assigned >= 0):,}/{hbca_assigned.size:,})",
        seed=19,
        point_size_background=0.04,
        point_size_high=0.14,
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=ADORA_COLORS[gene],
            markeredgecolor="none",
            label=f"{gene} high",
        )
        for gene in GENES
    ]
    fig.legend(handles=handles, loc="lower center", frameon=False, ncol=4, fontsize=9)
    fig.suptitle(
        f"Independent UMAPs with ADORA-high overlay (q{quantile:g} among nonzero cells per gene)",
        fontsize=14,
        y=0.985,
    )
    fig.text(
        0.5,
        0.035,
        "UMAP coordinates are not shared between datasets; compare enriched neighborhoods and labels, not x/y positions.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")

    return {
        "threshold_rule": f"per-gene q{quantile:g} among cells with expression > 0",
        "tabula_thresholds": {gene: float(tabula_thresholds[i]) for i, gene in enumerate(GENES)},
        "hbca_thresholds": {gene: float(hbca_thresholds[i]) for i, gene in enumerate(GENES)},
        "tabula_assigned_high": {gene: int(np.count_nonzero(tabula_assigned == i)) for i, gene in enumerate(GENES)},
        "hbca_assigned_high": {gene: int(np.count_nonzero(hbca_assigned == i)) for i, gene in enumerate(GENES)},
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tabula-h5ad", type=Path, default=DEFAULT_TABULA_H5AD)
    parser.add_argument("--hbca-h5ad", type=Path, default=DEFAULT_HBCA_H5AD)
    parser.add_argument("--tabula-expression-cache", type=Path, default=DEFAULT_TABULA_EXPR_CACHE)
    parser.add_argument("--hbca-expression-cache", type=Path, default=DEFAULT_HBCA_EXPR_CACHE)
    parser.add_argument("--cell-type-out", type=Path, default=DEFAULT_CELL_TYPE_OUT)
    parser.add_argument("--adora-out", type=Path, default=DEFAULT_ADORA_OUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--max-points", type=int, default=180_000)
    parser.add_argument("--top-n-cell-types", type=int, default=16)
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--point-size", type=float, default=0.055)
    parser.add_argument("--alpha", type=float, default=0.55)
    parser.add_argument("--adora-quantile", type=float, default=0.75)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary: dict[str, object] = {
        "note": "Tabula Sapiens X_umap and HBCA X_UMAP are independent UMAP coordinate systems, not a merged embedding.",
        "cell_type_comparison": plot_cell_type_comparison(
            tabula_h5ad=args.tabula_h5ad,
            hbca_h5ad=args.hbca_h5ad,
            output_path=args.cell_type_out,
            max_points=args.max_points,
            top_n_cell_types=args.top_n_cell_types,
            seed=args.seed,
            point_size=args.point_size,
            alpha=args.alpha,
            dpi=args.dpi,
        ),
    }
    adora_summary = plot_adora_comparison(
        tabula_h5ad=args.tabula_h5ad,
        hbca_h5ad=args.hbca_h5ad,
        tabula_expression_cache=args.tabula_expression_cache,
        hbca_expression_cache=args.hbca_expression_cache,
        output_path=args.adora_out,
        quantile=args.adora_quantile,
        dpi=args.dpi,
    )
    summary["adora_high_comparison"] = adora_summary
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {args.summary}")


if __name__ == "__main__":
    main()
