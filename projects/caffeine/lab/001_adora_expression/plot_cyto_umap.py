#!/usr/bin/env python
"""Plot Tabula Sapiens UMAP cells colored by dominant cytoskeleton gene family.

Each cell is colored by whichever family has the highest total raw expression.
Unlike the ADORA plot (which highlights a rare high-expressing minority), cytoskeleton
genes are ubiquitous, so this paints the entire atlas and shows the structural
identity of every region.

Gene families
-------------
Actin (ubiquitous)  ACTB, ACTG1                              — present in all cells
Actin (muscle)      ACTA1, ACTA2, ACTC1, ACTG2              — contractile tissue
Tubulin             TUBA1A/1B/1C, TUBB/2A/4B, TUBG1         — highway / centrosome
Neurofilaments      NEFL, NEFM, NEFH                         — long axons
Vimentin            VIM                                       — connective tissue
Keratins            KRT1/5/8/14/18/19                        — epithelial cells
"""

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
DEFAULT_EXPR_CACHE = LAB_DIR / "cache" / "tabula_sapiens_cyto_expression.npz"
DEFAULT_FIGURE = LAB_DIR / "figures" / "tabula_sapiens_cyto_umap.png"
DEFAULT_SUMMARY = LAB_DIR / "figures" / "tabula_sapiens_cyto_umap_summary.json"

GENE_GROUPS: dict[str, list[str]] = {
    "Actin (ubiquitous)": ["ACTB", "ACTG1"],
    "Actin (muscle)": ["ACTA1", "ACTA2", "ACTC1", "ACTG2"],
    "Tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBB", "TUBB2A", "TUBB4B", "TUBG1"],
    "Neurofilaments": ["NEFL", "NEFM", "NEFH"],
    "Vimentin": ["VIM"],
    "Keratins": ["KRT1", "KRT5", "KRT8", "KRT14", "KRT18", "KRT19"],
}

COLORS: dict[str, str] = {
    "Actin (ubiquitous)": "#2196F3",  # blue
    "Actin (muscle)": "#F44336",  # red
    "Tubulin": "#FF9800",  # orange
    "Neurofilaments": "#9C27B0",  # purple
    "Vimentin": "#009688",  # teal
    "Keratins": "#E91E63",  # pink
}

GROUPS = tuple(GENE_GROUPS.keys())


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_var_names(f: h5py.File) -> np.ndarray:
    node = f["var/feature_name"]
    if isinstance(node, h5py.Group):
        categories = decode_values(node["categories"][:])
        codes = node["codes"][:]
        return np.array([categories[code] if code >= 0 else "" for code in codes], dtype=object)
    return np.array(decode_values(node[:]), dtype=object)


def find_gene_indices(f: h5py.File, gene_groups: dict[str, list[str]]) -> dict[str, int]:
    """Return {gene: var_index} for all genes found in the dataset; warns on missing."""
    names = read_var_names(f)
    found: dict[str, int] = {}
    for genes in gene_groups.values():
        for gene in genes:
            matches = np.flatnonzero(names == gene)
            if len(matches) == 1:
                found[gene] = int(matches[0])
            elif len(matches) > 1:
                print(f"  Warning: {gene} has {len(matches)} matches — skipping")
            else:
                print(f"  Warning: {gene} not found in var/feature_name — skipping")
    return found


def extract_group_expression(
    h5ad_path: Path,
    matrix_group: str,
    gene_groups: dict[str, list[str]],
    chunk_nnz: int,
) -> tuple[np.ndarray, dict[str, dict[str, int]]]:
    """Return (n_obs, n_groups) summed expression and {group: {gene: var_idx}} index."""
    with h5py.File(h5ad_path, "r") as f:
        gene_idx_map = find_gene_indices(f, gene_groups)

        group_names = list(gene_groups.keys())
        target_var = []
        target_grp = []
        found_by_group: dict[str, dict[str, int]] = {g: {} for g in group_names}

        for gi, (group, genes) in enumerate(gene_groups.items()):
            for gene in genes:
                if gene in gene_idx_map:
                    target_var.append(gene_idx_map[gene])
                    target_grp.append(gi)
                    found_by_group[group][gene] = gene_idx_map[gene]

        target_var = np.array(target_var, dtype=np.int64)
        target_grp = np.array(target_grp, dtype=np.int64)

        matrix = f[matrix_group]
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{matrix_group} is not a CSR matrix")

        n_obs, _ = matrix.attrs["shape"]
        indptr = matrix["indptr"][:]
        indices_ds = matrix["indices"]
        data_ds = matrix["data"]
        expr = np.zeros((int(n_obs), len(group_names)), dtype=np.float32)

        order = np.argsort(target_var)
        sorted_var = target_var[order]

        n_nnz = int(indices_ds.shape[0])
        for start in range(0, n_nnz, chunk_nnz):
            stop = min(start + chunk_nnz, n_nnz)
            indices = indices_ds[start:stop]
            positions = np.searchsorted(sorted_var, indices)
            in_range = positions < len(sorted_var)
            mask = np.zeros(indices.shape, dtype=bool)
            mask[in_range] = sorted_var[positions[in_range]] == indices[in_range]
            if not np.any(mask):
                continue

            offsets = np.flatnonzero(mask)
            rows = np.searchsorted(indptr, start + offsets, side="right") - 1
            grp_cols = target_grp[order[positions[mask]]]
            values = data_ds[start:stop][mask].astype(np.float32)
            np.add.at(expr, (rows, grp_cols), values)
            print(f"  grouped {int(np.count_nonzero(expr.sum(axis=1))):,} cells after {stop:,}/{n_nnz:,} nnz")

    return expr, found_by_group


def load_or_extract(
    h5ad_path: Path,
    cache_path: Path,
    matrix_group: str,
    gene_groups: dict[str, list[str]],
    chunk_nnz: int,
    force: bool,
) -> tuple[np.ndarray, dict[str, dict[str, int]]]:
    if cache_path.exists() and not force:
        cached = np.load(cache_path, allow_pickle=True)
        if (
            tuple(cached["groups"].tolist()) == tuple(gene_groups.keys())
            and str(cached["matrix_group"]) == matrix_group
        ):
            found_by_group = {
                g: {gene: int(idx) for gene, idx in zip(genes, idxs)}
                for g, genes, idxs in zip(
                    cached["groups"].tolist(),
                    cached["found_genes"],
                    cached["found_indices"],
                )
            }
            return cached["expression"], found_by_group

    expr, found_by_group = extract_group_expression(h5ad_path, matrix_group, gene_groups, chunk_nnz)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    groups = list(found_by_group.keys())
    found_genes = [list(found_by_group[g].keys()) for g in groups]
    found_indices = [list(found_by_group[g].values()) for g in groups]
    max_len = max(len(x) for x in found_genes)
    found_genes_arr = np.array([x + [""] * (max_len - len(x)) for x in found_genes], dtype=object)
    found_indices_arr = np.array([x + [-1] * (max_len - len(x)) for x in found_indices], dtype=np.int64)

    np.savez_compressed(
        cache_path,
        expression=expr,
        groups=np.array(groups),
        matrix_group=np.array(matrix_group),
        found_genes=found_genes_arr,
        found_indices=found_indices_arr,
    )
    print(f"Wrote expression cache {cache_path}")
    return expr, found_by_group


def plot_umap(
    h5ad_path: Path,
    output_path: Path,
    expr: np.ndarray,
    assigned: np.ndarray,
    groups: tuple[str, ...],
    point_size: float,
    alpha: float,
    dpi: int,
) -> None:
    with h5py.File(h5ad_path, "r") as f:
        umap = f["obsm/X_umap"][:]

    fig, ax = plt.subplots(figsize=(9.0, 8.0))

    # Draw unassigned cells (all-zero expression) in gray
    no_expr = assigned == -1
    if np.any(no_expr):
        ax.scatter(
            umap[no_expr, 0], umap[no_expr, 1],
            c="#cccccc", s=point_size, alpha=alpha * 0.5, linewidths=0, rasterized=True,
        )

    # Draw each family in order of ascending median expression so rare families
    # render on top of dominant ones
    group_totals = [np.median(expr[assigned == i, i]) if np.any(assigned == i) else 0.0 for i in range(len(groups))]
    draw_order = np.argsort(group_totals)[::-1]  # highest median drawn last (on top)

    for i in draw_order:
        mask = assigned == i
        if not np.any(mask):
            continue
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            c=COLORS[groups[i]], s=point_size, alpha=alpha, linewidths=0, rasterized=True,
        )

    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=COLORS[g], markeredgecolor="none",
            label=f"{g} ({int(np.sum(assigned == i)):,} cells)",
        )
        for i, g in enumerate(groups)
        if np.any(assigned == i)
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=8)
    ax.set_title("Tabula Sapiens — dominant cytoskeleton gene family per cell", fontsize=12)
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
    expr: np.ndarray,
    assigned: np.ndarray,
    found_by_group: dict[str, dict[str, int]],
    groups: tuple[str, ...],
    matrix_group: str,
) -> None:
    summary = {
        "matrix_group": matrix_group,
        "assignment_rule": "winner-takes-all on summed group expression per cell",
        "n_cells_total": int(expr.shape[0]),
        "n_cells_assigned": int(np.sum(assigned >= 0)),
        "n_cells_unassigned": int(np.sum(assigned == -1)),
        "groups": {
            g: {
                "n_cells_dominant": int(np.sum(assigned == i)),
                "genes_found": list(found_by_group[g].keys()),
                "genes_missing": [
                    gene for gene in GENE_GROUPS[g] if gene not in found_by_group[g]
                ],
                "median_expr_among_dominant": float(
                    np.median(expr[assigned == i, i]) if np.any(assigned == i) else 0.0
                ),
            }
            for i, g in enumerate(groups)
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    p.add_argument("--matrix-group", default="X")
    p.add_argument("--expression-cache", type=Path, default=DEFAULT_EXPR_CACHE)
    p.add_argument("--out", type=Path, default=DEFAULT_FIGURE)
    p.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    p.add_argument("--chunk-nnz", type=int, default=25_000_000)
    p.add_argument("--force-extract", action="store_true")
    p.add_argument("--point-size", type=float, default=0.05)
    p.add_argument("--alpha", type=float, default=0.55)
    p.add_argument("--dpi", type=int, default=260)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    expr, found_by_group = load_or_extract(
        h5ad_path=args.h5ad,
        cache_path=args.expression_cache,
        matrix_group=args.matrix_group,
        gene_groups=GENE_GROUPS,
        chunk_nnz=args.chunk_nnz,
        force=args.force_extract,
    )

    # Winner-takes-all: each cell → group with highest summed expression
    # Cells with all-zero expression get assigned = -1
    any_expr = expr.sum(axis=1) > 0
    assigned = np.full(expr.shape[0], -1, dtype=np.int8)
    assigned[any_expr] = np.argmax(expr[any_expr], axis=1).astype(np.int8)

    plot_umap(
        h5ad_path=args.h5ad,
        output_path=args.out,
        expr=expr,
        assigned=assigned,
        groups=GROUPS,
        point_size=args.point_size,
        alpha=args.alpha,
        dpi=args.dpi,
    )
    write_summary(
        summary_path=args.summary,
        expr=expr,
        assigned=assigned,
        found_by_group=found_by_group,
        groups=GROUPS,
        matrix_group=args.matrix_group,
    )


if __name__ == "__main__":
    main()
