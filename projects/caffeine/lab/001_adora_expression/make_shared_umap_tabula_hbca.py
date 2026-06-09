#!/usr/bin/env python
"""Build a sampled shared UMAP for Tabula Sapiens and HBCA from shared genes.

This is an intentionally practical integration pass: it uses raw/count-like
matrices, shared gene symbols, Scanpy-style library-size normalization,
TruncatedSVD, simple dataset centering in latent space, and UMAP.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import umap
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_TABULA_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_HBCA_H5AD = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_all_non_neuronal_b165f033.h5ad"
DEFAULT_TABULA_ADORA = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
DEFAULT_HBCA_ADORA = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_adora_expression.npz"
DEFAULT_OUT_DIR = LAB_DIR / "figures" / "shared_umap_tabula_hbca"
ADORA_GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")
ADORA_COLORS = {
    "ADORA1": "#1b9e77",
    "ADORA2A": "#d95f02",
    "ADORA2B": "#7570b3",
    "ADORA3": "#e7298a",
}


def decode_array(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x) for x in values]


def read_categorical_values(parent: h5py.Group, name: str) -> np.ndarray:
    node = parent[name]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = node["codes"][:]
        categories = np.asarray(decode_array(node["categories"][:]), dtype=object)
        out = np.full(codes.shape, "NA", dtype=object)
        mask = codes >= 0
        out[mask] = categories[codes[mask]]
        return out
    values = node[:]
    return np.asarray(decode_array(values), dtype=object)


def read_feature_names(f: h5py.File, var_path: str) -> np.ndarray:
    return read_categorical_values(f[var_path], "feature_name")


def canonical_gene(name: str) -> str:
    return name.strip().upper()


def first_gene_indices(names: np.ndarray) -> dict[str, int]:
    out: dict[str, int] = {}
    for idx, name in enumerate(names):
        key = canonical_gene(str(name))
        if key and key not in out:
            out[key] = idx
    return out


def read_bool_or_default(f: h5py.File, path: str, n: int) -> np.ndarray:
    if path not in f:
        return np.zeros(n, dtype=bool)
    node = f[path]
    if isinstance(node, h5py.Group):
        vals = read_categorical_values(f["/".join(path.split("/")[:-1])], path.split("/")[-1])
        return np.asarray([str(x).lower() == "true" for x in vals], dtype=bool)
    return node[:].astype(bool)


def select_shared_genes(
    tabula: h5py.File,
    hbca: h5py.File,
    n_genes: int,
) -> tuple[list[str], np.ndarray, np.ndarray, int]:
    tabula_names = read_feature_names(tabula, "raw/var")
    hbca_names = read_feature_names(hbca, "var")
    tabula_idx_by_gene = first_gene_indices(tabula_names)
    hbca_idx_by_gene = first_gene_indices(hbca_names)
    shared = sorted(set(tabula_idx_by_gene).intersection(hbca_idx_by_gene))

    tabula_filtered = read_bool_or_default(tabula, "raw/var/feature_is_filtered", len(tabula_names))
    tabula_std = tabula["raw/var/std"][:]
    scores = []
    for gene in shared:
        idx = tabula_idx_by_gene[gene]
        if tabula_filtered[idx]:
            continue
        if gene.startswith(("MT-", "RPL", "RPS", "ERCC-")):
            continue
        score = tabula_std[idx]
        if np.isfinite(score):
            scores.append((float(score), gene))

    ranked = [gene for _, gene in sorted(scores, reverse=True)]
    selected = ranked[:n_genes]
    for gene in ADORA_GENES:
        if gene in tabula_idx_by_gene and gene in hbca_idx_by_gene and gene not in selected:
            selected.append(gene)
    selected = sorted(selected)
    tabula_cols = np.asarray([tabula_idx_by_gene[g] for g in selected], dtype=np.int64)
    hbca_cols = np.asarray([hbca_idx_by_gene[g] for g in selected], dtype=np.int64)
    return selected, tabula_cols, hbca_cols, len(shared)


def stratified_sample(codes: np.ndarray, n_total: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_obs = codes.shape[0]
    if n_total >= n_obs:
        return np.arange(n_obs, dtype=np.int64)
    valid_codes = np.unique(codes[codes >= 0])
    if valid_codes.size == 0:
        return np.sort(rng.choice(n_obs, n_total, replace=False))

    per_group = max(1, int(np.ceil(n_total / valid_codes.size)))
    chosen: list[np.ndarray] = []
    chosen_mask = np.zeros(n_obs, dtype=bool)
    for code in valid_codes:
        members = np.flatnonzero(codes == code)
        take = min(per_group, members.size)
        picked = rng.choice(members, take, replace=False)
        chosen.append(picked)
        chosen_mask[picked] = True

    initial = np.concatenate(chosen) if chosen else np.array([], dtype=np.int64)
    if initial.size > n_total:
        initial = rng.choice(initial, n_total, replace=False)
    elif initial.size < n_total:
        remaining = np.flatnonzero(~chosen_mask)
        fill = rng.choice(remaining, n_total - initial.size, replace=False)
        initial = np.concatenate([initial, fill])
    return np.sort(initial.astype(np.int64))


def csr_subset_rows_cols_with_totals(
    matrix_group: h5py.Group,
    row_idx: np.ndarray,
    selected_cols: np.ndarray,
    n_selected: int,
    progress_label: str,
) -> tuple[sp.csr_matrix, np.ndarray]:
    shape = tuple(int(x) for x in matrix_group.attrs["shape"])
    indptr = matrix_group["indptr"][:]
    data_ds = matrix_group["data"]
    indices_ds = matrix_group["indices"]
    remap = np.full(shape[1], -1, dtype=np.int32)
    remap[selected_cols] = np.arange(n_selected, dtype=np.int32)

    out_indptr = np.zeros(row_idx.size + 1, dtype=np.int64)
    out_indices_parts: list[np.ndarray] = []
    out_data_parts: list[np.ndarray] = []
    totals = np.zeros(row_idx.size, dtype=np.float64)

    for out_row, source_row in enumerate(row_idx):
        start = int(indptr[source_row])
        end = int(indptr[source_row + 1])
        row_data = data_ds[start:end]
        row_indices = indices_ds[start:end]
        totals[out_row] = float(np.sum(row_data))
        mapped = remap[row_indices]
        keep = mapped >= 0
        if np.any(keep):
            cols = mapped[keep].astype(np.int32, copy=False)
            vals = row_data[keep].astype(np.float32, copy=False)
            order = np.argsort(cols)
            out_indices_parts.append(cols[order])
            out_data_parts.append(vals[order])
            out_indptr[out_row + 1] = out_indptr[out_row] + int(np.count_nonzero(keep))
        else:
            out_indptr[out_row + 1] = out_indptr[out_row]
        if (out_row + 1) % 10000 == 0:
            print(f"{progress_label}: extracted {out_row + 1:,}/{row_idx.size:,} rows", flush=True)

    out_indices = (
        np.concatenate(out_indices_parts).astype(np.int32, copy=False)
        if out_indices_parts
        else np.array([], dtype=np.int32)
    )
    out_data = (
        np.concatenate(out_data_parts).astype(np.float32, copy=False)
        if out_data_parts
        else np.array([], dtype=np.float32)
    )
    matrix = sp.csr_matrix((out_data, out_indices, out_indptr), shape=(row_idx.size, n_selected))
    return matrix, totals


def normalize_log1p_from_totals(x: sp.csr_matrix, totals: np.ndarray, target_sum: float) -> sp.csr_matrix:
    x = x.tocsr(copy=True)
    good = totals > 0
    scale = np.zeros_like(totals, dtype=np.float32)
    scale[good] = (target_sum / totals[good]).astype(np.float32)
    x = sp.diags(scale).dot(x).tocsr()
    x.data = np.log1p(x.data)
    return x


def read_obs_frame(f: h5py.File, columns: list[str], idx: np.ndarray) -> pd.DataFrame:
    data: dict[str, np.ndarray] = {}
    for col in columns:
        if f"obs/{col}" in f:
            data[col] = read_categorical_values(f["obs"], col)[idx]
        else:
            data[col] = np.full(idx.shape, "NA", dtype=object)
    return pd.DataFrame(data)


def assign_adora(expr: np.ndarray, thresholds: dict[str, float]) -> np.ndarray:
    high = np.zeros(expr.shape, dtype=bool)
    for i, gene in enumerate(ADORA_GENES):
        high[:, i] = expr[:, i] >= thresholds[gene]
    masked = np.where(high, expr, -np.inf)
    best = np.argmax(masked, axis=1)
    labels = np.full(expr.shape[0], "none", dtype=object)
    any_high = np.any(high, axis=1)
    for i, gene in enumerate(ADORA_GENES):
        labels[any_high & (best == i)] = gene
    return labels


def correct_latent(latent: np.ndarray, datasets: np.ndarray, method: str) -> tuple[np.ndarray, str]:
    if method == "center":
        corrected = latent.copy()
        for dataset in np.unique(datasets):
            mask = datasets == dataset
            corrected[mask] -= corrected[mask].mean(axis=0, keepdims=True)
        return corrected, "subtract dataset mean from each SVD component before UMAP"
    if method == "none":
        return latent.copy(), "none"
    if method == "harmony":
        try:
            import harmonypy as hm
        except ImportError as exc:
            raise SystemExit(
                "Harmony correction requested, but harmonypy is not installed. "
                "Run with `uv run --with harmonypy python ... --batch-correction harmony`."
            ) from exc
        metadata = pd.DataFrame({"dataset": datasets})
        harmony = hm.run_harmony(latent, metadata, vars_use=["dataset"], random_state=0)
        z_corr = np.asarray(harmony.Z_corr)
        if z_corr.shape == latent.shape:
            corrected = z_corr
        elif z_corr.T.shape == latent.shape:
            corrected = z_corr.T
        else:
            raise ValueError(f"Unexpected Harmony output shape {z_corr.shape}; expected {latent.shape}")
        return corrected.astype(np.float32), "Harmony correction on SVD components using dataset as batch"
    raise ValueError(f"Unknown batch correction method: {method}")


def plot_dataset(embedding: np.ndarray, dataset: np.ndarray, output: Path, point_size: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 7.6))
    colors = {"Tabula Sapiens": "#2166ac", "HBCA": "#b2182b"}
    for name, color in colors.items():
        mask = dataset == name
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            s=point_size,
            c=color,
            alpha=0.28,
            linewidths=0,
            label=f"{name} ({np.count_nonzero(mask):,})",
            rasterized=True,
        )
    ax.set_title("Shared UMAP from common genes")
    ax.set_xlabel("shared UMAP 1")
    ax.set_ylabel("shared UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, markerscale=3)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_adora(embedding: np.ndarray, labels: np.ndarray, output: Path, point_size: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 7.6))
    none = labels == "none"
    ax.scatter(
        embedding[none, 0],
        embedding[none, 1],
        s=point_size * 0.65,
        c="#d3d3d3",
        alpha=0.08,
        linewidths=0,
        label=f"not high ({np.count_nonzero(none):,})",
        rasterized=True,
    )
    for gene in ADORA_GENES:
        mask = labels == gene
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            s=point_size * 1.5,
            c=ADORA_COLORS[gene],
            alpha=0.65,
            linewidths=0,
            label=f"{gene} high ({np.count_nonzero(mask):,})",
            rasterized=True,
        )
    ax.set_title("ADORA-high cells in the shared UMAP")
    ax.set_xlabel("shared UMAP 1")
    ax.set_ylabel("shared UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tabula-h5ad", type=Path, default=DEFAULT_TABULA_H5AD)
    parser.add_argument("--hbca-h5ad", type=Path, default=DEFAULT_HBCA_H5AD)
    parser.add_argument("--tabula-adora-cache", type=Path, default=DEFAULT_TABULA_ADORA)
    parser.add_argument("--hbca-adora-cache", type=Path, default=DEFAULT_HBCA_ADORA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--cells-per-dataset", type=int, default=40000)
    parser.add_argument("--n-genes", type=int, default=2500)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--target-sum", type=float, default=10000.0)
    parser.add_argument("--batch-correction", choices=("center", "harmony", "none"), default="center")
    parser.add_argument("--seed", type=int, default=19)
    parser.add_argument("--point-size", type=float, default=2.0)
    parser.add_argument("--dpi", type=int, default=220)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.tabula_h5ad, "r") as tabula, h5py.File(args.hbca_h5ad, "r") as hbca:
        selected_genes, tabula_cols, hbca_cols, n_shared = select_shared_genes(tabula, hbca, args.n_genes)
        print(f"Selected {len(selected_genes):,} genes from {n_shared:,} shared gene symbols", flush=True)

        tabula_codes = tabula["obs/cell_type/codes"][:]
        hbca_codes = hbca["obs/supercluster_term/codes"][:]
        tabula_idx = stratified_sample(tabula_codes, args.cells_per_dataset, args.seed)
        hbca_idx = stratified_sample(hbca_codes, args.cells_per_dataset, args.seed + 1)
        print(f"Sampled {tabula_idx.size:,} Tabula and {hbca_idx.size:,} HBCA cells", flush=True)

        tabula_x, tabula_totals = csr_subset_rows_cols_with_totals(
            tabula["raw/X"], tabula_idx, tabula_cols, len(selected_genes), "Tabula raw/X"
        )
        hbca_x, hbca_totals = csr_subset_rows_cols_with_totals(
            hbca["X"], hbca_idx, hbca_cols, len(selected_genes), "HBCA X"
        )
        tabula_obs = read_obs_frame(
            tabula,
            ["cell_type", "broad_cell_class", "tissue", "tissue_in_publication"],
            tabula_idx,
        )
        hbca_obs = read_obs_frame(
            hbca,
            ["cell_type", "supercluster_term", "tissue", "ROIGroupCoarse", "ROIGroupFine", "roi"],
            hbca_idx,
        )

    print("Normalizing, log-transforming, and fitting shared SVD", flush=True)
    tabula_norm = normalize_log1p_from_totals(tabula_x, tabula_totals, args.target_sum)
    hbca_norm = normalize_log1p_from_totals(hbca_x, hbca_totals, args.target_sum)
    combined = sp.vstack([tabula_norm, hbca_norm], format="csr")
    svd = TruncatedSVD(n_components=args.n_components, random_state=args.seed)
    latent = svd.fit_transform(combined)
    latent = StandardScaler().fit_transform(latent)
    datasets = np.asarray(["Tabula Sapiens"] * tabula_idx.size + ["HBCA"] * hbca_idx.size, dtype=object)
    corrected, correction_description = correct_latent(latent, datasets, args.batch_correction)

    print("Computing shared UMAP", flush=True)
    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.5,
        metric="cosine",
        random_state=args.seed,
        low_memory=True,
        verbose=True,
    )
    embedding = reducer.fit_transform(corrected)

    print("Adding ADORA labels and writing outputs", flush=True)
    tabula_adora = np.load(args.tabula_adora_cache)["expression"][tabula_idx]
    hbca_adora = np.load(args.hbca_adora_cache)["expression"][hbca_idx]
    tabula_thresholds = {
        "ADORA1": 0.684346,
        "ADORA2A": 1.080287,
        "ADORA2B": 0.779279,
        "ADORA3": 0.943180,
    }
    hbca_thresholds = {gene: 1.0 for gene in ADORA_GENES}
    adora_labels = np.concatenate(
        [assign_adora(tabula_adora, tabula_thresholds), assign_adora(hbca_adora, hbca_thresholds)]
    )

    meta = pd.concat(
        [
            pd.DataFrame({"dataset": "Tabula Sapiens", "source_index": tabula_idx}).join(tabula_obs),
            pd.DataFrame({"dataset": "HBCA", "source_index": hbca_idx}).join(hbca_obs),
        ],
        ignore_index=True,
    )
    meta["shared_umap_1"] = embedding[:, 0]
    meta["shared_umap_2"] = embedding[:, 1]
    meta["adora_high_gene"] = adora_labels
    for i, gene in enumerate(ADORA_GENES):
        meta[f"{gene}_expression"] = np.concatenate([tabula_adora[:, i], hbca_adora[:, i]])

    coords_path = args.out_dir / "shared_umap_tabula_hbca_cells.csv"
    genes_path = args.out_dir / "shared_umap_tabula_hbca_genes.txt"
    npz_path = args.out_dir / "shared_umap_tabula_hbca_arrays.npz"
    summary_path = args.out_dir / "shared_umap_tabula_hbca_summary.json"
    meta.to_csv(coords_path, index=False)
    genes_path.write_text("\n".join(selected_genes) + "\n")
    np.savez_compressed(
        npz_path,
        embedding=embedding,
        latent=latent,
        latent_dataset_centered=corrected,
        dataset=datasets,
        tabula_indices=tabula_idx,
        hbca_indices=hbca_idx,
        genes=np.asarray(selected_genes, dtype=object),
        adora_high_gene=adora_labels,
    )

    dataset_fig = args.out_dir / "shared_umap_tabula_hbca_dataset.png"
    adora_fig = args.out_dir / "shared_umap_tabula_hbca_adora_high.png"
    plot_dataset(embedding, datasets, dataset_fig, args.point_size, args.dpi)
    plot_adora(embedding, adora_labels, adora_fig, args.point_size, args.dpi)

    summary = {
        "tabula_h5ad": str(args.tabula_h5ad),
        "hbca_h5ad": str(args.hbca_h5ad),
        "tabula_matrix": "raw/X",
        "hbca_matrix": "X",
        "cells_per_dataset_requested": args.cells_per_dataset,
        "tabula_cells": int(tabula_idx.size),
        "hbca_cells": int(hbca_idx.size),
        "shared_gene_symbols": int(n_shared),
        "selected_genes": int(len(selected_genes)),
        "normalization": f"counts per cell scaled to {args.target_sum:g}, then log1p",
        "latent": {
            "method": "TruncatedSVD on concatenated normalized shared-gene matrix",
            "n_components": args.n_components,
            "explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
            "batch_correction": correction_description,
        },
        "umap": {
            "n_neighbors": 15,
            "min_dist": 0.5,
            "metric": "cosine",
            "random_state": args.seed,
        },
        "original_umap_provenance": {
            "tabula": {
                "neighbors_method": "umap",
                "neighbors_metric": "euclidean",
                "neighbors_n_neighbors": 15,
                "neighbors_use_rep": "X_scvi",
                "umap_a": 0.583030019901822,
                "umap_b": 1.3341669931033755,
            },
            "hbca": "H5AD stores obsm/X_UMAP but no uns/umap or uns/neighbors parameters.",
        },
        "adora_thresholds": {
            "tabula": tabula_thresholds,
            "hbca": hbca_thresholds,
        },
        "outputs": {
            "cells_csv": str(coords_path),
            "arrays_npz": str(npz_path),
            "genes_txt": str(genes_path),
            "dataset_figure": str(dataset_fig),
            "adora_figure": str(adora_fig),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
