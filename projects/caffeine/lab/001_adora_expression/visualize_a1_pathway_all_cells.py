#!/usr/bin/env python
"""Score A1-pathway expression in every atlas cell and visualize native UMAPs.

The all-cell raster panels include every cell. Interactive Plotly views use a
stratified sample so they remain responsive on laptops and mobile browsers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.sparse as sp

from score_a1_branches import BRANCHES


LAB = Path(__file__).resolve().parent
OUT = LAB / "figures" / "a1_pathway_all_cells"

ATLASES = {
    "tabula_sapiens": {
        "path": LAB / "cache" / "tabula_sapiens_all_cells.h5ad",
        "umap": "X_umap",
        "cell_type": "cell_type",
        "tissue": "tissue_in_publication",
        "counts": False,
        "label": "Tabula Sapiens",
    },
    "hbca_neurons": {
        "path": LAB / "cache" / "human_brain_cell_atlas" / "hbca_all_neurons_8e10f1c4.h5ad",
        "umap": "X_UMAP",
        "cell_type": "supercluster_term",
        "tissue": "roi",
        "counts": True,
        "label": "Human Brain Cell Atlas neurons",
    },
}


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def categorical(group: h5py.Group, name: str) -> np.ndarray:
    node = group[name]
    if isinstance(node, h5py.Group):
        categories = decode(node["categories"][:])
        codes = node["codes"][:]
        result = np.full(codes.shape, "unknown", dtype=object)
        valid = codes >= 0
        result[valid] = categories[codes[valid]]
        return result
    return decode(node[:])


def feature_names(handle: h5py.File) -> np.ndarray:
    node = handle["var/feature_name"]
    if isinstance(node, h5py.Group):
        return decode(node["categories"][:])[node["codes"][:]]
    return decode(node[:])


def score_cells(handle: h5py.File, count_matrix: bool, chunk_size: int) -> tuple[np.ndarray, list[str]]:
    names = feature_names(handle)
    first_index = {}
    for index, name in enumerate(names):
        first_index.setdefault(str(name).upper(), index)

    found = {branch: [first_index[g] for g in genes if g in first_index] for branch, genes in BRANCHES.items()}
    missing = [g for genes in BRANCHES.values() for g in genes if g not in first_index]
    columns = np.asarray(sorted({index for indices in found.values() for index in indices}), dtype=np.int64)
    local_index = {column: index for index, column in enumerate(columns)}
    branch_local = {branch: np.asarray([local_index[index] for index in indices]) for branch, indices in found.items() if indices}

    matrix = handle["X"]
    indptr = matrix["indptr"][:]
    n_cells, n_genes = (int(v) for v in matrix.attrs["shape"])
    scores = np.zeros(n_cells, dtype=np.float32)
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        data_start, data_end = int(indptr[start]), int(indptr[end])
        block = sp.csr_matrix(
            (matrix["data"][data_start:data_end], matrix["indices"][data_start:data_end], indptr[start:end + 1] - data_start),
            shape=(end - start, n_genes),
        )
        selected = block[:, columns].astype(np.float32)
        if count_matrix:
            totals = np.asarray(block.sum(axis=1)).ravel()
            scale = np.divide(10_000.0, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
            selected = sp.diags(scale).dot(selected).tocsr()
            selected.data = np.log1p(selected.data)
        branch_scores = np.column_stack(
            [np.asarray(selected[:, indices].mean(axis=1)).ravel() for indices in branch_local.values()]
        )
        scores[start:end] = branch_scores.mean(axis=1)
        print(f"  scored {end:,}/{n_cells:,}", flush=True)
    return scores, missing


def stratified_sample(labels: np.ndarray, size: int, seed: int) -> np.ndarray:
    if size >= len(labels):
        return np.arange(len(labels))
    rng = np.random.default_rng(seed)
    selected = []
    groups, counts = np.unique(labels, return_counts=True)
    allocation = np.maximum(15, np.floor(size * counts / counts.sum()).astype(int))
    for group, take in zip(groups, allocation, strict=True):
        members = np.flatnonzero(labels == group)
        selected.append(rng.choice(members, min(take, len(members)), replace=False))
    merged = np.concatenate(selected)
    if len(merged) > size:
        merged = rng.choice(merged, size, replace=False)
    elif len(merged) < size:
        remaining = np.setdiff1d(np.arange(len(labels)), merged, assume_unique=False)
        merged = np.concatenate([merged, rng.choice(remaining, size - len(merged), replace=False)])
    return np.sort(merged)


def all_cell_raster(coords: np.ndarray, scores: np.ndarray, label: str, output: Path) -> None:
    xlow, xhigh = np.quantile(coords[:, 0], [.001, .999])
    ylow, yhigh = np.quantile(coords[:, 1], [.001, .999])
    keep = (coords[:, 0] >= xlow) & (coords[:, 0] <= xhigh) & (coords[:, 1] >= ylow) & (coords[:, 1] <= yhigh)
    count, xedges, yedges = np.histogram2d(coords[keep, 0], coords[keep, 1], bins=1100, range=((xlow, xhigh), (ylow, yhigh)))
    weighted, _, _ = np.histogram2d(coords[keep, 0], coords[keep, 1], bins=1100, range=((xlow, xhigh), (ylow, yhigh)), weights=scores[keep])
    mean = np.divide(weighted, count, out=np.full_like(weighted, np.nan), where=count > 0)
    vmax = float(np.nanquantile(mean, .995))
    alpha = np.clip(np.log1p(count) / np.log(12), 0, 1)
    rgba = plt.get_cmap("magma")(np.clip(mean / max(vmax, 1e-8), 0, 1))
    rgba[..., 3] = np.where(np.isfinite(mean), alpha, 0)
    fig, ax = plt.subplots(figsize=(11, 9), facecolor="#08080c")
    ax.set_facecolor("#08080c")
    ax.imshow(rgba.transpose(1, 0, 2), origin="lower", extent=(xlow, xhigh, ylow, yhigh), aspect="auto")
    ax.set_title(f"{label}: A1 pathway expression across all {len(coords):,} cells", color="white", pad=12)
    ax.text(.01, .015, "Pixel color = mean branch-balanced pathway expression; opacity = cell density", transform=ax.transAxes, color="#dddddd", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.tight_layout(); fig.savefig(output, dpi=180, facecolor=fig.get_facecolor()); plt.close(fig)


def interactive(coords: np.ndarray, scores: np.ndarray, cell_types: np.ndarray, tissues: np.ndarray, label: str, output: Path, sample_size: int) -> None:
    idx = stratified_sample(cell_types, sample_size, seed=23)
    custom = np.column_stack([idx, cell_types[idx], tissues[idx], scores[idx]])
    fig = go.Figure(go.Scattergl(
        x=coords[idx, 0], y=coords[idx, 1], mode="markers",
        marker={"size": 3, "color": scores[idx], "colorscale": "Magma", "opacity": .68, "colorbar": {"title": "A1 pathway"}},
        customdata=custom,
        hovertemplate="<b>%{customdata[1]}</b><br>tissue/region: %{customdata[2]}<br>pathway score: %{customdata[3]:.4f}<br>cell index: %{customdata[0]}<extra></extra>",
    ))
    fig.update_layout(
        title=f"{label}: A1 pathway on native UMAP (interactive stratified sample {len(idx):,}/{len(coords):,})",
        template="plotly_dark", paper_bgcolor="#08080c", plot_bgcolor="#08080c", dragmode="pan",
        margin={"l": 20, "r": 20, "t": 65, "b": 20},
        xaxis={"visible": False}, yaxis={"visible": False, "scaleanchor": "x", "scaleratio": 1},
    )
    fig.write_html(output, include_plotlyjs=True, full_html=True, config={"responsive": True, "scrollZoom": True})


def summaries(scores: np.ndarray, cell_types: np.ndarray, tissues: np.ndarray, atlas: str, output: Path) -> None:
    frame = pd.DataFrame({"atlas": atlas, "cell_type": cell_types, "tissue": tissues, "score": scores})
    rows = []
    for dimension in ("cell_type", "tissue"):
        table = frame.groupby(dimension, observed=True)["score"].agg(["size", "mean", "median"])
        table["positive_fraction"] = frame.assign(positive=frame.score > np.quantile(scores, .9)).groupby(dimension, observed=True).positive.mean()
        table = table.reset_index().rename(columns={dimension: "group"})
        table.insert(1, "dimension", dimension)
        rows.append(table)
    pd.concat(rows, ignore_index=True).sort_values(["dimension", "mean"], ascending=[True, False]).to_csv(output, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", choices=(*ATLASES, "all"), default="all")
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--interactive-points", type=int, default=100_000)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    chosen = ATLASES if args.atlas == "all" else {args.atlas: ATLASES[args.atlas]}
    manifest = {}
    for slug, config in chosen.items():
        print(f"\n{config['label']}", flush=True)
        with h5py.File(config["path"], "r") as handle:
            scores, missing = score_cells(handle, config["counts"], args.chunk_size)
            coords = handle[f"obsm/{config['umap']}"][:].astype(np.float32)
            cell_types = categorical(handle["obs"], config["cell_type"])
            tissues = categorical(handle["obs"], config["tissue"])
        np.savez_compressed(OUT / f"{slug}_a1_pathway_scores.npz", score=scores)
        all_cell_raster(coords, scores, config["label"], OUT / f"{slug}_all_cells.png")
        interactive(coords, scores, cell_types, tissues, config["label"], OUT / f"{slug}_interactive.html", args.interactive_points)
        summaries(scores, cell_types, tissues, config["label"], OUT / f"{slug}_group_summary.csv")
        manifest[slug] = {"n_cells": len(scores), "missing_genes": missing, "score_p90": float(np.quantile(scores, .9))}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
