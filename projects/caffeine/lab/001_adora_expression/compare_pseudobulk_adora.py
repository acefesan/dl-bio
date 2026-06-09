#!/usr/bin/env python
"""Compare pseudobulk ADORA receptor summaries in Tabula Sapiens and HBCA."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
FIGURES_DIR = LAB_DIR / "figures"
TABULA_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
TABULA_ADORA_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
HBCA_H5AD = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_all_non_neuronal_b165f033.h5ad"
HBCA_ADORA_CACHE = (
    LAB_DIR
    / "cache"
    / "human_brain_cell_atlas"
    / "compare_pseudobulk_hbca_adora_expression.npz"
)

GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")


def decode_value(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def decode_array(values: np.ndarray) -> list[str]:
    return [decode_value(value) for value in values]


def read_categorical(h5ad_path: Path, obs_column: str) -> tuple[np.ndarray, list[str]]:
    with h5py.File(h5ad_path, "r") as f:
        node = f["obs"].get(obs_column)
        if node is None:
            raise KeyError(f"{h5ad_path} has no obs/{obs_column!r}")
        if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
            raise ValueError(f"obs/{obs_column!r} is not an AnnData categorical column")
        return node["codes"][:], decode_array(node["categories"][:])


def read_gene_names(h5ad_path: Path) -> list[str]:
    with h5py.File(h5ad_path, "r") as f:
        feature_name = f["var"].get("feature_name")
        if isinstance(feature_name, h5py.Group) and {"codes", "categories"} <= set(feature_name.keys()):
            codes = feature_name["codes"][:]
            categories = decode_array(feature_name["categories"][:])
            return [categories[code] if code >= 0 else "" for code in codes]
        if isinstance(feature_name, h5py.Dataset):
            return decode_array(feature_name[:])
        return decode_array(f["var"]["_index"][:])


def find_gene_indices(h5ad_path: Path, genes: tuple[str, ...]) -> np.ndarray:
    gene_names = np.array(read_gene_names(h5ad_path), dtype=object)
    indices = []
    for gene in genes:
        matches = np.flatnonzero(gene_names == gene)
        if len(matches) == 0:
            raise KeyError(f"Could not find {gene!r} in {h5ad_path}")
        indices.append(int(matches[0]))
    return np.array(indices, dtype=np.int64)


def extract_csr_columns(h5ad_path: Path, gene_indices: np.ndarray, chunk_rows: int) -> np.ndarray:
    with h5py.File(h5ad_path, "r") as f:
        x = f["X"]
        shape = tuple(int(v) for v in x.attrs["shape"])
        if x.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"Expected csr_matrix X in {h5ad_path}")

        expression = np.zeros((shape[0], len(gene_indices)), dtype=np.float32)
        col_to_gene = {int(col): i for i, col in enumerate(gene_indices)}
        indptr = x["indptr"]
        indices = x["indices"]
        data = x["data"]

        for start_row in range(0, shape[0], chunk_rows):
            end_row = min(start_row + chunk_rows, shape[0])
            ptr = indptr[start_row : end_row + 1]
            start_ptr = int(ptr[0])
            end_ptr = int(ptr[-1])
            chunk_indices = indices[start_ptr:end_ptr]
            chunk_data = data[start_ptr:end_ptr]

            for row_offset in range(end_row - start_row):
                row_start = int(ptr[row_offset] - start_ptr)
                row_end = int(ptr[row_offset + 1] - start_ptr)
                row_indices = chunk_indices[row_start:row_end]
                row_data = chunk_data[row_start:row_end]
                for col, value in zip(row_indices, row_data, strict=False):
                    target = col_to_gene.get(int(col))
                    if target is not None:
                        expression[start_row + row_offset, target] = value
    return expression


def load_or_create_hbca_expression_cache(h5ad_path: Path, cache_path: Path, chunk_rows: int) -> np.ndarray:
    if cache_path.exists():
        cached = np.load(cache_path, allow_pickle=False)
        genes = tuple(decode_array(cached["genes"]))
        if genes != GENES:
            raise ValueError(f"Unexpected genes in {cache_path}: {genes}")
        return cached["expression"]

    gene_indices = find_gene_indices(h5ad_path, GENES)
    expression = extract_csr_columns(h5ad_path, gene_indices, chunk_rows)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        expression=expression,
        gene_indices=gene_indices,
        genes=np.array(GENES),
        matrix_group=np.array("X"),
    )
    print(f"Wrote {cache_path}")
    return expression


def load_tabula_expression(cache_path: Path) -> np.ndarray:
    cached = np.load(cache_path, allow_pickle=False)
    genes = tuple(decode_array(cached["genes"]))
    if genes != GENES:
        raise ValueError(f"Unexpected genes in {cache_path}: {genes}")
    return cached["expression"]


def summarize_by_group(
    dataset: str,
    source_label: str,
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
                    "dataset": dataset,
                    "source_label": source_label,
                    "group": group,
                    "group_normalized": group.casefold(),
                    "gene": gene,
                    "n_cells": n_cells,
                    "pct_expressing": float(pct_expr[i]),
                    "mean_expression": float(mean_expr[i]),
                }
            )
    return pd.DataFrame(rows)


def build_all_summaries(
    tabula_expression: np.ndarray,
    hbca_expression: np.ndarray,
    min_cells: int,
) -> pd.DataFrame:
    specs = [
        ("Tabula Sapiens", TABULA_H5AD, tabula_expression, "cell_type"),
        ("Tabula Sapiens", TABULA_H5AD, tabula_expression, "broad_cell_class"),
        ("Tabula Sapiens", TABULA_H5AD, tabula_expression, "tissue_in_publication"),
        ("HBCA non-neuronal", HBCA_H5AD, hbca_expression, "cell_type"),
        ("HBCA non-neuronal", HBCA_H5AD, hbca_expression, "supercluster_term"),
        ("HBCA non-neuronal", HBCA_H5AD, hbca_expression, "tissue"),
        ("HBCA non-neuronal", HBCA_H5AD, hbca_expression, "roi"),
    ]
    summaries = []
    for dataset, h5ad_path, expression, source_label in specs:
        codes, categories = read_categorical(h5ad_path, source_label)
        if expression.shape[0] != codes.shape[0]:
            raise ValueError(
                f"{dataset} {source_label}: expression rows {expression.shape[0]} "
                f"do not match obs rows {codes.shape[0]}"
            )
        summaries.append(summarize_by_group(dataset, source_label, expression, codes, categories, min_cells))
    return pd.concat(summaries, ignore_index=True)


def compare_overlapping_cell_types(summary: pd.DataFrame) -> pd.DataFrame:
    cell_type = summary[summary["source_label"] == "cell_type"].copy()
    wide = cell_type.pivot_table(
        index=["group_normalized", "gene"],
        columns="dataset",
        values=["group", "n_cells", "pct_expressing", "mean_expression"],
        aggfunc="first",
    )
    wide.columns = [f"{metric}__{dataset}" for metric, dataset in wide.columns]
    wide = wide.reset_index()
    required = [
        "pct_expressing__Tabula Sapiens",
        "pct_expressing__HBCA non-neuronal",
        "mean_expression__Tabula Sapiens",
        "mean_expression__HBCA non-neuronal",
    ]
    overlap = wide.dropna(subset=required).copy()
    overlap["pct_expressing_delta_hbca_minus_tabula"] = (
        overlap["pct_expressing__HBCA non-neuronal"] - overlap["pct_expressing__Tabula Sapiens"]
    )
    overlap["mean_expression_delta_hbca_minus_tabula"] = (
        overlap["mean_expression__HBCA non-neuronal"] - overlap["mean_expression__Tabula Sapiens"]
    )
    overlap["abs_pct_expressing_delta"] = overlap["pct_expressing_delta_hbca_minus_tabula"].abs()
    overlap["abs_mean_expression_delta"] = overlap["mean_expression_delta_hbca_minus_tabula"].abs()
    return overlap.sort_values(["abs_pct_expressing_delta", "abs_mean_expression_delta"], ascending=False)


def rank_summaries(summary: pd.DataFrame) -> pd.DataFrame:
    ranked = summary.copy()
    ranked["score"] = ranked["pct_expressing"] * np.log1p(ranked["mean_expression"])
    ranked = ranked.sort_values(["dataset", "source_label", "gene", "score"], ascending=[True, True, True, False])
    ranked["rank_within_dataset_label_gene"] = (
        ranked.groupby(["dataset", "source_label", "gene"]).cumcount() + 1
    )
    return ranked


def plot_comparison(overlap: pd.DataFrame, output_path: Path) -> None:
    genes = list(GENES)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))

    for ax, gene in zip(axes.flat, genes, strict=True):
        gene_overlap = overlap[overlap["gene"] == gene].head(8).iloc[::-1]
        y = np.arange(len(gene_overlap))
        labels = gene_overlap["group__Tabula Sapiens"].fillna(gene_overlap["group_normalized"]).str.replace("_", " ")
        ax.barh(
            y - 0.18,
            gene_overlap["pct_expressing__Tabula Sapiens"],
            height=0.34,
            label="Tabula",
            color="#5b8fd9",
        )
        ax.barh(
            y + 0.18,
            gene_overlap["pct_expressing__HBCA non-neuronal"],
            height=0.34,
            label="HBCA",
            color="#c85a54",
        )
        ax.set_yticks(y, labels)
        ax.set_xlabel("% expressing")
        ax.set_title(f"{gene}: overlapping cell-type labels")
        ax.grid(axis="x", color="#e6e6e6", linewidth=0.7)
        ax.set_axisbelow(True)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.95))
    fig.suptitle("Largest ADORA prevalence differences in overlapping Tabula/HBCA cell types", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tabula-h5ad", type=Path, default=TABULA_H5AD)
    parser.add_argument("--tabula-expression-cache", type=Path, default=TABULA_ADORA_CACHE)
    parser.add_argument("--hbca-h5ad", type=Path, default=HBCA_H5AD)
    parser.add_argument("--hbca-expression-cache", type=Path, default=HBCA_ADORA_CACHE)
    parser.add_argument("--min-cells", type=int, default=50)
    parser.add_argument("--chunk-rows", type=int, default=25_000)
    parser.add_argument("--out-prefix", type=Path, default=FIGURES_DIR / "compare_pseudobulk_adora")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    global TABULA_H5AD, HBCA_H5AD
    TABULA_H5AD = args.tabula_h5ad
    HBCA_H5AD = args.hbca_h5ad

    tabula_expression = load_tabula_expression(args.tabula_expression_cache)
    hbca_expression = load_or_create_hbca_expression_cache(
        args.hbca_h5ad,
        args.hbca_expression_cache,
        args.chunk_rows,
    )

    summary = build_all_summaries(tabula_expression, hbca_expression, args.min_cells)
    overlap = compare_overlapping_cell_types(summary)
    ranked = rank_summaries(summary)

    all_table = args.out_prefix.with_name(f"{args.out_prefix.name}_all_summaries.csv")
    overlap_table = args.out_prefix.with_name(f"{args.out_prefix.name}_cell_type_overlap.csv")
    ranked_table = args.out_prefix.with_name(f"{args.out_prefix.name}_ranked.csv")
    hbca_supercluster_table = args.out_prefix.with_name(f"{args.out_prefix.name}_hbca_supercluster_top.csv")
    hbca_roi_table = args.out_prefix.with_name(f"{args.out_prefix.name}_hbca_roi_top.csv")
    figure = args.out_prefix.with_name(f"{args.out_prefix.name}_cell_type_overlap.png")

    all_table.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_values(["dataset", "source_label", "group", "gene"]).to_csv(all_table, index=False)
    overlap.to_csv(overlap_table, index=False)
    ranked.to_csv(ranked_table, index=False)
    ranked[
        (ranked["dataset"] == "HBCA non-neuronal")
        & (ranked["source_label"] == "supercluster_term")
        & (ranked["rank_within_dataset_label_gene"] <= 10)
    ].to_csv(hbca_supercluster_table, index=False)
    ranked[
        (ranked["dataset"] == "HBCA non-neuronal")
        & (ranked["source_label"] == "roi")
        & (ranked["rank_within_dataset_label_gene"] <= 15)
    ].to_csv(hbca_roi_table, index=False)
    plot_comparison(overlap, figure)

    print(f"Wrote {all_table}")
    print(f"Wrote {overlap_table}")
    print(f"Wrote {ranked_table}")
    print(f"Wrote {hbca_supercluster_table}")
    print(f"Wrote {hbca_roi_table}")
    print(f"Wrote {figure}")


if __name__ == "__main__":
    main()
