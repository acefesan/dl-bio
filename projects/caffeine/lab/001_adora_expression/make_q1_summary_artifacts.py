#!/usr/bin/env python
"""Build the Q1 summary artifacts that the lab entry.md promised.

Reads ADORA expression that was already extracted into compact npz files plus
obs metadata read directly from the downloaded H5ADs with h5py (no SOMA, no
network). Combines Tabula Sapiens (peripheral, 1.14M cells) and HBCA non-
neuronal (brain, 888k cells) so the pseudobulk table is body-wide.

Outputs:
    figures/ranked_top20_per_receptor.png
    figures/donor_stratified_dotplot.png            (Tabula only, donor_id)
    figures/assay_stratified_dotplot.png            (Tabula only, assay)
    cache/pseudobulk_by_cell_type.feather           (combined Tabula + HBCA)
    cache/cross_receptor_overlap.feather            (combined; per-cell-type)
    cache/q1_summary.json                           run metadata

What this script does NOT cover:
    HBCA neurons (30 GB H5AD not yet downloaded) — neurons are the canonical
    brain ADORA1/2A story. Anything striatal/cortical here is absent.
"""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GENES = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
LAB_DIR = Path(__file__).resolve().parent
CACHE_DIR = LAB_DIR / "cache"
FIG_DIR = LAB_DIR / "figures"

TABULA_H5AD = CACHE_DIR / "tabula_sapiens_all_cells.h5ad"
TABULA_EXPR = CACHE_DIR / "tabula_sapiens_adora_expression.npz"
HBCA_H5AD = CACHE_DIR / "human_brain_cell_atlas" / "hbca_all_non_neuronal_b165f033.h5ad"
HBCA_EXPR = CACHE_DIR / "human_brain_cell_atlas" / "hbca_adora_expression.npz"


def decode(values):
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def read_cat_column(group: h5py.Group) -> np.ndarray:
    """Read an AnnData-style categorical obs column."""
    if "categories" in group and "codes" in group:
        cats = decode(group["categories"][:])
        codes = group["codes"][:]
        return np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
    return np.array(decode(group[:]), dtype=object)


def load_obs(h5ad_path: Path, columns: list[str]) -> pd.DataFrame:
    """Load obs columns from an H5AD without instantiating an AnnData."""
    out = {}
    with h5py.File(h5ad_path, "r") as f:
        for col in columns:
            if col in f["obs"]:
                out[col] = read_cat_column(f["obs"][col])
            else:
                out[col] = np.array([""] * f["obs"][list(f["obs"].keys())[0]].shape[0])
    return pd.DataFrame(out)


def load_expression(npz_path: Path) -> np.ndarray:
    """Load the (n_cells, 4) ADORA expression matrix saved by the extract scripts."""
    d = np.load(npz_path, allow_pickle=True)
    expr = d["expression"]
    saved_genes = list(d["genes"]) if "genes" in d.files else GENES
    # Reorder to the canonical GENES order
    order = [saved_genes.index(g) for g in GENES]
    return expr[:, order]


def pseudobulk(expr: np.ndarray, groups: np.ndarray, source: str) -> pd.DataFrame:
    """Mean expression and percent-expressing per group, per gene."""
    df = pd.DataFrame(expr, columns=GENES)
    df["group"] = groups
    grouped = df.groupby("group", observed=True)
    rows = []
    for gname, g in grouped:
        if not gname or gname == "nan":
            continue
        n = len(g)
        row = {"source": source, "cell_type": gname, "n_cells": n}
        for gene in GENES:
            vals = g[gene].values
            row[f"{gene}_mean"] = float(vals.mean())
            row[f"{gene}_pct_expressing"] = float((vals > 0).mean() * 100)
            nonzero = vals[vals > 0]
            row[f"{gene}_mean_nonzero"] = float(nonzero.mean()) if len(nonzero) else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def ranked_top20_plot(pb: pd.DataFrame, out_path: Path, min_cells: int = 50) -> None:
    """Per receptor, plot the top 20 cell types by mean_nonzero expression
    (with a min-cells gate to suppress tiny groups). Annotate with source."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for ax, gene in zip(axes.flat, GENES):
        sub = pb[pb["n_cells"] >= min_cells].copy()
        sub = sub.sort_values(f"{gene}_mean_nonzero", ascending=False).head(20)
        if sub.empty:
            ax.set_visible(False)
            continue
        colors = ["#1f77b4" if s == "Tabula Sapiens" else "#d62728" for s in sub["source"]]
        labels = [
            f"{ct}  ({s[0]}, n={int(n)})"
            for ct, s, n in zip(sub["cell_type"], sub["source"], sub["n_cells"])
        ]
        ypos = np.arange(len(sub))
        ax.barh(ypos, sub[f"{gene}_mean_nonzero"], color=colors)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(f"{gene} mean expression in expressing cells")
        ax.set_title(f"{gene} — top 20 cell types ({sub['source'].nunique()} atlases)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color="#1f77b4", label="Tabula Sapiens"),
        plt.Rectangle((0, 0), 1, 1, color="#d62728", label="HBCA non-neuronal"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Top 20 cell types per ADORA receptor (mean expression among expressing cells)",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def cross_receptor_overlap(expr: np.ndarray, cell_types: np.ndarray, source: str,
                           threshold: float = 0.0) -> pd.DataFrame:
    """For each cell type, how many cells express each combination of receptors?

    A cell counts as expressing a receptor when expression > threshold.
    Rows are cell types; columns are combination labels.
    """
    positive = expr > threshold
    rows = []
    for ct in np.unique(cell_types):
        if not ct or ct == "nan":
            continue
        mask = cell_types == ct
        pos = positive[mask]
        n = int(mask.sum())
        row = {"source": source, "cell_type": ct, "n_cells": n}
        # Per-gene positives
        for j, g in enumerate(GENES):
            row[f"{g}_pos"] = int(pos[:, j].sum())
        # Multi-receptor combinations (subsets of size >= 2)
        for k in range(2, len(GENES) + 1):
            for combo in combinations(range(len(GENES)), k):
                label = "+".join(GENES[i] for i in combo)
                m = np.ones(pos.shape[0], dtype=bool)
                for i in combo:
                    m &= pos[:, i]
                row[label] = int(m.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def stratified_dotplot(expr: np.ndarray, strat: np.ndarray, ct: np.ndarray,
                       title: str, out_path: Path, top_n: int = 20) -> None:
    """Dotplot: rows = top cell types by total ADORA detection, cols = (gene, stratum).

    Dot size = pct expressing, dot color = mean expression (in stratum, in cell type).
    """
    df = pd.DataFrame(expr, columns=GENES)
    df["cell_type"] = ct
    df["stratum"] = strat
    df = df[(df["cell_type"] != "") & (df["stratum"] != "")]

    # Pick top N cell types by total positive ADORA cells (any receptor)
    df["any_pos"] = (df[GENES] > 0).any(axis=1)
    top_cells = df.groupby("cell_type")["any_pos"].sum().sort_values(ascending=False).head(top_n).index
    df = df[df["cell_type"].isin(top_cells)]

    strata = sorted(df["stratum"].unique())
    if len(strata) > 24:
        # Keep the strata with the most cells
        strata = (
            df["stratum"].value_counts().head(24).index.tolist()
        )
        df = df[df["stratum"].isin(strata)]
        strata = sorted(strata)

    grouped = df.groupby(["cell_type", "stratum"], observed=True)
    means = grouped[GENES].mean()
    pcts = grouped[GENES].apply(lambda g: (g > 0).mean() * 100)

    # Build a long-form dataframe with one row per (cell_type, stratum, gene)
    rec = []
    for (ct_, st_), m in means.iterrows():
        p = pcts.loc[(ct_, st_)]
        for gene in GENES:
            rec.append({"cell_type": ct_, "stratum": st_, "gene": gene,
                        "mean": float(m[gene]), "pct": float(p[gene])})
    long = pd.DataFrame(rec)

    cell_types_sorted = list(top_cells)
    col_keys = [(g, s) for g in GENES for s in strata]

    n_rows = len(cell_types_sorted)
    n_cols = len(col_keys)
    fig_w = max(10, 0.25 * n_cols)
    fig_h = max(6, 0.32 * n_rows)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    x_index = {(g, s): i for i, (g, s) in enumerate(col_keys)}
    y_index = {ct_: i for i, ct_ in enumerate(cell_types_sorted)}

    xs, ys, sizes, colors = [], [], [], []
    for _, r in long.iterrows():
        if r["cell_type"] not in y_index:
            continue
        if (r["gene"], r["stratum"]) not in x_index:
            continue
        xs.append(x_index[(r["gene"], r["stratum"])])
        ys.append(y_index[r["cell_type"]])
        sizes.append(max(r["pct"] * 4, 2))
        colors.append(r["mean"])

    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis",
                    vmin=0, vmax=max(0.5, np.percentile(colors, 95) if colors else 1.0),
                    edgecolor="none")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([f"{g}\n{s}" for g, s in col_keys], rotation=90, fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(cell_types_sorted, fontsize=8)
    # Vertical separators between gene blocks
    for i in range(1, len(GENES)):
        ax.axvline(i * len(strata) - 0.5, color="grey", linewidth=0.5)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.invert_yaxis()
    ax.set_title(title)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, label="Mean expression")
    # Size legend
    for pct in (1, 5, 15, 30):
        ax.scatter([], [], s=pct * 4, color="lightgrey", label=f"{pct}%")
    ax.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.18, 1.0),
              frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for p in (TABULA_H5AD, TABULA_EXPR, HBCA_H5AD, HBCA_EXPR):
        if not p.exists():
            sys.exit(f"missing input: {p}")
    FIG_DIR.mkdir(exist_ok=True)
    CACHE_DIR.mkdir(exist_ok=True)

    print("Loading Tabula Sapiens metadata...")
    tabula_obs = load_obs(TABULA_H5AD, ["cell_type", "tissue", "donor_id", "assay"])
    tabula_expr = load_expression(TABULA_EXPR)
    if tabula_expr.shape[0] != len(tabula_obs):
        sys.exit(f"Tabula expr/obs mismatch: {tabula_expr.shape[0]} vs {len(tabula_obs)}")

    print("Loading HBCA metadata...")
    hbca_obs = load_obs(HBCA_H5AD, ["cell_type", "tissue", "supercluster_term", "donor_id", "assay"])
    hbca_expr = load_expression(HBCA_EXPR)
    if hbca_expr.shape[0] != len(hbca_obs):
        sys.exit(f"HBCA expr/obs mismatch: {hbca_expr.shape[0]} vs {len(hbca_obs)}")

    print(f"Tabula: {len(tabula_obs):,} cells; HBCA: {len(hbca_obs):,} cells")

    # Pseudobulk per source, then concat
    print("Building pseudobulk_by_cell_type...")
    pb_tabula = pseudobulk(tabula_expr, tabula_obs["cell_type"].values, "Tabula Sapiens")
    pb_hbca = pseudobulk(hbca_expr, hbca_obs["cell_type"].values, "HBCA non-neuronal")
    pb = pd.concat([pb_tabula, pb_hbca], ignore_index=True)
    pb_out = CACHE_DIR / "pseudobulk_by_cell_type.feather"
    pb.to_feather(pb_out)
    print(f"  → {pb_out}  ({len(pb)} rows)")

    # Ranked top-20 figure
    print("Building ranked_top20_per_receptor.png...")
    fig_out = FIG_DIR / "ranked_top20_per_receptor.png"
    ranked_top20_plot(pb, fig_out)
    print(f"  → {fig_out}")

    # Cross-receptor overlap
    print("Building cross_receptor_overlap...")
    ovl_tabula = cross_receptor_overlap(tabula_expr, tabula_obs["cell_type"].values, "Tabula Sapiens")
    ovl_hbca = cross_receptor_overlap(hbca_expr, hbca_obs["cell_type"].values, "HBCA non-neuronal")
    ovl = pd.concat([ovl_tabula, ovl_hbca], ignore_index=True)
    ovl_out = CACHE_DIR / "cross_receptor_overlap.feather"
    ovl.to_feather(ovl_out)
    print(f"  → {ovl_out}  ({len(ovl)} rows)")

    # Donor- and assay-stratified dotplots (Tabula only; HBCA donor/assay heterogeneity differs)
    print("Building donor_stratified_dotplot.png...")
    stratified_dotplot(
        tabula_expr, tabula_obs["donor_id"].values, tabula_obs["cell_type"].values,
        "Tabula Sapiens — ADORA × top cell types × donor",
        FIG_DIR / "donor_stratified_dotplot.png",
    )
    print(f"  → {FIG_DIR / 'donor_stratified_dotplot.png'}")

    print("Building assay_stratified_dotplot.png...")
    stratified_dotplot(
        tabula_expr, tabula_obs["assay"].values, tabula_obs["cell_type"].values,
        "Tabula Sapiens — ADORA × top cell types × assay",
        FIG_DIR / "assay_stratified_dotplot.png",
    )
    print(f"  → {FIG_DIR / 'assay_stratified_dotplot.png'}")

    summary = {
        "tabula_cells": int(len(tabula_obs)),
        "hbca_cells": int(len(hbca_obs)),
        "genes": GENES,
        "pseudobulk_rows": int(len(pb)),
        "cross_receptor_rows": int(len(ovl)),
        "missing_coverage": "HBCA neurons (2.48M cells, 30 GB) not yet downloaded; "
                            "brain neuronal ADORA1/2A signal absent.",
    }
    (CACHE_DIR / "q1_summary.json").write_text(json.dumps(summary, indent=2))
    print("Done.")


if __name__ == "__main__":
    main()
