"""Simplified module-score (Tirosh/Seurat AddModuleScore-style) pathway
scoring, run at pseudobulk (per-cell-type) level rather than per-cell, for
compute efficiency. Reused across both atlases via ATLAS_CONFIG below.

Method:
1. Compute each gene's global mean expression across ALL cells (one pass).
2. Bin all genes into N expression bins by that global mean (quantile bins).
3. For each pathway gene, sample K control genes from its same bin
   (excluding pathway genes) -> control set is expression-matched, so the
   comparison isn't just "pathway genes vs. average gene" (biased by how
   highly the pathway genes are expressed in general).
4. Per cell-type group: score = mean(pathway gene group-means) -
   mean(control gene group-means).

Simplification vs. the textbook per-cell AddModuleScore: this scores each
CELL-TYPE GROUP once (using group-level pseudobulk means), not each cell
individually. Cheaper, and sufficient to answer "which cell types are
elevated for this pathway relative to background" -- the question this
was built for -- but not a substitute for true per-cell scoring if later
work needs cell-level (not group-level) values.
"""
from __future__ import annotations
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

LAB = Path("/home/acefsan/src/dl_bio/projects/caffeine/lab/001_adora_expression")

PANEL = [
    "ADORA1", "GNAI1", "GNAI2", "GNAI3", "GNAO1",
    "ADCY1", "ADCY5", "ADCY6", "ADCY8",
    "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKACA", "PRKACB",
    "CREB1",
    "KCNJ3", "KCNJ6", "KCNJ9", "KCNJ5",
    "CACNA1B", "CACNA1A",
    "PLCB1", "PLCB2", "PLCB3", "PLCB4", "PRKCA", "PRKCB", "PRKCD", "PRKCE",
    "CALM1", "CALM2", "CALM3",
    "NFKB1", "RELA",
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PREX1", "RAC1",
    "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
    "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK14",
    "AKT1", "AKT2", "AKT3", "VEGFA",
    "NOS1", "NOS2", "NOS3", "GUCY1A1", "GUCY1B1", "PRKG1", "PRKG2",
]

NBINS = 25
CTRL_PER_GENE = 20
RNG = np.random.default_rng(0)


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def run_atlas(h5ad_path: Path, group_col: str, label: str, out_csv: Path):
    print(f"\n=== {label} ({h5ad_path.name}, grouped by {group_col}) ===")
    with h5py.File(h5ad_path, "r") as f:
        fn = f["var/feature_name"]
        symbols = np.array(dec(fn["categories"][:]))[fn["codes"][:]]
        sym_index = {s: i for i, s in enumerate(symbols)}

        gc = f[f"obs/{group_col}"]
        cats = np.array(dec(gc["categories"][:]))
        codes = gc["codes"][:]
        group_labels = np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)

        X = f["X"]
        indptr = X["indptr"][:]
        n = indptr.shape[0] - 1
        ncols = int(X.attrs["shape"][1])
        data = X["data"]
        ind = X["indices"]

        # Pass 1: global per-gene mean expression (all genes), one sparse pass.
        gene_sum = np.zeros(ncols, dtype=np.float64)
        B = 200_000
        for s in range(0, n, B):
            e = min(s + B, n)
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[d0:d1], ind[d0:d1], indptr[s:e + 1] - d0), shape=(e - s, ncols))
            gene_sum += np.asarray(blk.sum(axis=0)).ravel()
            print(f"  pass1 rows {s:,}-{e:,}", flush=True)
        gene_mean = gene_sum / n

        # Bin all genes by global mean expression (quantile bins).
        order = np.argsort(gene_mean)
        bin_id = np.empty(ncols, dtype=int)
        bin_id[order] = np.floor(np.linspace(0, NBINS, ncols, endpoint=False)).astype(int)
        gene_bin = pd.Series(bin_id, index=symbols[np.arange(ncols)]) if False else None
        # symbols may repeat (rare) -- build bin lookup by column index instead
        col_bin = bin_id  # indexed by column position

        panel_cols = [sym_index[g] for g in PANEL if g in sym_index]
        panel_found = [g for g in PANEL if g in sym_index]
        missing = [g for g in PANEL if g not in sym_index]
        panel_col_set = set(panel_cols)

        # Build control set: for each panel gene, sample CTRL_PER_GENE genes
        # from the same expression bin, excluding panel genes themselves.
        control_cols = []
        bins_to_cols = {}
        for col in range(ncols):
            bins_to_cols.setdefault(col_bin[col], []).append(col)
        for col in panel_cols:
            candidates = [c for c in bins_to_cols[col_bin[col]] if c not in panel_col_set]
            if len(candidates) == 0:
                continue
            k = min(CTRL_PER_GENE, len(candidates))
            control_cols.extend(RNG.choice(candidates, size=k, replace=False).tolist())
        control_cols = sorted(set(control_cols))
        print(f"  panel genes found: {len(panel_found)}/{len(PANEL)} (missing: {missing})")
        print(f"  control genes sampled: {len(control_cols)}")

        want_cols = sorted(set(panel_cols) | set(control_cols))
        want_index = {c: i for i, c in enumerate(want_cols)}

        # Pass 2: per-group pseudobulk mean, restricted to panel+control columns.
        groups = pd.unique(group_labels)
        groups = [g for g in groups if g and g != "nan"]
        sums = {g: np.zeros(len(want_cols)) for g in groups}
        counts = {g: 0 for g in groups}
        for s in range(0, n, B):
            e = min(s + B, n)
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[d0:d1], ind[d0:d1], indptr[s:e + 1] - d0), shape=(e - s, ncols))
            sub = blk[:, want_cols].toarray()
            lbls = group_labels[s:e]
            df = pd.DataFrame(sub)
            df["g"] = lbls
            gsum = df.groupby("g", observed=True).sum()
            gcnt = df.groupby("g", observed=True).size()
            for g in gsum.index:
                if g not in sums:
                    continue
                sums[g] += gsum.loc[g].values
                counts[g] += int(gcnt.loc[g])
            print(f"  pass2 rows {s:,}-{e:,}", flush=True)

    panel_idx = [want_index[c] for c in panel_cols]
    control_idx = [want_index[c] for c in control_cols]

    rows = []
    for g in groups:
        if counts[g] == 0:
            continue
        means = sums[g] / counts[g]
        panel_score = means[panel_idx].mean()
        control_score = means[control_idx].mean()
        rows.append({
            "atlas": label,
            "group": g,
            "n_cells": counts[g],
            "panel_mean": panel_score,
            "control_mean": control_score,
            "module_score": panel_score - control_score,
        })
    out = pd.DataFrame(rows).sort_values("module_score", ascending=False)
    out.to_csv(out_csv, index=False)
    print(f"  saved {out_csv}")
    print(out.head(15).to_string(index=False))
    print("  ...")
    print(out.tail(5).to_string(index=False))
    return out


if __name__ == "__main__":
    run_atlas(
        LAB / "cache/human_brain_cell_atlas/hbca_all_neurons_8e10f1c4.h5ad",
        "supercluster_term",
        "HBCA neurons",
        LAB / "figures/a1_module_score_hbca_neurons.csv",
    )
    run_atlas(
        LAB / "cache/tabula_sapiens_all_cells.h5ad",
        "cell_type",
        "Tabula Sapiens",
        LAB / "figures/a1_module_score_tabula_sapiens.csv",
    )
