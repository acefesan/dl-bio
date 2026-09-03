#!/usr/bin/env python
"""Per-BRANCH module scoring of the A1 pathway panel, by cell type.

Same AddModuleScore-style method as score_a1_module_pathway.py (Tirosh et
al. 2016 control-gene binning, run at pseudobulk/cell-type level), but
scored separately for each of the six documented ADORA1 signaling
branches instead of collapsing all 61 genes into one number. Answers
"which cell types express WHICH branch," which the combined score can't.

Control genes are sampled per panel gene and tracked back to their gene,
so each branch is scored against the controls matched to that branch's
own genes (not a single global control pool).

Also dumps the full per-gene x per-cell-type mean table, which is what's
needed to explain outliers (e.g. why Thalamic excitatory scores ~0 on the
combined panel despite ranking #2 for ADORA1 alone).

Outputs:
  figures/a1_branch_scores_<atlas>.csv       (group x branch scores)
  figures/a1_pergene_means_<atlas>.csv       (group x gene mean expression)
"""
from __future__ import annotations
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

LAB = Path("/home/acefsan/src/dl_bio/projects/caffeine/lab/001_adora_expression")

BRANCHES = {
    "receptor": ["ADORA1"],
    "Gi_alpha": ["GNAI1", "GNAI2", "GNAI3", "GNAO1"],
    "adenylyl_cyclase": ["ADCY1", "ADCY5", "ADCY6", "ADCY8"],
    "PKA": ["PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKACA", "PRKACB"],
    "cAMP_downstream_CREB": ["CREB1"],
    "Gbg_GIRK": ["KCNJ3", "KCNJ6", "KCNJ9", "KCNJ5"],
    "Gbg_Ca_channels": ["CACNA1B", "CACNA1A"],
    "PLC_PKC_Ca": ["PLCB1", "PLCB2", "PLCB3", "PLCB4", "PRKCA", "PRKCB",
                    "PRKCD", "PRKCE", "CALM1", "CALM2", "CALM3"],
    "NFkB": ["NFKB1", "RELA"],
    "PI3K_MAPK": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PREX1", "RAC1",
                   "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK8", "MAPK9",
                   "MAPK10", "MAPK11", "MAPK14", "AKT1", "AKT2", "AKT3", "VEGFA"],
    "NO_cGMP": ["NOS1", "NOS2", "NOS3", "GUCY1A1", "GUCY1B1", "PRKG1", "PRKG2"],
}
PANEL = [g for gs in BRANCHES.values() for g in gs]

NBINS = 25
CTRL_PER_GENE = 20
RNG = np.random.default_rng(0)


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def run_atlas(h5ad_path: Path, group_col: str, label: str, slug: str):
    print(f"\n=== {label} ({group_col}) ===", flush=True)
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
        data, ind = X["data"], X["indices"]

        # Pass 1: global per-gene mean, for expression-matched binning.
        gene_sum = np.zeros(ncols, dtype=np.float64)
        B = 200_000
        for s in range(0, n, B):
            e = min(s + B, n)
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[d0:d1], ind[d0:d1], indptr[s:e+1] - d0), shape=(e-s, ncols))
            gene_sum += np.asarray(blk.sum(axis=0)).ravel()
            print(f"  pass1 {s:,}", flush=True)
        gene_mean = gene_sum / n

        order = np.argsort(gene_mean)
        col_bin = np.empty(ncols, dtype=int)
        col_bin[order] = np.floor(np.linspace(0, NBINS, ncols, endpoint=False)).astype(int)

        panel_found = [g for g in PANEL if g in sym_index]
        panel_cols = {g: sym_index[g] for g in panel_found}
        panel_col_set = set(panel_cols.values())

        bins_to_cols = {}
        for c in range(ncols):
            bins_to_cols.setdefault(col_bin[c], []).append(c)

        # Controls tracked per panel gene, so branches get matched controls.
        ctrl_for_gene = {}
        for g, c in panel_cols.items():
            cand = [x for x in bins_to_cols[col_bin[c]] if x not in panel_col_set]
            k = min(CTRL_PER_GENE, len(cand))
            ctrl_for_gene[g] = RNG.choice(cand, size=k, replace=False).tolist() if k else []

        all_ctrl = sorted({c for v in ctrl_for_gene.values() for c in v})
        want = sorted(set(panel_cols.values()) | set(all_ctrl))
        widx = {c: i for i, c in enumerate(want)}
        print(f"  panel {len(panel_found)}/{len(PANEL)}, controls {len(all_ctrl)}", flush=True)

        groups = [g for g in pd.unique(group_labels) if g and g != "nan"]
        sums = {g: np.zeros(len(want)) for g in groups}
        counts = {g: 0 for g in groups}
        for s in range(0, n, B):
            e = min(s + B, n)
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix((data[d0:d1], ind[d0:d1], indptr[s:e+1] - d0), shape=(e-s, ncols))
            sub = blk[:, want].toarray()
            df = pd.DataFrame(sub)
            df["g"] = group_labels[s:e]
            gs_ = df.groupby("g", observed=True).sum()
            gc_ = df.groupby("g", observed=True).size()
            for g in gs_.index:
                if g in sums:
                    sums[g] += gs_.loc[g].values
                    counts[g] += int(gc_.loc[g])
            print(f"  pass2 {s:,}", flush=True)

    means = {g: sums[g] / counts[g] for g in groups if counts[g] > 0}

    # Per-gene mean table
    pergene = pd.DataFrame(
        {g: [means[g][widx[panel_cols[gene]]] for gene in panel_found] for g in means},
        index=panel_found,
    ).T
    pergene.insert(0, "n_cells", [counts[g] for g in means])
    pergene.index.name = "group"
    pergene.to_csv(LAB / f"figures/a1_pergene_means_{slug}.csv")

    # Per-branch scores
    rows = []
    for g in means:
        row = {"atlas": label, "group": g, "n_cells": counts[g]}
        for bname, bgenes in BRANCHES.items():
            bg = [x for x in bgenes if x in panel_cols]
            if not bg:
                row[bname] = np.nan
                continue
            p = np.mean([means[g][widx[panel_cols[x]]] for x in bg])
            cc = sorted({c for x in bg for c in ctrl_for_gene[x]})
            c_ = np.mean([means[g][widx[c]] for c in cc]) if cc else 0.0
            row[bname] = p - c_
        rows.append(row)
    out = pd.DataFrame(rows)
    out["combined"] = out[list(BRANCHES)].mean(axis=1)
    out = out.sort_values("combined", ascending=False)
    out.to_csv(LAB / f"figures/a1_branch_scores_{slug}.csv", index=False)
    print(out.to_string(index=False), flush=True)


if __name__ == "__main__":
    run_atlas(LAB / "cache/human_brain_cell_atlas/hbca_all_neurons_8e10f1c4.h5ad",
              "supercluster_term", "HBCA neurons", "hbca_neurons")
    run_atlas(LAB / "cache/tabula_sapiens_all_cells.h5ad",
              "cell_type", "Tabula Sapiens", "tabula_sapiens")
