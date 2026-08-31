#!/usr/bin/env python
"""Co-expression check: exhaustive A1-pathway gene panel in Hippocampal CA4.

Groundwork for blogs/a1-receptor-pathway. Pulls a curated gene panel
spanning every documented ADORA1 signaling branch (Gi/cAMP/PKA, Gbg->GIRK,
Gbg->presynaptic Ca2+ channels, PLC/PKC/Ca2+, PI3K/MAPK, NO/cGMP) and
checks pct_expressing / mean_nonzero / mean_all for each gene, restricted
to cells labeled "Hippocampal CA4" in obs/supercluster_term of the HBCA
neuron atlas (2,480,956 cells total; CA4 is our #1 ADORA1 hit, n=10,654).

This answers "does CA4 have the machinery," not "is the pathway active in
CA4" — see the caveats in distribution-scratchpad.md and the earlier
RNA-seq interpretation discussion (mRNA presence != protein activity).
Also no baseline comparison cell type is included here yet, which matters:
most of these genes are broadly-expressed general signaling machinery, so
high presence in CA4 alone doesn't establish CA4-specificity.

Output: figures/a1_pathway_hippocampal_ca4_coexpression.csv
"""

from __future__ import annotations
from pathlib import Path
import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd

LAB = Path("/home/acefsan/src/dl_bio/projects/caffeine/lab/001_adora_expression")
H5AD = LAB / "cache/human_brain_cell_atlas/hbca_all_neurons_8e10f1c4.h5ad"

BRANCHES = {
    "Receptor": ["ADORA1"],
    "Gi/o alpha subunits": ["GNAI1", "GNAI2", "GNAI3", "GNAO1"],
    "Adenylyl cyclase (Gi-inhibited isoforms)": ["ADCY1", "ADCY5", "ADCY6", "ADCY8"],
    "PKA subunits": ["PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKACA", "PRKACB"],
    "cAMP downstream (CREB)": ["CREB1"],
    "Gbg -> GIRK channels": ["KCNJ3", "KCNJ6", "KCNJ9", "KCNJ5"],
    "Gbg -> presynaptic Ca2+ channels": ["CACNA1B", "CACNA1A"],
    "PLC/PKC/Ca2+ branch": ["PLCB1", "PLCB2", "PLCB3", "PLCB4", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "CALM1", "CALM2", "CALM3"],
    "NF-kB (downstream of PLC/PKC branch)": ["NFKB1", "RELA"],
    "PI3K/MAPK branch": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PREX1", "RAC1",
                          "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
                          "MAPK8", "MAPK9", "MAPK10", "MAPK14", "MAPK11",
                          "AKT1", "AKT2", "AKT3", "VEGFA"],
    "NO/cGMP branch": ["NOS1", "NOS2", "NOS3", "GUCY1A1", "GUCY1B1", "PRKG1", "PRKG2"],
}
ALL_GENES = [g for genes in BRANCHES.values() for g in genes]
gene_to_branch = {g: b for b, genes in BRANCHES.items() for g in genes}

def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

with h5py.File(H5AD, "r") as f:
    fn = f["var/feature_name"]
    symbols = np.array(dec(fn["categories"][:]))[fn["codes"][:]]
    sym_index = {s: i for i, s in enumerate(symbols)}

    found = [g for g in ALL_GENES if g in sym_index]
    missing = [g for g in ALL_GENES if g not in sym_index]
    cols = [sym_index[g] for g in found]
    print(f"Found {len(found)}/{len(ALL_GENES)} genes. Missing: {missing}")

    sc = f["obs/supercluster_term"]
    cats = np.array(dec(sc["categories"][:]))
    codes = sc["codes"][:]
    labels = np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
    mask = labels == "Hippocampal CA4"
    n_ca4 = int(mask.sum())
    print(f"Hippocampal CA4 cells: {n_ca4}")

    X = f["X"]
    indptr = X["indptr"][:]
    n = indptr.shape[0] - 1
    ncols = int(X.attrs["shape"][1])
    data = X["data"]
    ind = X["indices"]

    collected = []
    B = 200_000
    for s in range(0, n, B):
        e = min(s + B, n)
        blk_mask = mask[s:e]
        if not blk_mask.any():
            continue
        d0, d1 = int(indptr[s]), int(indptr[e])
        blk = sp.csr_matrix(
            (data[d0:d1], ind[d0:d1], indptr[s:e + 1] - d0),
            shape=(e - s, ncols),
        )
        sub = blk[blk_mask][:, cols].toarray()
        collected.append(sub)
        print(f"  rows {s:,}-{e:,}: kept {blk_mask.sum()}", flush=True)

    arr = np.concatenate(collected, axis=0)
    assert arr.shape[0] == n_ca4

df = pd.DataFrame(arr, columns=found)
rows = []
for g in found:
    vals = df[g].values
    nonzero = vals[vals > 0]
    rows.append({
        "branch": gene_to_branch[g],
        "gene": g,
        "n_ca4_cells": n_ca4,
        "pct_expressing": float((vals > 0).mean() * 100),
        "mean_nonzero": float(nonzero.mean()) if len(nonzero) else 0.0,
        "mean_all": float(vals.mean()),
    })
out = pd.DataFrame(rows).sort_values(["branch", "pct_expressing"], ascending=[True, False])
out_path = LAB / "figures/a1_pathway_hippocampal_ca4_coexpression.csv"
out.to_csv(out_path, index=False)
print(f"\nSaved {out_path}")
print(out.to_string(index=False))
if missing:
    print("\nMissing gene symbols (not found in var/feature_name):", missing)
