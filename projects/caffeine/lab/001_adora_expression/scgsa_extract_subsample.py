"""Subsample cells per cell-type group and export raw counts in 10x HDF5 form.

Why subsample: the rank-based scGSA methods (AUCell, UCell, JASMINE, ssGSEA)
score each cell against the full gene axis, which is not tractable across
millions of cells. Wang & Thakar benchmark at 20/50/200/500 cells per group
and use 200 per group for their gene-set-size and noise experiments, where
all methods had stabilised. We take 200 per group with the same rationale.

Why raw counts from both atlases: HBCA stores raw counts in X, Tabula stores
log-normalised values in X and raw counts in raw/X. Exporting counts from both
and letting Seurat::NormalizeData run identically removes the cross-atlas
scale mismatch that made earlier magnitudes non-comparable.

Output per atlas (10x-style HDF5 that Seurat::Read10X_h5 reads directly):
  <out>/<atlas>_counts.h5   genes x cells CSC
  <out>/<atlas>_cells.csv   barcode -> group label
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

LAB = Path("/home/acefsan/src/dl_bio/projects/caffeine/lab/001_adora_expression")
OUT = Path(os.environ.get("SCGSA_WORK", "scgsa_work")) / "data"
OUT.mkdir(parents=True, exist_ok=True)

N_PER_GROUP = 200
MIN_GROUP = 500
SEED = 0
BLOCK = 200_000

ATLASES = {
    "hbca": dict(
        path=LAB / "cache/human_brain_cell_atlas/hbca_all_neurons_8e10f1c4.h5ad",
        group_col="supercluster_term",
        x_root="X",
        var_root="var",
    ),
    "tabula": dict(
        path=LAB / "cache/tabula_sapiens_all_cells.h5ad",
        group_col="cell_type",
        x_root="raw/X",
        var_root="raw/var",
    ),
}


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def cat_column(f, path):
    g = f[path]
    cats = np.array(dec(g["categories"][:]))
    codes = g["codes"][:]
    return cats, codes


def run(name: str, cfg: dict) -> None:
    print(f"\n=== {name} ===", flush=True)
    with h5py.File(cfg["path"], "r") as f:
        var_root = cfg["var_root"]
        if f"{var_root}/feature_name" in f:
            cats, codes = cat_column(f, f"{var_root}/feature_name")
            symbols = cats[codes]
        else:  # raw/var may lack feature_name; fall back to main var
            cats, codes = cat_column(f, "var/feature_name")
            symbols = cats[codes]

        gcats, gcodes = cat_column(f, f"obs/{cfg['group_col']}")
        n_cells_total = gcodes.shape[0]

        counts = pd.Series(gcodes).value_counts()
        keep_codes = sorted(c for c, n in counts.items() if n >= MIN_GROUP)
        print(f"  groups kept (>= {MIN_GROUP} cells): "
              f"{len(keep_codes)} / {len(gcats)}")

        rng = np.random.default_rng(SEED)
        chosen = []
        for c in keep_codes:
            idx = np.flatnonzero(gcodes == c)
            take = min(N_PER_GROUP, idx.shape[0])
            chosen.append(rng.choice(idx, size=take, replace=False))
        chosen = np.sort(np.concatenate(chosen))
        print(f"  sampled {chosen.shape[0]:,} cells of {n_cells_total:,}")

        X = f[cfg["x_root"]]
        indptr = X["indptr"][:]
        data_ds, ind_ds = X["data"], X["indices"]
        n_genes = int(X.attrs["shape"][1]) if "shape" in X.attrs else len(symbols)

        chosen_set = 0
        blocks = []
        for s in range(0, n_cells_total, BLOCK):
            e = min(s + BLOCK, n_cells_total)
            sel = chosen[(chosen >= s) & (chosen < e)]
            if sel.size == 0:
                continue
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix(
                (data_ds[d0:d1], ind_ds[d0:d1], indptr[s:e + 1] - d0),
                shape=(e - s, n_genes),
            )
            blocks.append(blk[sel - s])
            chosen_set += sel.size
            print(f"    rows {s:,}-{e:,}  (+{sel.size})", flush=True)

        mat = sp.vstack(blocks, format="csr")
        assert mat.shape[0] == chosen.shape[0] == chosen_set

        groups = gcats[gcodes[chosen]]

    # CSR over (cells x genes) has identical arrays to CSC over (genes x cells),
    # which is exactly the layout Read10X_h5 expects.
    mat = mat.astype(np.float32)
    barcodes = np.array([f"cell_{i}" for i in range(mat.shape[0])])

    out_h5 = OUT / f"{name}_counts.h5"
    if out_h5.exists():
        out_h5.unlink()
    with h5py.File(out_h5, "w") as o:
        m = o.create_group("matrix")
        m.create_dataset("data", data=mat.data, compression="gzip",
                         compression_opts=1)
        m.create_dataset("indices", data=mat.indices.astype(np.int64),
                         compression="gzip", compression_opts=1)
        m.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array([n_genes, mat.shape[0]],
                                                dtype=np.int32))
        m.create_dataset("barcodes", data=barcodes.astype("S"))
        feat = m.create_group("features")
        feat.create_dataset("id", data=symbols.astype("S"))
        feat.create_dataset("name", data=symbols.astype("S"))
        feat.create_dataset("feature_type",
                            data=np.array(["Gene Expression"] * n_genes,
                                          dtype="S"))
        feat.create_dataset("genome",
                            data=np.array(["GRCh38"] * n_genes, dtype="S"))

    pd.DataFrame({"barcode": barcodes, "group": groups}).to_csv(
        OUT / f"{name}_cells.csv", index=False)

    nnz = int(mat.nnz)
    print(f"  wrote {out_h5}  ({n_genes:,} genes x {mat.shape[0]:,} cells, "
          f"nnz={nnz:,})")
    print(f"  integral counts: {bool(np.allclose(mat.data[:100000], np.round(mat.data[:100000])))}")


if __name__ == "__main__":
    which = sys.argv[1:] or list(ATLASES)
    for nm in which:
        run(nm, ATLASES[nm])
