#!/usr/bin/env python
"""Extract ADORA expression for all HBCA neurons.

The neuron H5AD (hbca_all_neurons_8e10f1c4.h5ad, 2,480,956 cells) has an X CSR
matrix (no raw) and var/feature_name as an AnnData categorical. We pull the four
ADORA genes for every cell in memory-safe CSR blocks and save a compact npz that
mirrors the non-neuronal npz layout (expression shape (n_cells, 4), genes).
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

GENES = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
LAB_DIR = Path(__file__).resolve().parent
CACHE_DIR = LAB_DIR / "cache" / "human_brain_cell_atlas"
H5AD = CACHE_DIR / "hbca_all_neurons_8e10f1c4.h5ad"
OUT = CACHE_DIR / "hbca_neurons_adora_expression.npz"
EXPECTED_ROWS = 2_480_956


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def main() -> None:
    with h5py.File(H5AD, "r") as f:
        fn = f["var/feature_name"]
        symbols = np.array(dec(fn["categories"][:]))[fn["codes"][:]]  # aligned to columns
        cols = [int(np.where(symbols == g)[0][0]) for g in GENES]
        print(f"gene -> column index: {dict(zip(GENES, cols))}")

        X = f["X"]
        indptr = X["indptr"][:]
        n = indptr.shape[0] - 1
        ncols = int(X.attrs["shape"][1])
        data = X["data"]
        ind = X["indices"]
        print(f"X shape: ({n}, {ncols})")

        out = np.zeros((n, 4), dtype=np.float32)
        B = 100_000
        for s in range(0, n, B):
            e = min(s + B, n)
            d0, d1 = int(indptr[s]), int(indptr[e])
            blk = sp.csr_matrix(
                (data[d0:d1], ind[d0:d1], indptr[s : e + 1] - d0),
                shape=(e - s, ncols),
            )
            out[s:e] = blk[:, cols].toarray()
            print(f"  rows {s:,}-{e:,}", flush=True)

    assert out.shape[0] == EXPECTED_ROWS, f"row mismatch: {out.shape[0]} != {EXPECTED_ROWS}"
    np.savez_compressed(OUT, expression=out, genes=np.array(GENES, dtype=object))
    print(f"Saved {OUT} with shape {out.shape}")


if __name__ == "__main__":
    main()
