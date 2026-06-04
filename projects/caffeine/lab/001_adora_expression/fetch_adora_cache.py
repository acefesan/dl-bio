"""Fetch ADORA cache, v3 — stratified sample via coord-based reads.

Background — read these first:
    wiki/labs/001-fetch-stall-postmortem.md
    wiki/concepts/tiledb-soma-storage.md
    wiki/concepts/census-source-h5ads.md

V2 (per-cell_type chunking with obs_value_filter) confirmed that value-filter
X reads have a ~10h-per-call floor regardless of cell count — every fragment
that might contain matching cells must be walked and decompressed. The 2026-06-01
brain probe took 10h 39m to return 11 cells × 4 genes with nnz=1.

V3 changes the access pattern:

    1. ONE global obs scan to get (soma_joinid, cell_type, tissue_general) for
       all primary normal human cells. Cached as parquet. Slow but one-time;
       the prior brain-only scan took ~6h, so a full-human scan is expected
       in the 1-2 day range.
    2. Stratified sample up to N cells per cell_type in pandas. Local, seconds.
    3. ONE X query with obs_coords=<sampled soma_joinids> and var_value_filter
       for the four ADORA genes. Hypothesis: coordinate lookup on the obs axis
       uses TileDB's dimension index rather than scanning fragments row by row,
       so this call should be dramatically faster than v2's value-filter calls.
    4. Save as cache/adora_stratified.h5ad.

Phase 3 is unproven on this layout. If it is still slow, the post-mortem will
tell us; the cached obs metadata (phase 1) remains reusable for any fallback
strategy (e.g., source H5AD downloads — see census-source-h5ads.md wiki page).

Usage:
    # Full run, default 1000 cells per cell_type
    python fetch_adora_cache.py

    # Smaller sample for sanity testing
    python fetch_adora_cache.py --cells-per-cell-type 100

    # Resume — reuses cached obs scan and final h5ad if present
    python fetch_adora_cache.py --resume

    # Force a fresh obs scan (otherwise the parquet cache is reused)
    python fetch_adora_cache.py --rescan

Output:
    cache/_obs_human_primary_normal.parquet   global obs scan cache
    cache/_sample_metadata.parquet            sampled cell metadata
    cache/adora_stratified.h5ad               final AnnData
    cache/fetch.log                           human log
    cache/stats.jsonl                         per-phase structured metrics
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import cellxgene_census
import pandas as pd
import tiledbsoma

GENES = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
OBS_COLS_SCAN = ["soma_joinid", "cell_type", "tissue_general"]
OBS_COLS_FETCH = [
    "soma_joinid", "cell_type", "tissue", "tissue_general",
    "assay", "dataset_id", "donor_id",
]

SCRIPT_DIR = Path(__file__).parent
CACHE_DIR = SCRIPT_DIR / "cache"
OBS_CACHE = CACHE_DIR / "_obs_human_primary_normal.parquet"
SAMPLE_CACHE = CACHE_DIR / "_sample_metadata.parquet"
FINAL_OUT = CACHE_DIR / "adora_stratified.h5ad"
STATS_FILE = CACHE_DIR / "stats.jsonl"


def setup_logging() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(CACHE_DIR / "fetch.log", mode="a"),
        ],
    )


def emit_stat(record: dict) -> None:
    record = {"ts": time.time(), **record}
    with open(STATS_FILE, "a") as f:
        f.write(json.dumps(record) + "\n")


def open_census():
    """Open Census with generous S3 timeouts and a large initial buffer."""
    ctx = tiledbsoma.SOMATileDBContext(
        tiledb_config={
            "vfs.s3.connect_timeout_ms": "60000",
            "vfs.s3.request_timeout_ms": "300000",
            "soma.init_buffer_bytes": str(512 * 1024 * 1024),
        }
    )
    return cellxgene_census.open_soma(context=ctx)


def scan_human_obs(census, *, use_cache: bool) -> pd.DataFrame:
    """Phase 1: global obs scan. Returns DataFrame[soma_joinid, cell_type, tissue_general]."""
    if use_cache and OBS_CACHE.exists():
        df = pd.read_parquet(OBS_CACHE)
        logging.info(
            "[scan] cached: %d cells across %d cell_types, %d tissues (from %s)",
            len(df), df["cell_type"].nunique(),
            df["tissue_general"].nunique(), OBS_CACHE.name,
        )
        return df

    logging.info(
        "[scan] global obs.read for primary normal human — expect many hours"
    )
    t0 = time.monotonic()
    obs = (
        census["census_data"]["homo_sapiens"]
        .obs.read(
            column_names=OBS_COLS_SCAN,
            value_filter="is_primary_data == True and disease == 'normal'",
        )
        .concat()
        .to_pandas()
    )
    dt = time.monotonic() - t0

    logging.info(
        "[scan] done: %d cells, %d cell_types, %d tissues, %.0fs",
        len(obs), obs["cell_type"].nunique(),
        obs["tissue_general"].nunique(), dt,
    )

    OBS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    obs.to_parquet(OBS_CACHE)
    logging.info("[scan] wrote %s (%.1f MB)",
                 OBS_CACHE.name, OBS_CACHE.stat().st_size / 1e6)

    emit_stat({
        "op": "scan", "cells": int(len(obs)),
        "cell_types": int(obs["cell_type"].nunique()),
        "tissues": int(obs["tissue_general"].nunique()),
        "seconds": dt,
    })
    return obs


def stratified_sample(obs: pd.DataFrame, n_per_ct: int, seed: int) -> pd.DataFrame:
    """Phase 2: sample up to n_per_ct cells per cell_type. Local, seconds."""
    t0 = time.monotonic()
    sampled = (
        obs.groupby("cell_type", group_keys=False, observed=True)
        .apply(lambda g: g.sample(n=min(len(g), n_per_ct), random_state=seed))
        .reset_index(drop=True)
    )
    dt = time.monotonic() - t0
    logging.info(
        "[sample] %d cells across %d cell_types (cap=%d/ct) in %.1fs",
        len(sampled), sampled["cell_type"].nunique(), n_per_ct, dt,
    )
    SAMPLE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    sampled.to_parquet(SAMPLE_CACHE)
    logging.info("[sample] wrote %s", SAMPLE_CACHE.name)
    emit_stat({
        "op": "sample", "cells": int(len(sampled)),
        "cell_types": int(sampled["cell_type"].nunique()),
        "cap_per_cell_type": n_per_ct, "seed": seed, "seconds": dt,
    })
    return sampled


def fetch_by_coords(census, soma_joinids: list[int]) -> ad.AnnData:
    """Phase 3: single coord-based X read. Hypothesis: dramatically faster than value filter."""
    gene_list = "[" + ", ".join(f"'{g}'" for g in GENES) + "]"
    logging.info(
        "[fetch] get_anndata with %d obs_coords + ADORA var filter — testing coord pushdown",
        len(soma_joinids),
    )
    t0 = time.monotonic()
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        measurement_name="RNA",
        obs_coords=soma_joinids,
        var_value_filter=f"feature_name in {gene_list}",
        obs_column_names=OBS_COLS_FETCH,
    )
    dt = time.monotonic() - t0
    nnz = int(adata.X.nnz) if hasattr(adata.X, "nnz") else int((adata.X != 0).sum())
    rate = adata.n_obs / dt if dt > 0 else float("inf")
    logging.info(
        "[fetch] done: %d cells × %d genes in %.1fs (%.0f cells/s, nnz=%d)",
        adata.n_obs, adata.n_vars, dt, rate, nnz,
    )
    emit_stat({
        "op": "fetch", "cells": int(adata.n_obs), "genes": int(adata.n_vars),
        "nnz": nnz, "seconds": dt,
    })
    return adata


def write_atomic(adata: ad.AnnData, path: Path) -> None:
    """Atomic write so a kill never leaves a half-written .h5ad behind."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    adata.write_h5ad(tmp)
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cells-per-cell-type", type=int, default=1000,
                    help="Stratification cap (default: 1000)")
    ap.add_argument("--seed", type=int, default=42,
                    help="Random seed for stratified sample (default: 42)")
    ap.add_argument("--resume", action="store_true",
                    help="Reuse cached obs scan and final h5ad if present")
    ap.add_argument("--rescan", action="store_true",
                    help="Force a fresh global obs scan (overrides --resume for phase 1)")
    args = ap.parse_args()

    setup_logging()
    logging.info("Cache dir:           %s", CACHE_DIR)
    logging.info("Cells per cell_type: %d", args.cells_per_cell_type)
    logging.info("Random seed:         %d", args.seed)
    logging.info("Resume:              %s", args.resume)
    logging.info("Rescan obs:          %s", args.rescan)

    if FINAL_OUT.exists() and args.resume:
        size_mb = FINAL_OUT.stat().st_size / 1e6
        logging.info(
            "[main] %s already exists (%.1f MB). Nothing to do; "
            "delete it or drop --resume to redo.", FINAL_OUT.name, size_mb,
        )
        return

    census = open_census()
    try:
        obs = scan_human_obs(census, use_cache=not args.rescan)
        sample = stratified_sample(obs, args.cells_per_cell_type, args.seed)
        coords = sorted(sample["soma_joinid"].astype(int).tolist())
        adata = fetch_by_coords(census, coords)
        write_atomic(adata, FINAL_OUT)
        size_mb = FINAL_OUT.stat().st_size / 1e6
        logging.info(
            "[main] wrote %s shape=%s (%.1f MB)",
            FINAL_OUT.name, adata.shape, size_mb,
        )
        emit_stat({
            "op": "final", "cells": int(adata.n_obs),
            "genes": int(adata.n_vars), "size_mb": size_mb,
        })
    finally:
        logging.info("Done")


if __name__ == "__main__":
    main()
