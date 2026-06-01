"""Fetch ADORA expression caches for Q1 — v2, post-stall rewrite.

Background — read these first:
    wiki/labs/001-fetch-stall-postmortem.md
    wiki/concepts/tiledb-soma-storage.md
    wiki/concepts/network-and-io-instrumentation.md

The v1 of this script chunked by dataset_id and never produced a usable cache
(see post-mortem). The hypothesis is that the Census X matrix is partitioned
cell-major, so a per-dataset_id slice still walks fragments containing every
gene, while the four-gene var_value_filter only shrinks the answer. We do not
yet have an instrumented run that proves a faster shape works, so this rewrite
prioritizes evidence over throughput:

    1. Chunk by (tissue, cell_type) instead of (tissue, dataset_id). Cell types
       are the analysis grouping we care about anyway, and a cell_type filter
       reduces both the answer size and the cell range the query walks.
    2. Checkpoint per cell_type, not per N. Each completed cell_type writes
       immediately to its own .h5ad. A kill loses at most one cell_type's work.
    3. Instrument every call. Per-call cells / nnz / seconds / cells-per-second
       go to fetch.log at INFO and to stats.jsonl as one JSON record per call.
    4. Cache the obs enumeration. The previous run spent 7h enumerating brain
       datasets; the same cost applies to cell_type enumeration. We pay it
       once per tissue and reuse a JSON cache.
    5. Probe mode. --probe fetches the smallest cell_type per tissue and
       projects full-tissue time before committing.

Usage:
    # Probe brain to see whether this shape is even viable (timing report)
    python fetch_adora_cache.py --tissues brain --probe

    # Full per-cell-type fetch (resumable)
    python fetch_adora_cache.py --tissues brain --resume

    # Cap individual cell_type queries (skip mega cell_types)
    python fetch_adora_cache.py --tissues brain --max-cells-per-chunk 200000 --resume

    # Multi-tissue overnight run after probe shows it is viable
    python fetch_adora_cache.py --tissues brain heart liver --resume

Output:
    cache/<tissue>/<cell_type_slug>.h5ad         per-cell-type cache (atomic)
    cache/<tissue>.h5ad                          concatenated per-tissue cache
    cache/adora_all_tissues.h5ad                 concatenated cross-tissue cache
    cache/_enumerate_<tissue>.json               cached cell_type enumeration
    cache/fetch.log                              human log (appended)
    cache/stats.jsonl                            per-call metrics (one JSON per line)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from pathlib import Path

import anndata as ad
import cellxgene_census

GENES = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
OBS_COLS = ["cell_type", "tissue", "tissue_general", "assay", "dataset_id", "donor_id"]
DEFAULT_TISSUES = [
    "brain",
    "heart",
    "liver",
    "adipose tissue",
    "kidney",
    "blood",
    "intestine",
    "lung",
]

SCRIPT_DIR = Path(__file__).parent
CACHE_DIR = SCRIPT_DIR / "cache"
STATS_FILE = CACHE_DIR / "stats.jsonl"
SLUG_RE = re.compile(r"[^a-z0-9]+")


def slugify(s: str) -> str:
    return SLUG_RE.sub("_", s.lower()).strip("_") or "unnamed"


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
    ctx = cellxgene_census.get_default_soma_context(
        tiledb_config={
            "vfs.s3.connect_timeout_ms": "60000",
            "vfs.s3.request_timeout_ms": "300000",
            "soma.init_buffer_bytes": str(512 * 1024 * 1024),
        }
    )
    return cellxgene_census.open_soma(context=ctx)


def list_cell_types_for_tissue(census, tissue: str, *, use_cache: bool) -> list[tuple[str, int]]:
    """Return [(cell_type, n_cells)] sorted ascending. Cached as JSON per tissue."""
    cache_path = CACHE_DIR / f"_enumerate_{slugify(tissue)}.json"

    if use_cache and cache_path.exists():
        data = json.loads(cache_path.read_text())
        logging.info(
            "[%s] enumerate cached: %d cell_types (from %s)",
            tissue, len(data), cache_path.name,
        )
        return [(ct, int(n)) for ct, n in data]

    logging.info("[%s] enumerate cell_types via obs.read — expect minutes to hours", tissue)
    t0 = time.monotonic()
    obs = (
        census["census_data"]["homo_sapiens"]
        .obs.read(
            column_names=["cell_type"],
            value_filter=(
                f"tissue_general == '{tissue}' "
                f"and is_primary_data == True "
                f"and disease == 'normal'"
            ),
        )
        .concat()
        .to_pandas()
    )
    dt = time.monotonic() - t0
    counts = obs["cell_type"].value_counts().sort_values(ascending=True)
    pairs = [(str(ct), int(n)) for ct, n in counts.items()]

    cache_path.write_text(json.dumps(pairs, indent=2))
    logging.info(
        "[%s] enumerate done: %d total cells, %d cell_types, %.0fs (cached to %s)",
        tissue, len(obs), len(pairs), dt, cache_path.name,
    )
    emit_stat({
        "op": "enumerate", "tissue": tissue, "cells": int(len(obs)),
        "cell_types": len(pairs), "seconds": dt,
    })
    return pairs


def fetch_one_cell_type(
    census, tissue: str, cell_type: str, retries: int, backoff: int
) -> ad.AnnData | None:
    """Fetch ADORA expression for one (tissue, cell_type). Retries on any exception."""
    safe_ct = cell_type.replace("'", "''")
    obs_filter = (
        f"tissue_general == '{tissue}' "
        f"and cell_type == '{safe_ct}' "
        f"and is_primary_data == True "
        f"and disease == 'normal'"
    )
    gene_list = "[" + ", ".join(f"'{g}'" for g in GENES) + "]"
    last_err: Exception | None = None

    for attempt in range(1, retries + 1):
        t0 = time.monotonic()
        try:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                measurement_name="RNA",
                var_value_filter=f"feature_name in {gene_list}",
                obs_value_filter=obs_filter,
                obs_column_names=OBS_COLS,
            )
            dt = time.monotonic() - t0
            n_cells = int(adata.n_obs)
            nnz = int(adata.X.nnz) if hasattr(adata.X, "nnz") else int((adata.X != 0).sum())
            rate = n_cells / dt if dt > 0 else float("inf")
            logging.info(
                "    [%s] %s: %d cells in %.1fs (%.0f cells/s, nnz=%d)",
                tissue, cell_type, n_cells, dt, rate, nnz,
            )
            emit_stat({
                "op": "fetch", "tissue": tissue, "cell_type": cell_type,
                "cells": n_cells, "nnz": nnz, "seconds": dt,
            })
            return adata
        except Exception as e:
            last_err = e
            wait = backoff * (2 ** (attempt - 1))
            logging.warning(
                "    [%s] %s attempt %d/%d failed (%s): sleeping %ds",
                tissue, cell_type, attempt, retries, type(e).__name__, wait,
            )
            emit_stat({
                "op": "fetch_retry", "tissue": tissue, "cell_type": cell_type,
                "attempt": attempt, "error": type(e).__name__,
            })
            time.sleep(wait)

    logging.error(
        "    [%s] %s: giving up after %d attempts: %s",
        tissue, cell_type, retries, last_err,
    )
    emit_stat({
        "op": "fetch_failed", "tissue": tissue, "cell_type": cell_type,
        "error": repr(last_err),
    })
    return None


def write_atomic(adata: ad.AnnData, path: Path) -> None:
    """Write h5ad to a tmp file then rename. Avoids half-written files on kill."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    adata.write_h5ad(tmp)
    os.replace(tmp, path)


def fetch_tissue(
    census,
    tissue: str,
    retries: int,
    backoff: int,
    max_cells_per_chunk: int | None,
    resume: bool,
) -> None:
    safe_tissue = slugify(tissue)
    tissue_dir = CACHE_DIR / safe_tissue
    tissue_dir.mkdir(parents=True, exist_ok=True)
    tissue_h5ad = CACHE_DIR / f"{safe_tissue}.h5ad"

    if tissue_h5ad.exists() and resume:
        logging.info("[%s] tissue cache exists, skipping: %s", tissue, tissue_h5ad)
        return

    cell_types = list_cell_types_for_tissue(census, tissue, use_cache=True)
    if not cell_types:
        logging.warning("[%s] no cell_types matched filter", tissue)
        return

    t_start = time.monotonic()
    for i, (ct, n_cells) in enumerate(cell_types, 1):
        ct_slug = slugify(ct)
        ct_path = tissue_dir / f"{ct_slug}.h5ad"

        if ct_path.exists() and resume:
            logging.info(
                "  [%s] [%d/%d] %s (n=%d) — cached, skip",
                tissue, i, len(cell_types), ct, n_cells,
            )
            continue

        if max_cells_per_chunk is not None and n_cells > max_cells_per_chunk:
            logging.info(
                "  [%s] [%d/%d] %s (n=%d) > --max-cells-per-chunk=%d, skip",
                tissue, i, len(cell_types), ct, n_cells, max_cells_per_chunk,
            )
            continue

        logging.info(
            "  [%s] [%d/%d] %s (n=%d) — fetching",
            tissue, i, len(cell_types), ct, n_cells,
        )
        adata = fetch_one_cell_type(census, tissue, ct, retries, backoff)
        if adata is None or adata.n_obs == 0:
            continue
        write_atomic(adata, ct_path)
        size_mb = ct_path.stat().st_size / (1024 * 1024)
        logging.info("    [%s] %s wrote %s (%.1f MB)", tissue, ct, ct_path.name, size_mb)

    t_total = time.monotonic() - t_start
    logging.info("[%s] all cell_types attempted in %.0fs", tissue, t_total)

    pieces = []
    for ct, _ in cell_types:
        p = tissue_dir / f"{slugify(ct)}.h5ad"
        if p.exists():
            pieces.append(ad.read_h5ad(p))
    if pieces:
        combined = ad.concat(pieces, axis=0, merge="same")
        write_atomic(combined, tissue_h5ad)
        logging.info(
            "[%s] wrote tissue cache: %s shape=%s",
            tissue, tissue_h5ad.name, combined.shape,
        )
        emit_stat({
            "op": "tissue_done", "tissue": tissue,
            "cells": int(combined.n_obs), "cell_types": len(pieces),
            "seconds": t_total,
        })
    else:
        logging.warning("[%s] no per-cell_type pieces to combine", tissue)


def build_combined_cache() -> None:
    paths = sorted(
        p for p in CACHE_DIR.glob("*.h5ad")
        if p.name != "adora_all_tissues.h5ad"
    )
    if not paths:
        logging.warning("No per-tissue caches to combine")
        return
    logging.info("Combining %d per-tissue caches", len(paths))
    combined = ad.concat([ad.read_h5ad(p) for p in paths], axis=0, merge="same")
    out = CACHE_DIR / "adora_all_tissues.h5ad"
    write_atomic(combined, out)
    logging.info("Combined: shape=%s → %s", combined.shape, out.name)


def probe(census, tissue: str, retries: int, backoff: int) -> None:
    """Fetch the smallest cell_type in a tissue and project full-tissue time."""
    cell_types = list_cell_types_for_tissue(census, tissue, use_cache=True)
    if not cell_types:
        logging.warning("[%s] probe: no cell_types", tissue)
        return
    # Smallest cell_type with at least 10 cells; fall back to absolute smallest
    candidates = [(ct, n) for ct, n in cell_types if n >= 10]
    ct, n = (candidates or cell_types)[0]
    total = sum(n for _, n in cell_types)

    logging.info(
        "[%s] probe: fetching smallest viable cell_type '%s' (n=%d of %d total)",
        tissue, ct, n, total,
    )
    t0 = time.monotonic()
    adata = fetch_one_cell_type(census, tissue, ct, retries, backoff)
    dt = time.monotonic() - t0
    if adata is None:
        logging.error("[%s] probe: fetch failed", tissue)
        return

    # Linear projection: assumes per-cell time is roughly constant. It will not be —
    # large cell_types should be cheaper per cell — but it gives an upper bound.
    proj_s = (total / n) * dt if n > 0 else float("inf")
    logging.info(
        "[%s] probe done: %d cells in %.1fs. "
        "Linear upper-bound projection for full tissue (%d cells) = %.0fs = %.2fh",
        tissue, adata.n_obs, dt, total, proj_s, proj_s / 3600,
    )
    emit_stat({
        "op": "probe", "tissue": tissue, "cell_type": ct,
        "cells": int(adata.n_obs), "seconds": dt,
        "total_cells": total, "projected_seconds": proj_s,
    })


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--tissues", nargs="+", default=DEFAULT_TISSUES)
    ap.add_argument("--retries", type=int, default=5)
    ap.add_argument(
        "--backoff", type=int, default=30,
        help="Initial backoff seconds; doubles per attempt",
    )
    ap.add_argument(
        "--max-cells-per-chunk", type=int, default=None,
        help="Skip cell_types larger than this — safety valve for mega cell_types",
    )
    ap.add_argument(
        "--resume", action="store_true",
        help="Skip cell_types and tissues that already have .h5ad on disk",
    )
    ap.add_argument(
        "--probe", action="store_true",
        help="Per tissue, fetch only the smallest cell_type and report projected time",
    )
    args = ap.parse_args()

    setup_logging()
    logging.info("Cache dir: %s", CACHE_DIR)
    logging.info("Tissues:   %s", args.tissues)
    logging.info("Mode:      %s", "probe" if args.probe else "full")
    if args.max_cells_per_chunk:
        logging.info("Cap:       --max-cells-per-chunk=%d", args.max_cells_per_chunk)

    census = open_census()
    try:
        for t in args.tissues:
            try:
                if args.probe:
                    probe(census, t, args.retries, args.backoff)
                else:
                    fetch_tissue(
                        census, t, args.retries, args.backoff,
                        args.max_cells_per_chunk, args.resume,
                    )
            except Exception:
                logging.exception("[%s] unhandled error", t)

        if not args.probe:
            build_combined_cache()
    finally:
        logging.info("Done")


if __name__ == "__main__":
    main()
