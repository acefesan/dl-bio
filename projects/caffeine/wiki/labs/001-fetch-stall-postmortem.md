# Lab 001: Fetch Stall Post-Mortem

## Summary

On 2026-05-31 the first overnight run of `fetch_adora_cache.py --tissues brain` completed two of brain's 186 source datasets in ~21 hours of wall-clock time and was killed before the first on-disk checkpoint fired. Projected at that rate, brain alone would have required ~54 days. The run produced no usable cache and lost ~21 hours of in-memory work. This page is the lessons-learned record.

The raw evidence (log, process snapshot, network snapshot, disk snapshot) lives in [raw lab-001 stall post-mortem](../raw/lab-001-stall-postmortem.md). The conceptual explanation lives in [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md). The measurement plan for the next run lives in [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md).

## What Happened

A single invocation of the fetch script for the brain tissue progressed as follows:

| Step | Elapsed | Outcome |
|---|---:|---|
| open Census, list brain datasets | 7h 02m | 186 datasets, 18.5M cells |
| dataset 1 | 11h 14m | completed in memory, never checkpointed |
| dataset 2 | 3h 07m | completed in memory, never checkpointed |
| dataset 3 | partial | started, killed before completion |
| any `.h5ad` written? | — | no |

No retry warnings appeared in the log. Each individual `cellxgene_census.get_anndata` call was simply slow, not failing-and-retrying.

## What This Confirms

1. The slowness is **inside the SOMA query path**, not in the script's exception handler or our local code.
2. The per-`dataset_id` `get_anndata` call has cost approximately proportional to the cells in the dataset (and possibly to the genes too), even though the result keeps only four genes. See [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) for why.
3. The script's checkpoint cadence (every 3 datasets) is too coarse for individual-dataset times measured in hours. Any kill, crash, or reboot before the third dataset completes loses everything.
4. The "low instantaneous downlink" signal observed during the run is **consistent with**, not contradicted by, a working download. SOMA fetches alternate between bursty S3 reads and CPU-heavy decompression and filtering. See [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md) for the categories of evidence that would have settled this faster.

## What This Does Not Confirm

- The actual ratio of "S3 bytes transferred" to "AnnData bytes returned". We did not instrument the run. The 50.3 GB lifetime `eth0` counter cannot be attributed to the fetch in isolation.
- Whether a network problem (DNS, transient S3 latency, packet loss) contributed on top of the layout cost. Possible but unproven.
- Whether other configurations (`init_buffer_bytes`, `vfs.s3.*` knobs, region) would have meaningfully reduced wall-clock time. The current script already sets generous timeouts and a 512 MB initial buffer.

## Decisions That Follow

The next iteration of the fetch script should change at least the following:

1. **Drop dataset_id chunking.** Replace it with finer-grained obs-axis chunking (cell_type within tissue, or batched `soma_joinid` ranges via the native iterator). See [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) and [001 data flow](001-data-flow.md).
2. **Checkpoint per chunk, not every N chunks.** Each completed chunk should immediately append to a per-tissue cache file (or write a separate chunk file) so a kill at any moment loses at most the work in flight.
3. **Instrument the run.** Per the [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md) plan: `pidstat` time series, `nethogs` for per-process bandwidth, `py-spy dump` if the process appears stuck, an in-script per-chunk timer and byte counter logged at INFO.
4. **Reconsider the scope.** If even cell-type-chunked Census queries are unworkable, fall back to pre-published atlas H5ADs that already materialize cells of interest. The smallest viable Lab 001 deliverable is not "all eight tissues" but "one tissue with a publishable dotplot", and brain is the right tissue if any.

These decisions update the lab plan but do not change the question. See [001 ADORA expression](001-adora-expression.md).

## Status

| Item | Status |
|---|---|
| First overnight brain fetch | aborted, no cache, no result |
| Current `cache/` directory contents | `fetch.log` only |
| `fetch_adora_cache.py` script | unchanged; flagged for rewrite |
| Lab entry `entry.md` | not yet updated to reflect this run |
| Wiki updated | yes (this page, [raw](../raw/lab-001-stall-postmortem.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)) |

Related pages: [001 ADORA expression](001-adora-expression.md), [001 data flow](001-data-flow.md), [001 CellxGene Census API](001-cellxgene-census-api.md), [001 H5AD and AnnData cache](001-h5ad-anndata-cache.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md), [raw post-mortem](../raw/lab-001-stall-postmortem.md)
