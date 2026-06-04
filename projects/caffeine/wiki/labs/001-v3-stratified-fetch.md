# Lab 001: V3 Stratified Fetch — Walkthrough

## Summary

V3 is the current implementation of Lab 001's data fetch. It builds a balanced cross-tissue cross-cell-type sample of human ADORA expression by paying one big obs scan, sampling locally, and making one large coordinate-based X read.

For *why* this shape, see [001 fetch stall post-mortem](001-fetch-stall-postmortem.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) (especially the Pushdown section), and [Census source H5ADs](../concepts/census-source-h5ads.md) (the alternative path).

This page walks through what v3 actually does, phase by phase, with the empirical numbers from the 2026-06-03 run.

## The Hypothesis V3 Is Testing

V1 chunked by `dataset_id` and v2 chunked by `cell_type` using `obs_value_filter`. Both forms filter on TileDB **attributes**, which have no index — the storage layer must walk every fragment and scan every row. The v2 probe (2026-06-01) confirmed this floor: a query returning 11 cells took 10h 39m.

V3 filters on the obs **dimension** (`soma_joinid`) via `obs_coords`. Dimensions are indexed; TileDB can in principle seek directly to the requested coordinates and skip fragments that contain none of them. See [pushdown](../concepts/tiledb-soma-storage.md#pushdown--why-dimensions-are-fast-and-attributes-are-slow) for the mechanism.

The catch: pushdown only helps if the coordinates cluster. V3's 735,583 sampled coordinates are spread across the entire 62.5M-cell range, so most fragments still overlap at least one coordinate. The hoped-for savings come mostly from intra-fragment seeking (within a fragment, skip non-matching rows) and from the var axis (4 of ~60,000 genes).

## Phase 1 — Global obs scan

### Query

```python
census["census_data"]["homo_sapiens"].obs.read(
    column_names=["soma_joinid", "cell_type", "tissue_general"],
    value_filter="is_primary_data == True and disease == 'normal'",
).concat().to_pandas()
```

| What it reads | Why |
|---|---|
| obs DataFrame for `homo_sapiens` | the cell roster |
| value_filter on `is_primary_data` and `disease` | exclude duplicate cells across studies; keep normal cells only |
| three columns: `soma_joinid`, `cell_type`, `tissue_general` | the minimum needed for stratified sampling. Note: `dataset_id` is intentionally not in this scan — adding it would not change cost but would inflate the output by tens of MB |

### Cache file

`cache/_obs_human_primary_normal.parquet` — paid once per Census version. Future runs reuse it via `--resume` (the script defaults to reading from cache unless `--rescan` is passed).

### Observed cost (2026-06-03)

| | Value |
|---|---|
| wall-clock | **5h 51m 10s** (21,070s) |
| cells returned | 62,464,283 |
| cell_types | 856 |
| tissues | 61 |
| parquet size on disk | 302.7 MB |

Notable: the brain-only v2 scan (18.5M cells, same value-filter shape) took 5h 59m — *the same time* as the full-human scan. Obs scan cost is dominated by the fragment walk, not by how many cells survive the filter. Adding the global filter saved nothing.

## Phase 2 — Stratified sample

### Logic

```python
sampled = (
    obs.groupby("cell_type", group_keys=False, observed=True)
       .apply(lambda g: g.sample(n=min(len(g), 1000), random_state=42))
       .reset_index(drop=True)
)
```

For each `cell_type`, take up to 1000 cells (random, seed 42). Cell_types with fewer than 1000 contribute all their cells. The default cap is configurable via `--cells-per-cell-type`.

### Cache file

`cache/_sample_metadata.parquet` — the row identities of the cells we are about to fetch. Useful as a checkpoint and as the provenance record for downstream analysis.

### Observed cost (2026-06-03)

| | Value |
|---|---|
| wall-clock | 3.7s (pure pandas) |
| cells sampled | **735,583** |
| cell_types represented | 856 (all) |
| parquet size on disk | 3.8 MB |

### Tissue distribution of the sample

Top 15 of 61 tissues, by sampled cell count:

| tissue_general | cells |
|---|---:|
| brain | 100,988 |
| eye | 66,356 |
| lung | 58,510 |
| blood | 53,544 |
| kidney | 41,720 |
| bone marrow | 33,379 |
| small intestine | 32,469 |
| liver | 24,914 |
| colon | 23,588 |
| heart | 21,356 |
| endocrine gland | 19,514 |
| breast | 19,338 |
| respiratory system | 17,744 |
| skin of body | 16,875 |
| musculature | 16,063 |

Brain dominates because it has the most cell_types — the 1000/cell_type cap fills more slots when there are more cell_types to fill.

## Phase 3 — Coordinate-based X read

### Query

```python
cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    measurement_name="RNA",
    obs_coords=[12, 47, 891, ...],          # the 735,583 sampled soma_joinids
    var_value_filter="feature_name in ['ADORA1','ADORA2A','ADORA2B','ADORA3']",
    obs_column_names=[
        "soma_joinid", "cell_type", "tissue", "tissue_general",
        "assay", "dataset_id", "donor_id",
    ],
)
```

### Why this is the architectural test

- `obs_coords` is the **dimension pushdown** form. TileDB knows the row IDs and can use the index.
- `var_value_filter` is still attribute-based — TileDB cannot push it to the fragment selector; it shrinks the AnnData returned but not the column tiles read.

If we wanted symmetric pushdown on both axes, we would resolve gene symbols to `soma_joinid` integers first and pass `var_coords`. That is the obvious follow-up if v3 Phase 3 turns out to be column-read-bound.

### What `obs_column_names` brings along

The richer obs metadata is included so the resulting AnnData supports per-cell provenance:

| Column | Why we keep it |
|---|---|
| `soma_joinid` | join key back to `_obs_human_primary_normal.parquet` |
| `cell_type` | primary grouping variable for the downstream analysis |
| `tissue` | finer than `tissue_general` |
| `tissue_general` | what we stratified on |
| `assay` | which scRNA assay produced the cell |
| `dataset_id` | provenance back to the source study and its H5AD |
| `donor_id` | for later donor-aware QC |

### Observed cost

*Pending — Phase 3 still running at the time of writing. To be filled in once the fetch returns or is aborted.*

## Phase 4 — Atomic write

### Output file

`cache/adora_stratified.h5ad` — the deliverable.

The script writes to a `.tmp` first then `os.replace()`s into place so a kill at any moment leaves either the old file or nothing at all, never a half-written file.

### Expected schema

| AnnData slot | Shape | Contents |
|---|---|---|
| `adata.X` | up to 735,583 × 4 | sparse expression for the four ADORA receptors |
| `adata.obs` | up to 735,583 × 7 | the columns listed above |
| `adata.var` | 4 × N | gene metadata; index is `feature_id` (Ensembl), `feature_name` is the symbol |

Expected on-disk size is modest — sparse 735k × 4 with low nonzero rate compresses to well under 100 MB.

## What This Dataset Enables

- ADORA expression by cell_type across all human tissues at single-cell resolution, with 1000 cells per cell_type ceiling.
- Cross-tissue contrast for cell_types that appear in multiple tissues (e.g., macrophages in liver vs lung vs brain).
- Provenance traceable per cell via `dataset_id` joined back to `census["census_info"]["datasets"]`.

## What It Does Not Cover

| Limitation | Mitigation |
|---|---|
| **Receptors only.** No downstream caffeine pathway genes (CYP1A2, AHR, CREB1, NFAT family, PDE4B/D, RYR1/2/3, HDAC4/5, etc.) | Phase 3 can be re-run with an expanded `GENES` list. The obs scan and sample caches are reusable; only the X read repeats. |
| **One sample per cell_type.** Cells from many studies and donors get pooled before sampling. | The sample metadata records `dataset_id` and `donor_id`. A per-donor or per-dataset breakdown is possible downstream. |
| **No per-dataset attribution at sample time.** `dataset_id` is not in the obs-scan cache. | Either re-scan obs with `dataset_id` (another ~6h), or accept that attribution comes from the resulting AnnData's `obs.dataset_id` column. |
| **Census 2025-11-08 only.** Pinned by the active stable release at run time. | Future runs against a new release reuse the script but must redo Phase 1. |

## File Map

After a complete v3 run, `lab/001_adora_expression/cache/` contains:

| File | Source | When to delete |
|---|---|---|
| `_obs_human_primary_normal.parquet` | Phase 1, ~300 MB | only on Census version change or schema change |
| `_sample_metadata.parquet` | Phase 2, ~4 MB | only on sample-config change (cap or seed) |
| `adora_stratified.h5ad` | Phase 4, < 100 MB | only when re-fetching with a different gene set |
| `fetch.log` | every phase, human | never; cumulative across runs |
| `stats.jsonl` | every phase, structured | never; cumulative across runs |
| `_enumerate_brain.json` | leftover from v2 brain probe | safe to delete; v3 does not use it |

## Related Pages

[001 ADORA expression](001-adora-expression.md), [001 data flow](001-data-flow.md), [001 CellxGene Census API](001-cellxgene-census-api.md), [001 H5AD and AnnData cache](001-h5ad-anndata-cache.md), [001 fetch stall post-mortem](001-fetch-stall-postmortem.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [Census source H5ADs](../concepts/census-source-h5ads.md), [SOMA axes and X](../concepts/soma-axes-and-x.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)
