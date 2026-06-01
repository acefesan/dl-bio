# Lab 001: Notebook Guide

## Summary

`explore_adora_expression.ipynb` is a learning notebook. Its job is to teach the `cellxgene_census` API and the shape of the data before turning the workflow into a more robust script.

The notebook is not just "analysis." It is a guided tour of:

- opening the Census,
- inspecting metadata,
- filtering cells,
- filtering genes,
- loading the cached `AnnData` subset,
- aggregating expression by cell type,
- plotting receptor expression.

## Notebook Sections

| Section | What it teaches |
|---|---|
| 0. Install dependencies | packages needed for Census access and plotting |
| 1. What is the Census? | how to open the Census object with `open_soma()` |
| 2. Peek at the summary | how to inspect global Census metadata |
| 3. The cell metadata `obs` | how cell-level fields like tissue and cell type are stored |
| 4. The gene metadata `var` | how gene symbols and feature IDs are stored |
| 5. Loading the cached AnnData subset | how to load `.h5ad` files created by the fetch script |
| 6. Understanding the AnnData | how to inspect `.X`, `.obs`, and `.var` |
| 7. Normalize + map Ensembl to symbol | how to make expression easier to compare and label |
| 8. Per-cell-type aggregation | how to summarize receptor expression by cell type |
| 9. Dotplot | how to visualize mean expression and percent expressing |
| 10. Save intermediates | how to write reusable outputs |

## Core Notebook Objects

| Object | Meaning |
|---|---|
| `census` | open handle to the remote Census data object |
| `summary` | high-level Census summary table |
| `datasets` | dataset metadata table |
| `human` | human SOMA experiment inside the Census |
| `var_df` | gene metadata table |
| `adora` | gene metadata rows for ADORA genes |
| `adata` | selected cells x selected genes expression object |

## Why the Script Exists Too

The notebook is great for learning, but live Census queries can be fragile over large tissues. `fetch_adora_cache.py` exists because:

- tissue-wide queries can stream many cells from cloud storage,
- network/DNS failures can interrupt long reads,
- per-dataset chunking makes retries practical,
- `.h5ad` caches make the notebook faster after the first run.

So the notebook teaches the API and analyzes cached data; the script makes the workflow repeatable and keeps large downloads out of the interactive VS Code session.

## Cache-First Rule

Section 5 now reads local `.h5ad` cache files with `anndata.read_h5ad()`. It checks for `adora_all_tissues.h5ad` first, then falls back to `adora_brain.h5ad`. Because VS Code notebooks can start kernels from different working directories, it checks both `cache/` and `lab/001_adora_expression/cache/`.

If neither file exists, the notebook raises a `FileNotFoundError` telling you to run `fetch_adora_cache.py`. That is intentional. The notebook should not silently launch a large live Census download.

## What to Read Before Running

Read in this order:

1. [CellxGene Census API](001-cellxgene-census-api.md)
2. [Data flow](001-data-flow.md)
3. [H5AD and AnnData cache](001-h5ad-anndata-cache.md)
4. This page
5. `../../lab/001_adora_expression/fetch_adora_cache.py`
6. `../../lab/001_adora_expression/explore_adora_expression.ipynb`

Related pages: [Lab 001 overview](001-adora-expression.md), [H5AD and AnnData cache](001-h5ad-anndata-cache.md), [ADORA receptors](../concepts/adenosine-receptors.md), [cAMP signaling](../concepts/camp-signaling.md)
