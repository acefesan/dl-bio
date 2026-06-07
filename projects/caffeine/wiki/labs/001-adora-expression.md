# Lab 001: ADORA Expression

## Question

Which human cell types express ADORA1, ADORA2A, ADORA2B, and ADORA3 at the highest levels?

Here, **[ADORA](../concepts/adenosine-receptors.md)** means the human adenosine receptor gene family. These genes encode the A1, A2A, A2B, and A3 receptor proteins that caffeine can block.

This corresponds to Q1 in [research questions](../concepts/research-questions.md).

## Why It Matters

This is the first filter in the [cell type response model](../concepts/cell-type-response-model.md). A cell type without receptor expression is less likely to have a direct physiological caffeine response through adenosine receptor antagonism. A cell type with strong receptor expression becomes a candidate for chromatin accessibility, motif, and GRN follow-up.

## Lab Reading Path

If you are new to the API and data model, read:

1. [CellxGene Census API](001-cellxgene-census-api.md)
2. [Data flow](001-data-flow.md)
3. [H5AD and AnnData cache](001-h5ad-anndata-cache.md)
4. [Notebook guide](001-notebook-guide.md)
5. [Single-cell RNA-seq measurement](../concepts/single-cell-rna-seq-measurement.md)
6. This page
7. [001 ADORA interpretation](001-adora-interpretation.md)

If you are picking the lab back up after the 2026-05-31 fetch stall, also read:

8. [001 fetch stall post-mortem](001-fetch-stall-postmortem.md)
9. [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md)
10. [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)
11. [001 v3 stratified fetch](001-v3-stratified-fetch.md) — phase-by-phase walkthrough of what the current script actually does, with the 2026-06-03 numbers

## Inputs

Primary:

- CellxGene Census via `cellxgene_census`.
- Normal primary human data.
- Genes: ADORA1, ADORA2A, ADORA2B, ADORA3.

Secondary:

- GTEx v8 bulk TPM for tissue-level direction checks.

## Current Files

- `../../lab/001_adora_expression/entry.md`
- `../../lab/001_adora_expression/metadata.json`
- `../../lab/001_adora_expression/fetch_adora_cache.py`
- `../../lab/001_adora_expression/explore_adora_expression.ipynb`

## Script Behavior

`fetch_adora_cache.py` is the v3 stratified-sample fetcher. It builds a balanced dataset of up to N cells per `cell_type` across all primary normal human cells, in three phases:

1. **Global obs scan** (one-time, slow — many hours). Reads `(soma_joinid, cell_type, tissue_general)` for every primary normal human cell. Cached as `cache/_obs_human_primary_normal.parquet`. Reused on every subsequent run.
2. **Stratified sample** (seconds). Sample up to N (default 1000) cells per `cell_type` in pandas. Result cached as `cache/_sample_metadata.parquet`.
3. **Coord-based X read** (hypothesis: minutes; unproven). One `get_anndata` call with `obs_coords=<sampled soma_joinids>` and `var_value_filter` for the four ADORA genes. The coord pushdown should bypass the fragment-walk that crushed v1 and v2 — see [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) and [001 fetch stall post-mortem](001-fetch-stall-postmortem.md).

Final output: `cache/adora_stratified.h5ad`.

Why this shape: v2's probe (2026-06-01) confirmed that per-`cell_type` `obs_value_filter` queries pay ~10h per call regardless of result size, so chunking that way is non-viable. The stratified sample plus single coord-based read trades that for one obs scan paid once, then one X read that — *if* the dimension index works as hoped — does not walk every fragment.

## Planned Outputs

Cache artifact (from `fetch_adora_cache.py`):

- `cache/adora_stratified.h5ad` — sparse cells × ADORA receptors AnnData with per-cell `cell_type`, `tissue`, `tissue_general`, `assay`, `dataset_id`, `donor_id` provenance.

Downstream analysis artifacts (from the notebook, once the cache exists):

- [ADORA](../concepts/adenosine-receptors.md) dotplot across cell types.
- Top 20 cell types per receptor.
- Pseudobulk mean expression and percent-expressing tables.
- Cross-receptor overlap table.

## Interpretation Template

When results exist, fill in:

- top cell types per receptor,
- receptor co-expression patterns,
- tissue-specific surprises,
- disagreement with GTEx bulk,
- candidates for Q9 accessibility follow-up.

Related pages: [001 ADORA interpretation](001-adora-interpretation.md), [single-cell RNA-seq measurement](../concepts/single-cell-rna-seq-measurement.md), [001 v3 stratified fetch](001-v3-stratified-fetch.md), [CellxGene Census API](001-cellxgene-census-api.md), [data flow](001-data-flow.md), [H5AD and AnnData cache](001-h5ad-anndata-cache.md), [notebook guide](001-notebook-guide.md), [adenosine receptors](../concepts/adenosine-receptors.md), [public data landscape](../concepts/public-data-landscape.md), [lab 001 source](../raw/lab-001-source.md), [001 fetch stall post-mortem](001-fetch-stall-postmortem.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)
