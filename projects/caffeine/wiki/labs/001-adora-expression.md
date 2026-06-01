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
5. This page

If you are picking the lab back up after the 2026-05-31 fetch stall, also read:

6. [001 fetch stall post-mortem](001-fetch-stall-postmortem.md)
7. [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md)
8. [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)

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

`fetch_adora_cache.py` fetches expression data by tissue and dataset ID, writes partial checkpoints, retries failed requests, and combines per-tissue caches.

Default tissues:

- brain
- heart
- liver
- adipose tissue
- kidney
- blood
- intestine
- lung

## Planned Outputs

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

Related pages: [CellxGene Census API](001-cellxgene-census-api.md), [data flow](001-data-flow.md), [H5AD and AnnData cache](001-h5ad-anndata-cache.md), [notebook guide](001-notebook-guide.md), [adenosine receptors](../concepts/adenosine-receptors.md), [public data landscape](../concepts/public-data-landscape.md), [lab 001 source](../raw/lab-001-source.md), [001 fetch stall post-mortem](001-fetch-stall-postmortem.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)
