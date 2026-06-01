# Lab 001 Source

Source files:

- `../../lab/001_adora_expression/entry.md`
- `../../lab/001_adora_expression/metadata.json`
- `../../lab/001_adora_expression/fetch_adora_cache.py`
- `../../lab/001_adora_expression/explore_adora_expression.ipynb`

## Summary

Lab 001 implements Q1 from the proposal: identify which human cell types express ADORA1, ADORA2A, ADORA2B, and ADORA3 at the highest levels. [ADORA](../concepts/adenosine-receptors.md) is the human adenosine receptor gene-symbol prefix.

The lab uses CellxGene Census as the primary single-cell source and GTEx v8 bulk expression as a secondary sanity check. The current script fetches [ADORA](../concepts/adenosine-receptors.md) expression from normal primary human data, chunked by tissue and dataset, with retries and per-tissue `.h5ad` caches.

## Current Status

Status is `in_progress`. Planned outputs include:

- `figures/dotplot_adora_cell_type.png`
- `figures/ranked_top20_per_receptor.png`
- `pseudobulk_by_cell_type.feather`
- `cross_receptor_overlap.feather`

## Implementation Notes

`fetch_adora_cache.py` is designed for flaky network conditions:

- It queries dataset IDs per tissue.
- It fetches each tissue one dataset at a time.
- It writes partial `.h5ad` checkpoints.
- It combines per-tissue caches into `cache/adora_all_tissues.h5ad`.

Default tissues are brain, heart, liver, adipose tissue, kidney, blood, intestine, and lung.

## Links

Related pages: [001 adora expression](../labs/001-adora-expression.md), [adenosine receptors](../concepts/adenosine-receptors.md), [public data landscape](../concepts/public-data-landscape.md)
