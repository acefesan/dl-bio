# 001 ADORA Interpretation

## Summary

Lab 001 asks which human cell populations express the four adenosine receptor genes:

```text
ADORA1
ADORA2A
ADORA2B
ADORA3
```

The current local work uses the full Tabula Sapiens H5AD plus a compact extracted ADORA matrix:

```text
lab/001_adora_expression/cache/tabula_sapiens_all_cells.h5ad
lab/001_adora_expression/cache/tabula_sapiens_adora_expression.npz
```

Interpret the results as **single-cell RNA evidence for receptor transcript**, not as direct receptor-protein abundance or caffeine-response magnitude.

## What The Current Tabula Sapiens Run Measures

The source data are single-cell transcriptomes. Each row is a cell, each gene column is an RNA feature, and each ADORA value is the atlas's processed expression value for that gene in that cell.

Important local facts:

| Item | Value |
|---|---|
| Cells | 1,136,218 |
| Genes | 60,606 |
| ADORA genes extracted | ADORA1, ADORA2A, ADORA2B, ADORA3 |
| UMAP dimensions | 2 |
| Broad tissue labels | 28 |
| Cell-type labels | 180 |
| Broad tissue `Brain` present? | No |

Because brain is absent from `tissue_in_publication` and `tissue`, this file is useful for many peripheral tissues but is not a brain ADORA atlas.

## Current ADORA Signal

From the full Tabula Sapiens H5AD, the current extracted nonzero counts are:

| Gene | Cells with expression > 0 |
|---|---:|
| ADORA1 | 8,416 |
| ADORA2A | 3,223 |
| ADORA2B | 48,403 |
| ADORA3 | 22,399 |

This already tells the main visual story: ADORA expression is sparse. A regular UMAP expression layer will be visually dominated by cells with zero or near-zero values unless those cells are drawn faintly.

## Thresholds Used So Far

For the high-expression UMAP, "high" means the top quartile among nonzero-expressing cells for each gene:

| Gene | Threshold |
|---|---:|
| ADORA1 | 0.684 |
| ADORA2A | 1.080 |
| ADORA2B | 0.779 |
| ADORA3 | 0.943 |

This rule asks:

```text
among cells where this receptor was detected, which cells are relatively high?
```

It intentionally does not ask "which cells are nonzero?" because nonzero cells are already rare and vary by receptor.

For the interactive UMAP, the ADORA expression toggles draw expression below `0.04` as a very faint context layer. That threshold is a visualization floor, not a biological claim. Its job is to stop zero-like cells from hiding the sparse receptor-positive cells.

## How To Read The UMAP

UMAP is a two-dimensional projection of high-dimensional transcriptomes. It helps answer:

```text
Are receptor-positive cells concentrated in recognizable transcriptional neighborhoods?
```

It does not prove:

- clusters are discrete biological types,
- distance between far-apart islands is meaningful,
- a visible blob is driven by ADORA,
- ADORA transcript implies receptor protein,
- receptor-positive cells are caffeine-responsive in vivo.

Use UMAP to find questions. Use tables and dotplots to answer them.

## How To Read The Dotplots

The ADORA dotplots are the main interpretation tool.

Each dot answers one `group x gene` question:

| Dot channel | Meaning |
|---|---|
| Dot size | Percent of cells in that group with expression > 0 |
| Dot color | Mean expression in that group |
| Row | Cell type or tissue |
| Column | ADORA receptor gene |

This separates two different ideas:

| Pattern | Interpretation |
|---|---|
| Big pale dot | Many cells weakly express the gene |
| Small bright dot | Few cells strongly express the gene |
| Big bright dot | Common and strong signal in that group |
| Missing/tiny dark dot | Little observed evidence in this atlas layer |

For sparse receptor genes, dotplots are usually more interpretable than expression UMAPs because they name the cell populations and show prevalence directly.

For the broader set of plot types available — stacked violin, matrixplot, tracksplot, heatmap, density UMAP — and the standard analysis steps that produced the embeddings and labels this file uses (HVG selection, PCA, scVI integration, Leiden clustering, marker finding), see [scRNA visualization and analysis](../concepts/scrna-visualization-and-analysis.md).

## What The Current Dotplot Suggests

The top cell-type dotplot suggests:

| Receptor | Main visible pattern |
|---|---|
| ADORA1 | Pancreatic ductal cells stand out; additional smaller signals appear in eye/cardiac/testis-related cell types |
| ADORA2A | Sparse overall in this Tabula Sapiens file; enriched immune/thymic signals are modest |
| ADORA2B | Broad epithelial/progenitor-like signal, including basal, urothelial, goblet, and other epithelial groups |
| ADORA3 | Macrophage, microglia-like, monocyte, basophil, mast-cell, and other myeloid-associated populations |

Treat these as atlas-derived hypotheses. The next checks are tissue, donor, assay, and source-dataset stratification. For the immune labels in this table, including myeloid dendritic cells, see [immune cell types](../concepts/immune-cell-types.md).

## Tongue ADORA2B Example

The broad `Tongue` tissue dot is mostly an ADORA2B signal:

| Gene | Percent expressing in Tongue |
|---|---:|
| ADORA1 | 0.46% |
| ADORA2A | 0.04% |
| ADORA2B | 19.91% |
| ADORA3 | 0.62% |

Breaking Tongue down by cell type shows that the signal is mainly epithelial:

| Tongue cell type | Cells | ADORA2B percent expressing |
|---|---:|---:|
| basal cell | 16,263 | 35.5% |
| stratified squamous epithelial cell | 12,323 | 14.2% |
| taste receptor cell | 32 | 6.3% |

So the safer interpretation is:

```text
ADORA2B is visible in tongue barrier epithelium,
especially basal and stratified squamous epithelial cells.
```

It is not primarily a taste-cell result. For what basal, stratified squamous, and epithelial cells mean, see [epithelial cell types](../concepts/epithelial-cell-types.md).

## Why Brain Expectations Do Not Match This File

The project proposal correctly expects strong brain relevance for ADORA biology, especially ADORA1 and ADORA2A. The local Tabula Sapiens H5AD does not include a broad or fine brain tissue label.

That means:

- lack of `Brain` in the tissue dotplot is an atlas-coverage fact,
- ADORA2A being sparse here does not contradict striatal A2A biology,
- eye/retina or microglia-like labels should not be treated as a brain substitute,
- a dedicated brain single-cell or single-nucleus atlas is needed for brain-specific ADORA interpretation.

## Interpretation Checklist

For any surprising ADORA-positive blob or cell type:

1. Check `cell_type`, `broad_cell_class`, `tissue_in_publication`, and `tissue`.
2. Check how many cells and donors support the signal.
3. Check whether one assay or donor dominates.
4. Compare percent-expressing and mean expression separately.
5. Check whether the gene was present in source datasets before interpreting zeros.
6. Compare raw, normalized, and decontaminated layers if the result matters.
7. Validate against a tissue-specific atlas if the expected biology is absent from Tabula Sapiens.
8. Treat ADORA mRNA as candidate receptor machinery, not as direct response evidence.

## Follow-Up Analyses

Useful next artifacts:

- receptor x `tissue + cell_type` dotplot,
- donor-stratified ADORA summary,
- assay-stratified ADORA summary,
- source-dataset feature-presence audit,
- dedicated brain atlas ADORA extraction,
- chromatin accessibility check at ADORA loci for top cell types,
- receptor co-expression table: ADORA1-only, ADORA2A-only, ADORA2B-only, ADORA3-only, multi-receptor.

Related pages: [scRNA visualization and analysis](../concepts/scrna-visualization-and-analysis.md), [single-cell RNA-seq measurement](../concepts/single-cell-rna-seq-measurement.md), [epithelial cell types](../concepts/epithelial-cell-types.md), [immune cell types](../concepts/immune-cell-types.md), [001 ADORA expression](001-adora-expression.md), [Census X layers and feature presence](../concepts/census-x-layers-and-feature-presence.md), [Census source H5ADs](../concepts/census-source-h5ads.md), [adenosine receptors](../concepts/adenosine-receptors.md), [cell type response model](../concepts/cell-type-response-model.md)
