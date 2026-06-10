# 001 — The Joint Tabula × HBCA UMAP

## What This Is

Tabula Sapiens and the Human Brain Cell Atlas (HBCA) are two separate downloads, sequenced by different labs with different protocols. Each ships its own UMAP, but those two UMAPs live in **different coordinate systems** — you cannot put a Tabula cell and an HBCA cell on the same axes just because both files contain an `X_umap`. To ask "do HBCA brain glia land near their closest Tabula matches?" you have to build **one shared embedding from scratch** that sees both atlases at once.

That is what `make_shared_umap_tabula_hbca.py` does. This page is the narrative of how it works, so the figures in `figures/shared_umap_tabula_hbca/` are readable as *method*, not just pictures. For UMAP fundamentals underneath this, read [UMAP and dimensionality reduction](../concepts/umap-and-dimensionality-reduction.md) first.

**Script:** `lab/001_adora_expression/make_shared_umap_tabula_hbca.py`
**Outputs:** `figures/shared_umap_tabula_hbca/` (and `_harmony/`, plus `_smoketest` variants)

## Why You Can't Just Stack Two UMAPs

The two atlases share genes but were processed independently. Stacking their published UMAPs would be meaningless — Tabula's UMAP-coordinate `(3.2, -1.1)` and HBCA's `(3.2, -1.1)` have no relationship. The integration problem is: find a *common feature space*, project both atlases into it, remove the atlas-of-origin signal, and only then run one UMAP over the union. Every step below exists to make the two atlases comparable.

## The Pipeline, Step By Step

### 1. Find shared genes, keep the informative ones
Intersect Tabula's `raw/var/feature_name` with HBCA's `var/feature_name`: **56,999 gene symbols** are common to both. You don't want all of them. The script:
- drops Tabula-flagged `feature_is_filtered` genes,
- drops mitochondrial (`MT-`), ribosomal (`RPL`/`RPS`), and spike-in (`ERCC-`) genes — these dominate variance for technical, not biological, reasons,
- ranks the rest by Tabula's stored per-gene standard deviation and keeps the top **2,500**,
- then **force-adds the four ADORA genes** so the receptors of interest are guaranteed in the feature set even if they aren't highly variable.

Result: **2,504 genes** used as the shared feature space. (This is a fixed-panel stand-in for proper HVG selection — see [scRNA visualization and analysis](../concepts/scrna-visualization-and-analysis.md#highly-variable-gene-selection).)

### 2. Subsample, stratified by cell type
A million-cell UMAP is slow and the brain atlas would be outvoted by Tabula's size. The script takes **40,000 cells from each atlas** (80k total), sampled *stratified by cell label* — Tabula by `cell_type`, HBCA by `supercluster_term` — so rare populations are represented rather than swamped by the common ones. Equal 40k/40k draws also keep the two atlases balanced going into batch correction.

### 3. Normalize identically
Both subsampled count matrices get the same treatment: scale each cell's counts to a fixed library size of **10,000** (counts-per-10k), then `log1p`. Identical normalization is what makes a Tabula cell and an HBCA cell numerically comparable. (Tabula is read from `raw/X`, HBCA from `X` — the script pulls raw-count-like layers from each and normalizes them itself rather than trusting two differently-normalized published layers.)

### 4. One linear reduction over the union
The two normalized matrices are `vstack`ed into a single 80k × 2,504 sparse matrix, and **`TruncatedSVD` reduces it to 50 components** (those 50 components capture ~46% of the variance). SVD rather than PCA because the matrix is sparse and centering it would destroy that sparsity — see [the SVD-vs-PCA note](../concepts/umap-and-dimensionality-reduction.md#pca--svd-vs-t-sne-vs-umap). The 50-D latent is then z-scored (`StandardScaler`) so no single component dominates the distance metric.

### 5. Remove the atlas signal (this is the load-bearing step)
After SVD, the biggest axis of variation is often *which atlas a cell came from*, not what kind of cell it is. The script supports three correction modes via `--batch-correction`:

| Mode | What it does | Output dir |
|---|---|---|
| `center` (default) | subtract each atlas's mean from every SVD component | `shared_umap_tabula_hbca/` |
| `harmony` | run Harmony on the SVD components with `dataset` as the batch variable | `shared_umap_tabula_hbca_harmony/` |
| `none` | no correction — cells separate by atlas | (diagnostic) |

`center` is a cheap shift that removes first-order atlas offset. `harmony` is the real integrator: it iteratively pulls matching populations from the two atlases into shared neighborhoods. **The Harmony version is the one to trust** for "do these two atlases agree." The uncorrected/centered versions are useful precisely because they show what *failure* looks like — two atlases drifting into two blobs.

### 6. One UMAP over the corrected latent
`umap.UMAP(n_neighbors=15, min_dist=0.5, metric="cosine", random_state=19)` turns the corrected 50-D latent into 2D. Cosine metric is the usual choice on integrated latents. The seed is fixed (19) so the picture is reproducible — but recall the layout is still only meaningful *locally* ([why](../concepts/umap-and-dimensionality-reduction.md#what-a-umap-picture-means--and-doesnt)).

### 7. Label the ADORA-high cells
Pre-extracted ADORA expression (`*_adora_expression.npz`) is attached per cell. A cell is "high" for a receptor if it clears that receptor's threshold — Tabula uses the per-gene nonzero-quartile thresholds from the single-atlas analysis (e.g. ADORA1 ≥ 0.684), HBCA uses a flat 1.0. Each cell gets one `adora_high_gene` label (its strongest above-threshold receptor, else `none`) for the overlay coloring.

## The Outputs And How To Read Them

Each output directory contains:

| File | Content |
|---|---|
| `shared_umap_tabula_hbca_dataset.png` | joint UMAP colored by **source atlas** (Tabula blue, HBCA red) |
| `shared_umap_tabula_hbca_adora_high.png` | joint UMAP with **ADORA-high cells** drawn over a faint grey background |
| `shared_umap_tabula_hbca_cells.csv` | per-cell row: coordinates, both atlases' labels, ADORA expression + high-gene |
| `shared_umap_tabula_hbca_arrays.npz` | raw 2D embedding, the 50-D latent (pre- and post-correction), indices, gene list |
| `shared_umap_tabula_hbca_genes.txt` | the 2,504 genes used |
| `shared_umap_tabula_hbca_summary.json` | full run config + SVD variance + UMAP params + thresholds |

**Reading the `dataset` figure:** this is the integration QC plot. In the **Harmony** version, well-mixed regions (blue and red interleaved) mean the two atlases found shared cell populations there; regions that stay single-color are populations unique to one atlas (e.g. brain-only glia HBCA has and Tabula doesn't). In the **centered/uncorrected** version, heavy blue/red separation is the expected "batch not removed" look — that contrast is the teaching point.

**Reading the `adora_high` figure:** look for whether ADORA-high cells of a given color *cluster* rather than scatter. Clustering says the receptor marks a coherent population in the shared space; scatter says it's diffuse. Cross-reference the cluster against `*_cells.csv` to name it.

## The `_smoketest` Variants

`shared_umap_tabula_hbca_smoketest/` and `_harmony_smoketest/` are tiny fast runs (far fewer cells) used to confirm the pipeline executes end-to-end before committing to the full 80k run. They are not analysis outputs — ignore them when reading results.

## Caveats (Read Before Concluding Anything)

- **Subsample, not the full atlas.** 40k+40k of ~2M cells. Rare populations are stratified-in but counts are not exhaustive. For prevalence claims, use the full-atlas pseudobulk table, not this UMAP.
- **Fixed-panel, not HVG.** The 2,504 genes are a top-variance panel from Tabula's stored std, not a proper joint HVG selection. Good enough to position cells; not a publication-grade integration.
- **No neurons.** HBCA here is non-neuronal only. The brain's neuronal populations — exactly where ADORA1/2A biology is expected — are absent until the 30 GB neuron H5AD is downloaded. The joint UMAP cannot show what isn't loaded.
- **This is a deliberately practical integration**, as the script's own docstring says: raw-ish counts, library-size norm, SVD, simple centering or Harmony, UMAP. It answers "do the atlases broadly agree and where do ADORA cells sit," not "here is the definitive integrated reference."

## Provenance vs The Native UMAPs

Don't confuse this joint UMAP with the per-atlas published ones. Tabula's native UMAP was built on its `X_scvi` latent with Euclidean metric and scVI batch correction; HBCA ships a `X_UMAP` with no recorded parameters. This page's embedding is a *new* one computed here over both atlases with SVD + cosine UMAP. The `summary.json` records both this run's parameters and the originals' under `original_umap_provenance` so the distinction is auditable.

Related pages: [UMAP and dimensionality reduction](../concepts/umap-and-dimensionality-reduction.md), [scRNA visualization and analysis](../concepts/scrna-visualization-and-analysis.md), [001 ADORA interpretation](001-adora-interpretation.md), [001 ADORA expression](001-adora-expression.md), [ARTIFACTS.md](../../lab/001_adora_expression/ARTIFACTS.md).
