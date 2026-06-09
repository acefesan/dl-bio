# Single-Cell RNA-Seq Visualization And Analysis

## Summary

Single-cell RNA-seq starts with a sparse cell × gene matrix and ends with a story about cell populations and genes. The plot types and analysis steps between those two endpoints are mostly the same across labs because the community has converged on a small set of tools:

- **Scanpy** (Python, AnnData-backed) and **Seurat** (R) are the two foundational frameworks. For million-cell datasets like ours, Scanpy is the practical default.
- **scvi-tools** for probabilistic models (scVI, totalVI), **Harmony** for batch correction, **CellBender** for ambient RNA removal — these are the standard add-ons.
- Almost every plot you see in a single-cell paper comes out of one of ~7 Scanpy/Seurat plotting functions.

This page is a guided tour of the technique landscape with pointers to the canonical Scanpy function for each one, the bias each plot tends to introduce, and the Lab 001 figure that illustrates the technique where one exists.

For what the measurement itself actually is, see [single-cell RNA-seq measurement](single-cell-rna-seq-measurement.md). For how this lab interprets its current figures, see [001 ADORA interpretation](../labs/001-adora-interpretation.md).

## The Standard Analysis Pipeline

Practitioner consensus on the conventional steps, in order:

| Step | Standard tool | What it produces |
|---|---|---|
| 1. **Quality control** | `scanpy.pp.calculate_qc_metrics`, scrublet | per-cell filters: min genes, max mito %, doublet calls |
| 2. **Normalization** | `sc.pp.normalize_total` + `sc.pp.log1p` | library-size adjusted, log-transformed counts |
| 3. **Highly variable gene (HVG) selection** | `sc.pp.highly_variable_genes` | shortlist of 2k–4k genes that vary across cells |
| 4. **Dimensionality reduction (linear)** | `sc.tl.pca` | top 50 PCs as a denoised representation |
| 5. **Batch correction / integration** | Harmony, scVI, scanorama, BBKNN | embeddings aligned across donors/assays |
| 6. **Neighbor graph** | `sc.pp.neighbors` | kNN graph in PCA / integrated space |
| 7. **Clustering** | `sc.tl.leiden` (Scanpy) or Louvain (Seurat) | per-cell cluster label |
| 8. **2D embedding for visualization** | `sc.tl.umap`, sometimes `tsne` | 2D coords for plotting only |
| 9. **Marker / DE gene finding** | `sc.tl.rank_genes_groups` | genes that distinguish each cluster |
| 10. **Cell-type annotation** | CellTypist, manual, ontology lookup | biological labels for clusters |
| 11. **Downstream**: trajectories, communication, perturbation | scVelo, CellRank, CellChat, CellOracle | dynamics and inter-cell signaling |

Tabula Sapiens delivers the file already past step 10 — `obs.cell_type` is annotated, multiple embeddings are precomputed in `obsm/X_pca`, `obsm/X_scvi`, `obsm/X_umap`, and the neighbor graph and Leiden labels are in `obsp/` and `uns/`. That's why Lab 001 can drop straight into "plot ADORA expression" without redoing the upstream pipeline.

## Visualization Techniques

These are the plot types people actually make. Each one has a specific question it answers well and a way it can mislead.

### Embedding scatterplots (UMAP, PCA, t-SNE)

```python
sc.pl.embedding(adata, basis="X_umap", color="ADORA1")
sc.pl.umap(adata, color=["ADORA1", "ADORA2A", "cell_type"])
```

Project each cell from 60k-gene space down to 2D so you can color every cell by something and see whether the colored cells cluster.

**Practitioner consensus on the three flavors:**

| Embedding | When | Quirks |
|---|---|---|
| **PCA** | sanity check; the first few PCs often capture the dominant biology | linear, no neighborhood preservation |
| **UMAP** | default for visualization | distances between distant clusters are not meaningful; cluster shape is not biology |
| **t-SNE** | alternative to UMAP, sometimes nicer for highly structured data | slow at large N; tends to over-tighten clusters |

UMAP is the genre default. Seurat uses it. Scanpy uses it. Tabula Sapiens publishes it.

**Lab 001 figure that uses this:** `figures/tabula_sapiens_embeddings_by_tissue.png` shows the same 1.1M cells projected six different ways (PCA, scVI, UMAP, scVI-corrected UMAP, uncorrected UMAP variants) — useful for seeing how the embedding choice changes what's visible.

**The trap:** coloring 1.1M cells by a sparse gene like `ADORA1` (only 8,416 cells positive out of 1.14M) gives you a near-white plot where the 99.3% zero cells overdraw the 0.7% positive ones. The standard fix is to render zero-cells faintly and only draw nonzero or above-threshold cells brightly — exactly what `figures/tabula_sapiens_adora_high_umap_four_color.png` does.

### Dotplot

```python
sc.pl.dotplot(adata, var_names=["ADORA1","ADORA2A","ADORA2B","ADORA3"],
              groupby="cell_type")
```

**The single most useful plot in single-cell biology.** Each dot answers one `(cell_group × gene)` question:

| Channel | Meaning |
|---|---|
| dot size | percent of cells in that group with expression > 0 |
| dot color | mean expression in that group (often only among expressing cells) |
| rows | cell types or tissues |
| columns | genes |

The reason it dominates: it separates *prevalence* (how many cells in this group express the gene) from *magnitude* (when they do express it, how much). Both matter, and they often disagree.

**Lab 001 figures that use this:** `figures/tabula_sapiens_adora_dotplot_cell_type.png` (35 cell types × 4 receptors) and `figures/tabula_sapiens_adora_dotplot_tissue.png` (28 tissues × 4 receptors). Reading them: ADORA3 in microglia shows a big bright dot — common and strong. ADORA1 in pancreatic ductal cells shows a smaller bright dot — fewer cells, strong when present. ADORA2B in trachea/bladder shows large dots with moderate color — broad weak signal across the tissue.

### Stacked violin

```python
sc.pl.stacked_violin(adata, var_names=["ADORA1","ADORA2A","ADORA2B","ADORA3"],
                     groupby="cell_type")
```

Stacks compact violin plots — one per (gene × group) — so you can see the full per-cell distribution, not just the mean. Useful when "mean expression" hides a bimodal pattern (some cells very high, most cells zero). Dotplot summarizes; stacked violin shows the shape.

### Matrixplot

```python
sc.pl.matrixplot(adata, var_names=[...], groupby="cell_type")
```

A heatmap of mean expression by group. Same info as the *color* of a dotplot, without the size channel. Use when you don't care about prevalence — e.g., comparing magnitudes across already-known-expressing populations.

### Tracksplot

```python
sc.pl.tracksplot(adata, var_names=[...], groupby="cell_type")
```

Same info as a heatmap but height-encoded instead of color-encoded. Some readers find height easier to compare than color saturation. Niche choice.

### Heatmap

```python
sc.pl.heatmap(adata, var_names=[...], groupby="cell_type")
```

Every cell is a row; every gene is a column; cells grouped by category. Useful for showing within-group heterogeneity. With 1.1M cells, downsample or summarize first or it's unreadable.

### Violin / ridge plots

```python
sc.pl.violin(adata, keys=["ADORA1"], groupby="cell_type")
```

One full distribution per group, expanded. Cleanest single-gene comparison across cell types. Loses compactness when you have >20 cell types.

### Tissue / per-tissue breakdowns

Take one tissue, regroup cells by that tissue's cell types, dotplot it. This is how you go from "ADORA2B looks high in tongue" to "ADORA2B is mainly in tongue basal cells, not taste cells." Lab 001's `figures/tabula_sapiens_tongue_adora_cell_type_breakdown.png` is the textbook version of this drill-down.

### Less common, sometimes the right answer

| Plot | When |
|---|---|
| **Bubble plot of regulons** (SCENIC) | TF activity per cell type instead of gene expression |
| **Sankey / alluvial** | Showing how cells move between annotation schemes |
| **Density UMAP** (`sc.pl.embedding_density`) | Where is condition A vs B concentrated on the UMAP |
| **Trajectory plots** (PAGA, scVelo streamlines) | Continuous transitions instead of discrete clusters |
| **Cell-cell communication chord diagrams** (CellChat) | Inferred ligand→receptor signaling between cell types |

## Analysis Techniques Beyond Plotting

### Highly variable gene selection

Of 60k genes, most don't vary across cells in informative ways. HVG selection picks the ~2k–4k that do. Standard implementation in both Scanpy (`sc.pp.highly_variable_genes`) and Seurat: bin genes by mean expression, take the top variance-to-mean ratio per bin. This is what feeds PCA and ultimately UMAP.

### Clustering

The Leiden algorithm (in Scanpy) and Louvain (in Seurat) are graph-based community detection on the kNN cell graph. Both have a `resolution` parameter. Higher resolution → more, smaller clusters. Practitioners universally agree: try multiple resolutions, then merge or relabel. Picking one number and committing is a common rookie move.

### Batch correction and integration

Different donors, assays, labs, or tissue preparations produce technical variance that looks like biological variance in the embedding. Three families of methods are widely used:

| Method | Approach | When |
|---|---|---|
| **Harmony** | iteratively corrects PCA coordinates | fast, often the first thing to try |
| **scVI** (scvi-tools) | trains a variational autoencoder, gives a corrected latent space | when you need a probabilistic model or to handle many batches |
| **scanorama / BBKNN** | match nearest neighbors across batches | lightweight, sometimes preferred for integration without a learned model |

Tabula Sapiens used scVI to produce the `X_scvi` embedding visible in its multi-embedding plot. The "uncorrected" UMAPs in that same figure are what the data looks like *without* batch correction — clusters are smeared by donor and assay rather than separated by cell type.

### Differential expression / marker finding

`sc.tl.rank_genes_groups(adata, "cell_type")` ranks genes that distinguish each group. Backends: t-test, Wilcoxon (most common), logistic regression. The output is a per-cluster ranked gene list, used both for cell-type annotation and downstream interpretation. Subtle but important: significance with millions of cells is almost meaningless — effect sizes matter more than p-values.

### Cell-type annotation

Three styles, often combined:

1. **Manual** — look at marker gene rankings per cluster, name the cluster.
2. **Reference-based** — `CellTypist`, `Azimuth`, etc., transfer labels from a reference atlas. Fast, but blind to anything the reference didn't see.
3. **Ontology lookup** — map labels to the Cell Ontology (CL) and Uberon terms. This is what gives Tabula Sapiens cells their `cell_type_ontology_term_id` column.

### Trajectory / pseudotime inference

`scVelo`, `CellRank`, `PAGA`, `Slingshot` — methods that order cells along inferred developmental or state-transition trajectories. Useful when you suspect a continuum (stem → progenitor → differentiated) instead of discrete types. Caffeine biology doesn't obviously need this for Q1, but later questions about signaling cascades might.

### Cell-cell communication

`CellChat`, `LIANA`, `NicheNet` — infer ligand-receptor interactions between cell types from gene expression. Naturally relevant when you care about whether ADORA-expressing cells receive caffeine-relevant signaling input from neighbors.

### Perturbation simulation

`CellOracle` builds gene regulatory networks and lets you simulate "what would expression look like if I knocked down gene X?" The proposal calls for this as a way to predict caffeine response in cell types where no perturbation experiment was done.

## What Bioinformaticians Actually Reach For First

Distilled from practitioner discussions and the standard tutorials:

1. **First plot of any new dataset:** UMAP colored by `cell_type` and by `batch / donor`. Sanity check.
2. **First plot when you have a target gene set:** dotplot, cell_type on rows, genes on columns. This is what tells you where the signal is.
3. **First plot when you want to compare conditions:** stacked violin or split violin by condition.
4. **When the dotplot looks confusing:** drill down with a per-tissue or per-broader-cell-class dotplot. The tongue → basal cell example in this lab is the canonical pattern.
5. **When clusters look weird:** color the UMAP by `donor`, `assay`, `n_genes_by_counts`, `pct_counts_mt`. The weird cluster is usually one batch or one quality-control outlier.

## The Recurring Traps

| Trap | What goes wrong | Mitigation |
|---|---|---|
| Reading UMAP geometry as biology | Big distance between clusters doesn't mean big biological difference. Cluster shape isn't meaningful. | Confirm with dotplots, marker gene lists, and quantitative tests |
| Sparse-gene UMAP coloring | Zero-cells overdraw rare positive cells | Filter to nonzero, threshold, or render zeros faintly |
| Picking one clustering resolution | The "right" cluster count is task-dependent | Try multiple, justify the choice |
| Skipping batch correction | Donor or assay effects masquerade as cell types | Always plot UMAP colored by batch before annotation |
| Believing percent-expressing without checking feature presence | A gene's `0` can mean "not measured in this source dataset" | Use the [feature presence matrix](census-x-layers-and-feature-presence.md) for Census-derived data |
| Inferring protein from RNA | Single-cell RNA is evidence of transcript, not receptor on the membrane | State the limit; flag candidates for orthogonal validation |
| Significance with N=millions | Tiny effects become "highly significant" | Report effect sizes, not just p-values |

## What Lab 001 Has Already Done

Mapping the techniques above onto the current figure inventory:

| Figure | Technique | Notes |
|---|---|---|
| `tabula_sapiens_embeddings_by_tissue.png` | embedding scatter (six variants) | shows why scVI-corrected UMAP separates cell types better than uncorrected |
| `tabula_sapiens_embeddings_by_tissue_fine.png` | embedding scatter, finer label granularity | same technique, finer tissue labels |
| `tabula_sapiens_adora_high_umap.png` | thresholded UMAP coloring | renders only above-threshold cells to escape the zero-overdraw trap |
| `tabula_sapiens_adora_high_umap_four_color.png` | four-color thresholded UMAP | same idea, one color per receptor — shows the four spatial niches |
| `tabula_sapiens_adora_dotplot_cell_type.png` | dotplot (cell type × gene) | the main inference tool for Q1 |
| `tabula_sapiens_adora_dotplot_tissue.png` | dotplot (tissue × gene) | broad anatomical direction |
| `tabula_sapiens_tongue_adora_cell_type_breakdown.png` | per-tissue dotplot drilldown | textbook follow-up: tissue → cell-type composition |

What's not yet plotted (and might be useful for Q1):

- Stacked violin of ADORA per cell type — shows the within-group distribution
- Donor-stratified dotplot — does any single donor drive a signal
- Assay-stratified dotplot — does the receptor look very different between 10x and Smart-seq2
- `sc.pl.embedding_density` per receptor — where is the signal concentrated independent of cell-type labels
- Marker gene re-derivation (`sc.tl.rank_genes_groups`) on ADORA-positive vs ADORA-negative cells — what *else* characterizes the positive populations
- Receptor co-expression bar chart (ADORA1-only, ADORA2A-only, …, multi-receptor) — which cells express more than one

## Sources

- [scanpy core plotting functions](https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html) — canonical reference for `dotplot`, `matrixplot`, `stacked_violin`, `tracksplot`, `heatmap`.
- [Current best practices in single-cell RNA-seq analysis: a tutorial (Luecken & Theis 2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6582955/) — the standard reference workflow.
- [The impact of package selection and versioning on single-cell RNA-seq analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC11014608/) — Seurat vs Scanpy in practice.
- [Practical bioinformatics pipelines for single-cell RNA-seq data analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10189648/) — workflow-stage reference.
- [Top bioinformatics tools for scRNA-seq in 2025](https://neovarsity.org/blogs/top-bioinformatics-tools-for-scrna-seq) — current tool landscape (Scanpy, Seurat, scvi-tools, Harmony, CellBender, Squidpy).
- [Over 1000 tools reveal trends in the single-cell RNA-seq analysis landscape](https://www.biorxiv.org/content/10.1101/2021.08.13.456196.full.pdf) — long-tail tool ecosystem.

Related pages: [single-cell RNA-seq measurement](single-cell-rna-seq-measurement.md), [001 ADORA interpretation](../labs/001-adora-interpretation.md), [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [Census obs columns](census-obs-columns.md), [Census var columns](census-var-columns.md), [001 ADORA expression](../labs/001-adora-expression.md), [computational pipeline](computational-pipeline.md)
