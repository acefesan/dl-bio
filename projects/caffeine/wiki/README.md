# Caffeine Epigenome Wiki

This is a Karpathy-style LLM wiki for the caffeine epigenome project: persistent Markdown pages that compile project knowledge once, cross-link the concepts, and give future humans or agents a stable place to continue from.

Start here:

- [reading itinerary](reading-itinerary.md): a single front-to-back data-and-methods read in one sitting — download → UMAP fundamentals → the joint UMAP → every figure → open questions.
- [learning path](learning-path.md): the ramp-up order if you are new to the project.
- [concept map](concept-map.md): the graph of concepts and how they depend on each other.
- [research questions](concepts/research-questions.md): the project backlog, ordered by difficulty.
- [computational pipeline](concepts/computational-pipeline.md): the proposed end-to-end analysis pipeline.
- [maintenance](maintenance.md): rules for extending the wiki without turning it into notes soup.

## Core Thesis

Caffeine has sparse direct epigenomic perturbation data, so the tractable path is computational imputation: combine baseline cell atlases, receptor expression, chromatin accessibility, transcription factor activity, direct caffeine datasets where available, and caffeine GWAS/EWAS signals to predict cell-type-specific regulatory responses.

The broader frame is that caffeine is not only an adenosine-receptor ligand in humans. It is also a convergently evolved plant [purine alkaloid](concepts/purine-alkaloids.md), an ecological defense molecule, a pollinator-behavior molecule, a family member of related [methylxanthines](concepts/xanthines.md), and a culturally/industrially diverse exposure.

## Page Index

### Raw Sources

- [readme source](raw/readme-source.md): project README summary.
- [proposal source](raw/proposal-source.md): full proposal summary.
- [lab 001 source](raw/lab-001-source.md): [ADORA](concepts/adenosine-receptors.md) expression lab summary.
- [lab 001 stall post-mortem (raw)](raw/lab-001-stall-postmortem.md): observed log and process state of the 2026-05-31 brain fetch.
- [hank caffeine video](raw/hank-caffeine-video.md): claim map and source dossier for the vlogbrothers caffeine video.
- [karpathy wiki pattern](raw/karpathy-wiki-pattern.md): summary of the wiki style used here.

### Concepts

- [pharmacology vocabulary](concepts/pharmacology-vocabulary.md): plain-language glossary of the basic terms (ligand, receptor, agonist, antagonist) used everywhere else — start here if you are not a biologist.
- [caffeine molecular targets](concepts/caffeine-molecular-targets.md)
- [purine alkaloids](concepts/purine-alkaloids.md)
- [xanthines](concepts/xanthines.md)
- [caffeine in plants](concepts/caffeine-in-plants.md)
- [secondary caffeine targets](concepts/secondary-caffeine-targets.md)
- [cAMP signaling and ADORA cascades](concepts/camp-signaling.md)
- [methylxanthines and tea effects](concepts/methylxanthines-and-tea-effects.md)
- [caffeine cultural history](concepts/caffeine-cultural-history.md)
- [adenosine receptors](concepts/adenosine-receptors.md)
- [GPCR](concepts/gpcr.md)
- [G-protein](concepts/g-protein.md)
- [G-protein switching mechanics](concepts/g-protein-switching.md)
- [G-protein coupling](concepts/g-protein-coupling.md)
- [kinase](concepts/kinase.md)
- [signaling to transcription](concepts/signaling-to-transcription.md)
- [epigenomics vocabulary](concepts/epigenomics-vocabulary.md)
- [public data landscape](concepts/public-data-landscape.md)
- [single-cell RNA-seq measurement](concepts/single-cell-rna-seq-measurement.md): what scRNA-seq experiments measure, how an RNA matrix is produced, and how to interpret sparse receptor values.
- [scRNA visualization and analysis](concepts/scrna-visualization-and-analysis.md): the standard tools (Scanpy, Seurat, scvi-tools, Harmony), the plot landscape (UMAP, dotplot, stacked violin, matrixplot, tracksplot, heatmap, density UMAP), downstream analyses (clustering, marker finding, batch correction, trajectories, cell-cell communication, perturbation), and the recurring traps.
- [epithelial cell types](concepts/epithelial-cell-types.md): plain-language guide to epithelial, basal, and stratified squamous cells, with the tongue ADORA2B example.
- [immune cell types](concepts/immune-cell-types.md): plain-language guide to myeloid, dendritic, macrophage, monocyte, T-cell, B-cell, mast-cell, and basophil labels.
- [direct caffeine epigenomics](concepts/direct-caffeine-epigenomics.md)
- [gwas ewas pharmacogenomics](concepts/gwas-ewas-pharmacogenomics.md)
- [caffeine sensitivity genetics](concepts/caffeine-sensitivity-genetics.md)
- [cell type response model](concepts/cell-type-response-model.md)
- [unexpected responsive cell types](concepts/unexpected-responsive-cell-types.md)
- [computational pipeline](concepts/computational-pipeline.md)
- [research questions](concepts/research-questions.md)
- [Census core objects](concepts/census-core-objects.md): the object vocabulary behind CellxGene Census: Collection, DataFrame, TableReadIter, Arrow table, Experiment, Measurement, X array, and `soma_joinid`.
- [Census experiment tree](concepts/census-experiment-tree.md): the human `Experiment` object tree: `obs`, `ms["RNA"]`, `var`, `X`, and `feature_dataset_presence_matrix`.
- [Census obs columns](concepts/census-obs-columns.md): glossary for human cell metadata columns such as assay, cell type, tissue, disease, donor, and QC summaries.
- [Census var columns](concepts/census-var-columns.md): glossary for RNA feature metadata columns such as `feature_id`, `feature_name`, `nnz`, and `n_measured_obs`.
- [Census X layers and feature presence](concepts/census-x-layers-and-feature-presence.md): explains `X["raw"]`, `X["normalized"]`, sparse triples, and the feature presence matrix.
- [SOMA axes and X](concepts/soma-axes-and-x.md): how `obs.soma_joinid` and `var.soma_joinid` define the cell and gene coordinates of the sparse expression array.
- [TileDB-SOMA storage](concepts/tiledb-soma-storage.md): what the Census actually stores on S3 and why "few genes, all cells" is the worst-case query.
- [Census source H5ADs](concepts/census-source-h5ads.md): the *other* physical projection of the same data — one materialized H5AD per source dataset, no fragment walking.
- [network and I/O instrumentation](concepts/network-and-io-instrumentation.md): tools and signals for measuring what a slow fetch is really doing.

### Labs

- [001 adora expression](labs/001-adora-expression.md)
- [001 ADORA interpretation](labs/001-adora-interpretation.md): how to read the local ADORA UMAPs, dotplots, thresholds, and Tabula Sapiens coverage limits.
- [001 CellxGene Census API](labs/001-cellxgene-census-api.md)
- [001 data flow](labs/001-data-flow.md)
- [001 H5AD and AnnData cache](labs/001-h5ad-anndata-cache.md)
- [001 notebook guide](labs/001-notebook-guide.md)
- [001 fetch stall post-mortem](labs/001-fetch-stall-postmortem.md): why the 2026-05-31 brain fetch was unworkable and what to change next.
- [001 v3 stratified fetch](labs/001-v3-stratified-fetch.md): phase-by-phase walkthrough of the current fetch pipeline (obs scan → stratified sample → coord-based X read → atomic write) with empirical numbers from the 2026-06-03 run.

## Current Project State

The project is at proposal plus first lab. The active lab is Q1: map expression of the [ADORA](concepts/adenosine-receptors.md) adenosine receptor genes (ADORA1, ADORA2A, ADORA2B, and ADORA3) across human cell types using CellxGene Census, with GTEx v8 as a bulk-tissue sanity check.
