# Concept Map

## Project Logic

[caffeine molecular targets](concepts/caffeine-molecular-targets.md)
is downstream of
[purine alkaloids](concepts/purine-alkaloids.md),
[xanthines](concepts/xanthines.md),
[caffeine in plants](concepts/caffeine-in-plants.md),
[methylxanthines and tea effects](concepts/methylxanthines-and-tea-effects.md),
and
[caffeine cultural history](concepts/caffeine-cultural-history.md).

[pharmacology vocabulary](concepts/pharmacology-vocabulary.md)
defines the words needed for
[caffeine molecular targets](concepts/caffeine-molecular-targets.md),
[adenosine receptors](concepts/adenosine-receptors.md),
and
[signaling to transcription](concepts/signaling-to-transcription.md).

[secondary caffeine targets](concepts/secondary-caffeine-targets.md)
defines the higher-dose target systems behind the PDE, RyR, ATM/ATR, and mTOR branches.

[cAMP signaling and ADORA cascades](concepts/camp-signaling.md)
is the receptor-by-receptor signaling schematic for A1, A2A, A2B, and A3.

[GPCR](concepts/gpcr.md)
is the receptor architecture the [adenosine receptors](concepts/adenosine-receptors.md) belong to: how a ligand outside the cell becomes a signal inside it.

[G-protein](concepts/g-protein.md)
is the molecular switch a [GPCR](concepts/gpcr.md) flips: a GDP/GTP-gated protein that turns enzymes up or down once activated.

[G-protein switching mechanics](concepts/g-protein-switching.md)
is the optional deep dive on *how* the [GPCR](concepts/gpcr.md) toggles the switch — it acts as a GEF that pries GDP loose so GTP can replace it.

[G-protein coupling](concepts/g-protein-coupling.md)
decodes the Gi/o, Gs, Golf, and Gq labels that determine which branch a [G-protein](concepts/g-protein.md) activates.

[kinase](concepts/kinase.md)
is the enzyme step that turns the [cAMP](concepts/camp-signaling.md) signal into gene regulation: cAMP activates PKA, which phosphorylates transcription factors like CREB.

[signaling to transcription](concepts/signaling-to-transcription.md)
connects caffeine exposure to
[epigenomics vocabulary](concepts/epigenomics-vocabulary.md)
through transcription factors, histone modifiers, and chromatin accessibility.

[public data landscape](concepts/public-data-landscape.md)
provides the observable data needed for
[cell type response model](concepts/cell-type-response-model.md)
and
[computational pipeline](concepts/computational-pipeline.md).

[single-cell RNA-seq measurement](concepts/single-cell-rna-seq-measurement.md)
explains what the single-cell RNA experiments actually measure: cells or nuclei, captured transcripts, sparse expression matrices, raw/normalized/scaled layers, dropout, and why ADORA expression should be interpreted as transcript evidence rather than receptor protein.

[scRNA visualization and analysis](concepts/scrna-visualization-and-analysis.md)
is the downstream counterpart to the measurement page: the standard Scanpy/Seurat/scvi-tools/Harmony pipeline, every plot type practitioners actually use (UMAP/PCA/t-SNE, dotplot, stacked violin, matrixplot, tracksplot, heatmap, violin, density UMAP, trajectories), and the recurring traps. Each technique is anchored in a Lab 001 figure where one exists.

[UMAP and dimensionality reduction](concepts/umap-and-dimensionality-reduction.md)
zooms into the embedding mechanics behind those figures: the normalize → select-genes → PCA/SVD → batch-correct → UMAP chain, what UMAP's two-phase neighbor-graph-plus-force-layout actually computes, the `n_neighbors`/`min_dist` knobs, and which features of a UMAP picture are real versus artifacts not to over-read.

[001 joint UMAP](labs/001-joint-umap.md)
is the cross-atlas method built on those fundamentals: how `make_shared_umap_tabula_hbca.py` puts Tabula Sapiens and HBCA on one set of axes (shared-gene intersection, stratified sampling, TruncatedSVD, Harmony vs centering batch correction, cosine UMAP) and how to read the resulting dataset-colored QC plot and ADORA-high overlay.

[epithelial cell types](concepts/epithelial-cell-types.md)
defines epithelial, basal, and stratified squamous cells so ADORA2B-rich tissue dots such as `Tongue` can be read as barrier/surface-epithelium signals rather than vague tissue names.

[immune cell types](concepts/immune-cell-types.md)
defines myeloid, dendritic, macrophage, monocyte, T-cell, B-cell, mast-cell, and basophil labels so ADORA3-rich immune dots can be interpreted upstream of any one lab result.

[direct caffeine epigenomics](concepts/direct-caffeine-epigenomics.md)
anchors the model in actual perturbation evidence, especially HUVEC accessibility and reporter-assay work.

[gwas ewas pharmacogenomics](concepts/gwas-ewas-pharmacogenomics.md)
connects population-level caffeine traits to genes, enhancers, methylation sites, candidate causal tissues, and exposure definitions such as coffee, tea, plasma caffeine, and synthetic caffeine.

[caffeine sensitivity genetics](concepts/caffeine-sensitivity-genetics.md)
splits sensitivity into pharmacokinetic loci such as CYP1A2/AHR and pharmacodynamic loci such as ADORA2A.

[hank caffeine video](raw/hank-caffeine-video.md)
is the source dossier that broadened the wiki beyond the original [ADORA](concepts/adenosine-receptors.md)-first frame.

[001 adora expression](labs/001-adora-expression.md)
is the first concrete implementation step: receptor expression mapping.

[001 ADORA interpretation](labs/001-adora-interpretation.md)
is the companion page for reading the current Tabula Sapiens outputs: thresholds, dotplots, sparse UMAP overlays, missing brain coverage, and the checklist for turning blobs into cell-type hypotheses.

[001 CellxGene Census API](labs/001-cellxgene-census-api.md),
[001 data flow](labs/001-data-flow.md),
[001 H5AD and AnnData cache](labs/001-h5ad-anndata-cache.md),
and
[001 notebook guide](labs/001-notebook-guide.md)
explain the API, data object, download slice, HDF5 cache files, and notebook structure behind Lab 001.

[Census core objects](concepts/census-core-objects.md)
defines the TileDB-SOMA object vocabulary behind the API: Collection, DataFrame, TableReadIter, Arrow Table, Experiment, Measurement, X array, and `soma_joinid`.

[Census experiment tree](concepts/census-experiment-tree.md)
maps the concrete human `Experiment` observed in Census `2025-11-08`: `obs`, `ms["RNA"]`, `var`, `X["raw"]`, `X["normalized"]`, and `feature_dataset_presence_matrix`.

[Census obs columns](concepts/census-obs-columns.md),
[Census var columns](concepts/census-var-columns.md),
and
[Census X layers and feature presence](concepts/census-x-layers-and-feature-presence.md)
are the column-level glossaries for the object tree. They explain labels such as assay, cell type, tissue, disease, donor, feature name, nonzero counts, raw expression, normalized expression, and feature presence.

[SOMA axes and X](concepts/soma-axes-and-x.md)
explains the logical alignment between `obs.soma_joinid`, `var.soma_joinid`, and the sparse `X` array before the discussion drops down to physical storage.

[TileDB-SOMA storage](concepts/tiledb-soma-storage.md)
explains how the Census is laid out as cell-major sparse fragments on S3, and why that layout makes Lab 001's "few genes, many cells" query the worst-case access pattern.

[Census source H5ADs](concepts/census-source-h5ads.md)
documents the *other* physical projection of the same data — one materialized H5AD per source dataset on S3. The escape hatch from the fragment-walk problem when the query shape is "narrow gene set across many datasets", which is exactly Lab 001's shape.

[network and I/O instrumentation](concepts/network-and-io-instrumentation.md)
documents the tools and signals needed to measure what a slow fetch is actually doing (bandwidth, requests, syscalls, stacks). It pairs with [TileDB-SOMA storage](concepts/tiledb-soma-storage.md) for diagnosis.

[001 fetch stall post-mortem](labs/001-fetch-stall-postmortem.md)
is the lessons-learned record of the 2026-05-31 brain fetch and the bridge from those two concept pages back to Lab 001 design decisions. Its source observations are in [raw lab-001 stall post-mortem](raw/lab-001-stall-postmortem.md).

[001 v3 stratified fetch](labs/001-v3-stratified-fetch.md)
is the phase-by-phase walkthrough of the current fetch pipeline. It applies the lessons from the post-mortem and the [pushdown](concepts/tiledb-soma-storage.md) section: stop filtering on attributes, start filtering on the obs dimension via `obs_coords`. Empirical numbers from the 2026-06-03 run.

## Supporting Sources

Provenance: which raw dossier each page is built on. Verify a claim by following the page back to its source here.

- [hank caffeine video](raw/hank-caffeine-video.md) — *vlogbrothers, "Caffeine is Very, Very Strange", 2025-09-05* — is the dossier behind [caffeine in plants](concepts/caffeine-in-plants.md), [methylxanthines and tea effects](concepts/methylxanthines-and-tea-effects.md), [caffeine cultural history](concepts/caffeine-cultural-history.md), and the plant/invertebrate sections of [caffeine molecular targets](concepts/caffeine-molecular-targets.md). Its claim map names the primary literature behind each video claim.
- [lab 001 source](raw/lab-001-source.md) backs [001 adora expression](labs/001-adora-expression.md).
- [lab 001 stall post-mortem (raw)](raw/lab-001-stall-postmortem.md) backs [001 fetch stall post-mortem](labs/001-fetch-stall-postmortem.md) and grounds the worked example in [TileDB-SOMA storage](concepts/tiledb-soma-storage.md) and [network and I/O instrumentation](concepts/network-and-io-instrumentation.md).
- [proposal source](raw/proposal-source.md) and [readme source](raw/readme-source.md) are the project-level framing behind the receptor-and-epigenomics pages: [adenosine receptors](concepts/adenosine-receptors.md), [signaling to transcription](concepts/signaling-to-transcription.md), [epigenomics vocabulary](concepts/epigenomics-vocabulary.md), [public data landscape](concepts/public-data-landscape.md), [direct caffeine epigenomics](concepts/direct-caffeine-epigenomics.md), [gwas ewas pharmacogenomics](concepts/gwas-ewas-pharmacogenomics.md), [caffeine sensitivity genetics](concepts/caffeine-sensitivity-genetics.md), [cell type response model](concepts/cell-type-response-model.md), [computational pipeline](concepts/computational-pipeline.md), and [research questions](concepts/research-questions.md).
- [karpathy wiki pattern](raw/karpathy-wiki-pattern.md) documents the wiki style itself, not caffeine biology.

## Dependency Order

1. Read [hank caffeine video](raw/hank-caffeine-video.md).
2. Learn [purine alkaloids](concepts/purine-alkaloids.md).
3. Learn [xanthines](concepts/xanthines.md).
4. Learn [caffeine in plants](concepts/caffeine-in-plants.md).
5. Learn [caffeine cultural history](concepts/caffeine-cultural-history.md).
6. Learn [methylxanthines and tea effects](concepts/methylxanthines-and-tea-effects.md).
7. Learn [pharmacology vocabulary](concepts/pharmacology-vocabulary.md).
8. Learn [caffeine molecular targets](concepts/caffeine-molecular-targets.md).
9. Learn [adenosine receptors](concepts/adenosine-receptors.md).
10. Learn [GPCR](concepts/gpcr.md).
11. Learn [G-protein](concepts/g-protein.md).
12. (Optional deep dive) Learn [G-protein switching mechanics](concepts/g-protein-switching.md).
13. Learn [G-protein coupling](concepts/g-protein-coupling.md).
14. Learn [cAMP signaling and ADORA cascades](concepts/camp-signaling.md).
15. Learn [kinase](concepts/kinase.md).
16. Learn [secondary caffeine targets](concepts/secondary-caffeine-targets.md).
17. Learn [signaling to transcription](concepts/signaling-to-transcription.md).
18. Learn [epigenomics vocabulary](concepts/epigenomics-vocabulary.md).
19. Learn [public data landscape](concepts/public-data-landscape.md).
20. Learn [single-cell RNA-seq measurement](concepts/single-cell-rna-seq-measurement.md).
21. Learn [scRNA visualization and analysis](concepts/scrna-visualization-and-analysis.md).
21a. Learn [UMAP and dimensionality reduction](concepts/umap-and-dimensionality-reduction.md), then [001 joint UMAP](labs/001-joint-umap.md).
22. Learn [epithelial cell types](concepts/epithelial-cell-types.md).
23. Learn [immune cell types](concepts/immune-cell-types.md).
24. Read [direct caffeine epigenomics](concepts/direct-caffeine-epigenomics.md).
25. Read [cell type response model](concepts/cell-type-response-model.md).
26. Run or inspect [001 adora expression](labs/001-adora-expression.md).
27. Read [001 ADORA interpretation](labs/001-adora-interpretation.md).
28. Use [research questions](concepts/research-questions.md) to pick the next lab.
