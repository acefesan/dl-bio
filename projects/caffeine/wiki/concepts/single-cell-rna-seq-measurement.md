# Single-Cell RNA-Seq Measurement

## Summary

Single-cell RNA-seq measures **RNA molecules captured from individual cells or nuclei**. For this project, the key output is a sparse cell x gene matrix:

```text
cell barcode x gene -> expression value
```

An ADORA value in that matrix is not a direct receptor-protein measurement. It is evidence that RNA from `ADORA1`, `ADORA2A`, `ADORA2B`, or `ADORA3` was captured and assigned to a cell. That makes it useful for finding candidate caffeine-responsive cell populations, but it must be interpreted with single-cell sparsity, tissue coverage, and assay chemistry in mind.

For where the matrix lives in Census and H5AD files, see [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [SOMA axes and X](soma-axes-and-x.md), and [Census source H5ADs](census-source-h5ads.md).

## What The Experiment Does

Most single-cell RNA-seq pipelines follow the same conceptual path:

1. **Collect tissue.** A donor tissue sample is removed, preserved, and processed.
2. **Dissociate tissue.** The sample is broken into a suspension of cells or nuclei. This step can enrich or lose cell types.
3. **Capture cells.** Droplet methods such as 10x Genomics isolate many cells with barcoded beads. Plate methods such as Smart-seq2 sort cells into wells.
4. **Capture RNA.** Cellular RNA is reverse-transcribed into cDNA and tagged with cell barcodes and molecular identifiers when the chemistry supports them.
5. **Sequence reads.** The cDNA library is sequenced.
6. **Map reads to genes.** Reads are aligned or pseudoaligned to a reference transcriptome/genome.
7. **Count molecules.** Reads are collapsed into a gene-count matrix.
8. **Annotate cells.** Cells are assigned labels such as tissue, cell type, donor, assay, and dataset.

The Tabula Sapiens publication describes a cross-tissue atlas built with both 10x droplet sequencing and plate-based Smart-seq2, with expert cell-type annotation regularized to public ontologies. The local Lab 001 file is a processed H5AD projection of that kind of output, not raw FASTQ reads.

## Droplet Versus Plate Measurements

| Assay style | Strength | Trade-off |
|---|---|---|
| 10x droplet scRNA-seq | Many cells, good population coverage | Shallower per cell; sparse detection for low-expression genes |
| Smart-seq2 / plate scRNA-seq | Deeper transcript coverage per sorted cell | Fewer cells; sorting choices can bias cell populations |
| single-nucleus RNA-seq | Works for frozen or hard-to-dissociate tissues such as brain | Nuclear RNA differs from whole-cell RNA; some cytoplasmic transcripts are underrepresented |

For ADORA, this matters because GPCR transcripts can be low and sparse. A receptor may be biologically important even if only a small fraction of cells have captured mRNA in a single-cell matrix.

## What `Expression > 0` Means

In a sparse single-cell matrix, zero and nonzero have asymmetric meanings:

| Matrix value | Safer interpretation |
|---|---|
| `> 0` | The assay captured reads/molecules assigned to this gene in this cell. |
| `0` with gene present in the source dataset | No expression was observed for this cell. This may be true zero or technical dropout. |
| `0` with gene absent from the source dataset | Not interpretable as biological absence. |

For Census-derived work, the [feature presence matrix](census-x-layers-and-feature-presence.md) is the guardrail against calling a gene absent when the source dataset did not measure it.

## Raw, Normalized, And Scaled Values

Different layers answer different questions:

| Layer/value | Meaning | Good use |
|---|---|---|
| Raw counts | Molecule/read counts after pipeline processing | Detection fraction; "was this gene observed?" |
| Normalized expression | Counts adjusted for library size | Comparing rough expression magnitude across cells or groups |
| Scaled expression | Centered/scaled values for modeling or visualization | Embeddings, clustering, model features; not a direct abundance unit |

The local Tabula Sapiens H5AD has `X`, `raw/X`, `layers/decontXcounts`, and `layers/scale_data`. Lab 001's current ADORA extraction uses `X` from the cached H5AD. The interpretation should therefore be framed as "expression values in this processed atlas layer", with raw and decontaminated layers available for follow-up checks.

## Why Low Expression Looks Patchy

Single-cell receptor plots often look like sparse constellations rather than smooth tissue fields. Several things create that pattern:

- **Biology:** only some cell types or cell states express the receptor.
- **Bursting:** transcription happens in pulses, so cells of the same type can differ at capture time.
- **Depth:** a shallow cell may miss low-abundance transcripts.
- **Dropout:** mRNA was present but not captured or sequenced.
- **Ambient RNA:** free RNA in the suspension can contaminate droplets.
- **Dissociation stress:** processing can induce or suppress genes.
- **Batch and donor composition:** one tissue/donor may dominate a visual cluster.

That is why a UMAP blob is a clue, not a conclusion. The next move is to summarize by cell type, tissue, donor, assay, and source dataset.

## What This Means For ADORA

For caffeine biology, the single-cell RNA matrix answers the first filter:

```text
Which cells have evidence of ADORA receptor transcript?
```

It does not by itself answer:

- Is receptor protein present on the membrane?
- Is the receptor functional?
- How much caffeine reaches this cell type in vivo?
- What happens downstream after receptor blockade?
- Is the ADORA promoter or enhancer accessible in this cell type?

The Lab 001 expression results should seed follow-up analyses, especially [cell type response model](cell-type-response-model.md), [cAMP signaling and ADORA cascades](camp-signaling.md), and future chromatin accessibility work at ADORA loci.

## Reading Plots

| Plot | What it is good for | Main trap |
|---|---|---|
| UMAP colored by expression | Finding spatial neighborhoods with receptor signal | UMAP geometry is a projection; visual blobs are not proof of biology |
| Dotplot by cell type | Separating prevalence from magnitude | Cell-type labels can hide tissue/donor mixtures |
| Tissue dotplot | Broad anatomical direction | Missing tissues are absent from the atlas, not necessarily receptor-negative |
| Interactive UMAP | Connecting expression to cell metadata by hover/click | Sparse genes need thresholded/faint zero layers |

For the full landscape of plot types (UMAP/PCA/t-SNE, dotplot, stacked violin, matrixplot, tracksplot, heatmap, violin, density UMAP), the downstream analyses (clustering, marker finding, batch correction, trajectories, cell-cell communication, perturbation simulation), and the recurring traps practitioners flag, see [scRNA visualization and analysis](scrna-visualization-and-analysis.md).

For the local Tabula Sapiens file, broad tissue labels do **not** include brain. Brain ADORA biology needs a dedicated brain atlas rather than extrapolation from this H5AD.

## Sources

- Tabula Sapiens Consortium, "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans" — describes the multi-organ atlas, 10x and Smart-seq2 measurements, and expert ontology-based annotations.
- CELLxGENE Census documentation on normalized layers — describes the Census library-size normalized expression layer.
- Local Lab 001 processing notes: `../../lab/001_adora_expression/DATASET_PROCESSING.md`.

Related pages: [scRNA visualization and analysis](scrna-visualization-and-analysis.md), [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [Census source H5ADs](census-source-h5ads.md), [Census obs columns](census-obs-columns.md), [Census var columns](census-var-columns.md), [001 ADORA interpretation](../labs/001-adora-interpretation.md), [001 ADORA expression](../labs/001-adora-expression.md)
