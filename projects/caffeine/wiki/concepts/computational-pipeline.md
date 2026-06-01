# Computational Pipeline

## Summary

The proposed pipeline turns public atlases and direct caffeine datasets into cell-type-specific caffeine response predictions. It can be orchestrated with Snakemake or Nextflow once the analysis leaves notebook exploration.

## Phase 1: Data Collection

Collect:

- CellxGene/HCA single-cell RNA and ATAC data.
- ENCODE and Roadmap epigenomic tracks.
- GTEx expression and eQTLs.
- GEO caffeine treatment datasets.
- GWAS/EWAS summary statistics.

## Phase 2: Expression Mapping

Use Scanpy, scvi-tools, CellTypist, and related tools to build cell-type x gene expression matrices for [ADORA](adenosine-receptors.md) genes, CYP1A2, PDE genes, CREB1, RYR genes, and downstream signaling genes.

## Phase 3: Chromatin Accessibility

Use ArchR, SnapATAC2, Signac, MACS3, Cicero, and chromVAR to map accessible [ADORA](adenosine-receptors.md) regulatory elements, motif activity, and gene activity scores.

## Phase 4: Epigenetic Marks

Use ChIP-seq, ChromHMM, methylation tools, and super-enhancer calling to contextualize caffeine-relevant loci.

## Phase 5: Motifs and Footprints

Use HOMER, FIMO, MEME Suite, TOBIAS, and chromVAR to connect accessible DNA to TFs such as CREB, AP-1, NF-kB, Nrf2, AHR, HNF4A, MEF2, and NFAT.

## Phase 6: GRN Inference

Use pySCENIC, SCENIC+, CellOracle, or LINGER to infer regulatory networks and simulate perturbations.

## Phase 7: Integration and Visualization

Build unified matrices and visualizations:

- heatmaps,
- genome tracks,
- network graphs,
- enrichment summaries,
- browser track hubs.

## First Concrete Pipeline Slice

The first slice is:

1. Run [001 adora expression](../labs/001-adora-expression.md).
2. Rank cell types by [ADORA](adenosine-receptors.md) receptor expression.
3. Pull matched scATAC data for the top cell types.
4. Test whether [ADORA](adenosine-receptors.md) regulatory accessibility predicts receptor expression.

Related pages: [public data landscape](public-data-landscape.md), [cell type response model](cell-type-response-model.md), [epigenomics vocabulary](epigenomics-vocabulary.md)
