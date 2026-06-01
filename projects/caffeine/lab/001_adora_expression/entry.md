# 001 — ADORA receptor expression across human cell types

**Date:** 2026-04-22
**Question:** Q1 from [PROPOSAL.md](../../PROPOSAL.md) — which human cell types express ADORA1/2A/2B/3 at highest levels?
**Status:** in progress

## Hypothesis

Adenosine receptors are expressed non-uniformly across cell types. Based on bulk-tissue data summarized in the proposal, we expect:

- **ADORA1:** brain (cortex, hippocampus, cerebellum), adipose, heart atria, kidney, testis
- **ADORA2A:** striatal medium spiny neurons, T cells, macrophages, NK cells, endothelial cells, platelets
- **ADORA2B:** intestinal epithelium, endothelial, cardiac fibroblasts, astrocytes, dendritic cells
- **ADORA3:** mast cells, neutrophils, macrophages, eosinophils, lung

Q1 sharpens this to **cell-type resolution** — which specific T-cell subtype, which striatal neuron class, which macrophage polarization state. Bulk RNA-seq can't see this.

## Setup

**Datasets**
- Primary: CellxGene Census (unified HCA + Tabula Sapiens + others, SOMA-backed)
- Secondary (validation): GTEx v8 bulk TPM for the same tissues

**Tools**
- `cellxgene_census` Python API
- `scanpy` for AnnData manipulation + plotting
- `pandas` + `matplotlib`

**Notebook:** [explore_adora_expression.ipynb](explore_adora_expression.ipynb) — pedagogical walkthrough of the CellxGene Census API. Run this first.

## Results

*To be filled in as the notebook runs.*

Planned outputs:
- `figures/dotplot_adora_cell_type.png` — cell type × 4 receptors matrix
- `figures/ranked_top20_per_receptor.png` — top 20 cell types per receptor
- `pseudobulk_by_cell_type.feather` — mean expression + % expressing per cell type
- `cross_receptor_overlap.feather` — cell types with high expression of multiple receptors

## Interpretation

*To be filled in.*

## Next steps

- Compare against GTEx bulk pseudobulk to confirm direction
- Feed top cell types into Q9 (chromatin accessibility at ADORA loci)
- Feed top cell types into Q13 (CellOracle perturbation targets)
