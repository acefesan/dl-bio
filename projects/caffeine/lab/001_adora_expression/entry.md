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

Initial atlas sanity check:
- `plot_tabula_embeddings_by_tissue.py` reads the Tabula Sapiens H5AD through `h5py` and plots every `obsm` embedding without materializing `obs` as a pandas frame.
- `figures/tabula_sapiens_embeddings_by_tissue.png` colors all 1,136,218 cells by `obs/tissue_in_publication` (28 broad tissue labels).
- `figures/tabula_sapiens_embeddings_by_tissue_fine.png` colors all cells by `obs/tissue` (75 finer anatomical labels).
- `figures/tabula_sapiens_adora_dotplot_tissue.png` summarizes ADORA expression by broad tissue (`obs/tissue_in_publication`).
- `figures/tabula_sapiens_adora_dotplot_cell_type.png` summarizes ADORA expression by cell type.
- Interactive Plotly/WebGL explorers can be regenerated locally with `make_tabula_interactive_embedding.py`; the HTML outputs are intentionally not committed because they are large.
- `cache/tabula_sapiens_donor_summary.csv` summarizes the 24 donor IDs represented in the H5AD.

Q1 summary artifacts (built by `make_q1_summary_artifacts.py` on 2026-06-08):
- `figures/ranked_top20_per_receptor.png` — top 20 cell types per receptor (mean expression among expressing cells), now combined across **three atlases**: Tabula Sapiens, HBCA non-neuronal, and HBCA neurons (2,480,956 cells, grouped by `supercluster_term`).
- `figures/donor_stratified_dotplot.png` — Tabula × donor sanity check.
- `figures/assay_stratified_dotplot.png` — Tabula × assay sanity check.
- `cache/pseudobulk_by_cell_type.feather` — (mean, mean_nonzero, pct_expressing) per gene per cell type per source, across all three atlases (HBCA neurons grouped by `supercluster_term`).
- `cache/cross_receptor_overlap.feather` — per cell type, count of cells positive for each subset of the four receptors (now also three-atlas).

For a guided walkthrough of every artifact, see [ARTIFACTS.md](ARTIFACTS.md).

## Interpretation

**Receptor-by-receptor headline reads (from the combined pseudobulk):**

- **ADORA1.** The neuronal expectation is now **confirmed**: with HBCA neurons folded in, the top ADORA1 neuron superclusters are hippocampal CA4 (mean 2.58 among expressing cells, 68.4% expressing, n=10,654), thalamic excitatory (2.41, 60.2%, n=85,546), deep-layer corticothalamic and 6b (2.03, 56.2%, n=78,396), deep-layer near-projecting (2.01, 62.3%), and deep-layer intratelencephalic (1.86, 53.2%, n=228,467). This matches the proposal's prior of strong ADORA1 in cortex, hippocampus, and thalamus, both in mean transcript level and in prevalence (over half of these neurons carry ADORA1 transcript). The glial half of the brain ADORA1 story still holds: HBCA oligodendrocytes (n=494,966, 20% expressing, mean 1.31 in expressing cells), oligodendrocyte precursor cells (n=105,734, 23% expressing), astrocytes (n=155,025, 14% expressing). Tabula non-brain hits (platelets, vestibular supporting cells) sit at very low percent-expressing (<1%) and are likely small-sample noise rather than real biology.
- **ADORA2A.** The canonical striatal site is now **confirmed and dominant**: medium spiny neurons are the #1 ADORA2A neuron supercluster (mean 2.62 among expressing cells, 36.5% expressing, n=152,189) — by far the strongest ADORA2A signal in any atlas. The next neuronal hits drop off sharply (hippocampal CA4 1.14/12.9%, thalamic excitatory 1.12/8.9%, hippocampal CA1-3 1.09/7.6%, amygdala excitatory 1.08/6.6%), so the striatal MSN signal is specific, not a generic neuronal background. On the non-neuronal side, HBCA pericytes (n=3,693, 25% expressing) remain the cleanest vascular signal, and Tabula immune signals (thymocytes, T cells) appear but at low prevalence.
- **ADORA2B.** Distributed across glia (HBCA Bergmann glial cells in cerebellum 18% expressing, astrocytes 11% expressing) and Tabula barrier epithelium (tongue basal cells, bronchial smooth muscle, bladder, trachea). Matches the proposal's expectation: ADORA2B in intestinal/airway epithelium and astrocytes.
- **ADORA3.** The most coherent peripheral signal. Tabula myeloid populations — macrophages, microglia, monocytes, neutrophils, mast cells — and HBCA microglia all rank high. Plus glial cells outside the CNS (Mueller cells in retina, Schwann cells in peripheral nerve). This matches the textbook ADORA3 = "myeloid + mast cell + glia" picture.

**Cross-receptor co-expression:**

The most striking pattern is that **myeloid cells, especially microglia and macrophages, are the multi-receptor population.** Of 410 microglia in Tabula, 2 cells express *all four* ADORA receptors and 3 cells express three receptors at once. Macrophages: 29 cells co-express ADORA1+2B+3 out of 69,072. No other broad cell class shows quadruple-positive cells. This is biologically interesting — myeloid cells appear primed for adenosine-system integration rather than relying on one receptor subtype.

**Donor- and assay-stratified dotplots:** the patterns are consistent across donors (TSP1–30) and across assays (10x 3' v3 dominates, but Smart-seq2 sees comparable expression where overlap exists). No single donor or assay is driving the top-cell-type signals.

**Atlas-coverage note:** broad tissue labels in the local Tabula H5AD do not include brain, but HBCA is now represented by **both halves** — the non-neuronal release (glia/vascular) and the neuronal release (2,480,956 cells, grouped by `supercluster_term`). The strongest ADORA1/2A biological prior — neuronal expression in striatum (A2A medium spiny neurons), cortex/hippocampus/thalamus (A1) — is therefore now **confirmed** rather than untestable. The main residual caveat is assay homogeneity: the HBCA neuron file is 100% 10x 3' v3 (single-assay), so the neuronal numbers carry no cross-assay corroboration and inherit the dropout characteristics of droplet 3'-tag sequencing for sparse genes like the ADORAs.

**GTEx bulk pseudobulk comparison:** still TODO. The pseudobulk_by_cell_type feather is the right input for that join.

**Candidates for Q9 (chromatin accessibility at ADORA loci):** the cell types with the largest *and* most prevalent ADORA signal — HBCA oligodendrocytes, OPCs, astrocytes, microglia; Tabula macrophages, mast cells, bronchial smooth muscle, Bergmann glial cells — are the cells where matched scATAC would actually have a chance of showing ADORA promoter/enhancer accessibility.

## Next steps

- Compare against GTEx bulk pseudobulk to confirm direction
- Feed top cell types into Q9 (chromatin accessibility at ADORA loci)
- Feed top cell types into Q13 (CellOracle perturbation targets)
