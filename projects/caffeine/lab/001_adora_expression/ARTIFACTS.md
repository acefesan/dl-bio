# Lab 001 — Artifact Walkthrough

A guided tour of every figure and cache file in this lab. Read top-to-bottom and you have the full Q1 story; read by section and you can jump to whatever you want to understand.

For the lab question, hypothesis, and per-receptor interpretation, see [entry.md](entry.md). For the broader project goals, see [PROPOSAL.md](../../PROPOSAL.md).

## How To Read This Guide

Each entry is shaped like:

> **`<filename>`** — one-line gist
>
> **What it shows.** What the axes/columns/colors mean.
>
> **How to read it.** What pattern to look for.
>
> **What it's good for.** Why it exists.
>
> **Caveats.** What to *not* read into it.

Color legend used throughout:
- **Blue** in our combined figures = Tabula Sapiens (peripheral atlas, 1.14M cells, 75 tissues, no brain)
- **Red** = HBCA non-neuronal (Human Brain Cell Atlas v1.0, brain glia + endothelial/etc., 888k cells, no neurons)

---

## 1. Source Data On Disk

These are the inputs to everything below.

### `cache/tabula_sapiens_all_cells.h5ad`

The 45 GB Tabula Sapiens "All Cells" H5AD from Quake lab via CELLxGENE Census 2025-11-08. 1,136,218 human cells × 60,606 genes across 24 donors and 75 tissues. Healthy tissue only. Includes precomputed PCA, scVI latent, and UMAP embeddings, plus harmonized cell-type and tissue ontology labels.

This file is too big to load with `anndata.read_h5ad` without OOM. The scripts in this folder use `h5py` to read targeted columns instead.

### `cache/tabula_sapiens_adora_expression.npz`

The four ADORA columns extracted from the H5AD above. Shape `(1136218, 4)`, float32, CSR-densified. About 500 KB. Has keys `expression`, `gene_indices`, `genes`, `matrix_group`. This is the compact reusable form — load it instead of re-reading the 45 GB H5AD.

### `cache/human_brain_cell_atlas/hbca_all_non_neuronal_b165f033.h5ad`

The HBCA non-neuronal release, ~4.4 GB. 888,263 brain cells × 58,232 genes. Contains glia (oligodendrocytes, astrocytes, OPCs, microglia, Bergmann glia), endothelial, pericyte, ependymal, choroid plexus, fibroblast, vascular smooth muscle. **No neurons.**

### `cache/human_brain_cell_atlas/hbca_adora_expression.npz`

Same shape as the Tabula one but for HBCA, 888,263 × 4. About 300 KB.

### `cache/human_brain_cell_atlas/README.md`

Provenance for the HBCA download — dataset_id, source URI, file size, plus a note about the matching neurons release (`8e10f1c4-…`, 30 GB) that is the obvious next download.

### `cache/tabula_sapiens_donor_summary.csv`

24-row table of TSP donor metadata (age, sex, ethnicity, cell count, tissue coverage). Useful for sanity checks like "is the top hit driven by one donor?"

---

## 2. The Story-Telling Figures (Q1 Pseudobulk Outputs)

These are the figures that answer "which cell types express ADORA receptors?" Built by `make_q1_summary_artifacts.py`.

### `figures/ranked_top20_per_receptor.png`

> Top 20 cell types per ADORA receptor by mean expression among expressing cells, combined across Tabula Sapiens (blue) and HBCA non-neuronal (red).

**What it shows.** 4-panel grid, one per receptor. Each bar = one cell type. Bar length = mean ADORA value in cells that actually express the gene (`expression > 0`). Cell type label includes the source (T or H) and cell count.

**How to read it.** Top of the list is "strongest expression among expressing cells" — not the same as "highest fraction expressing." A small-n cell type with one strong outlier can climb the list; those are noted in the entry.md caveats. For real signals, look for entries with large `n=` and big bars (HBCA oligodendrocytes / astrocytes are the cleanest examples for ADORA1).

**What it's good for.** Spotting candidates for Q9 (chromatin accessibility) follow-up. A cell type that ranks high here is worth chasing in matched scATAC data.

**Caveats.** Mean-among-expressing is sensitive to outliers when prevalence is tiny. The dotplot (Section 3) is the better all-in-one read. Brain neuron signals are absent because HBCA neurons aren't downloaded.

### `figures/donor_stratified_dotplot.png`

> Tabula Sapiens ADORA expression by top 20 cell types × 24 donors.

**What it shows.** Dotplot. Columns are (gene × donor); rows are cell types. Dot size = percent expressing in that cell-type-and-donor cohort; dot color = mean expression. Vertical grey lines separate the four ADORA gene blocks.

**How to read it.** For each cell type, scan horizontally across one gene's block. If the dot pattern is consistent across donors, the signal is donor-independent (real). If one donor has a fat dot and the rest are empty, that donor is driving the cell-type-level signal alone (sample-effect risk).

**What it's good for.** "Is one donor responsible for the macrophage ADORA3 signal?" answerable in seconds.

**Caveats.** Cell types with very few cells per donor will show tiny or absent dots — small-cell-type bias is real. Not all 24 donors contribute to all cell types because each donor contributed different tissues.

### `figures/assay_stratified_dotplot.png`

> Same as above, but stratified by assay (10x 3' v3, 10x 5' v2, Smart-seq2, Smart-seq3) instead of donor.

**What it shows / how to read it.** Same channel encoding as donor dotplot. The interesting check: do plate-based assays (Smart-seq2/3) and droplet assays (10x 3' v3, 5' v2) agree on which cell types express ADORA? Disagreement = the receptor signal is assay-dependent.

**What it's good for.** Sanity check against assay-specific dropout. Adenosine receptors are GPCRs and GPCRs often have low transcript counts — Smart-seq2's deeper per-cell coverage may detect them where 10x droplets miss them.

**Caveats.** Most cells (1.03M of 1.14M) are 10x 3' v3, so other assay columns are visually faint by sample size alone.

---

## 3. The Original Tabula Sapiens Figures

These predate the Q1 summary pass. They cover only Tabula Sapiens (peripheral) and use Tabula's precomputed embeddings.

### `figures/tabula_sapiens_adora_dotplot_cell_type.png` + `.csv`

> ADORA expression by cell type, top 35 ADORA-enriched cell types from Tabula Sapiens.

**What it shows.** Dotplot. Rows = cell types (ranked by total ADORA detection); columns = four ADORA genes. Dot size = percent expressing; dot color = mean expression.

**How to read it.** This is the canonical Q1 plot. Big bright dot = common AND strong. Small dark dot = rare AND weak. Big pale dot = many cells weakly express it. Small bright dot = few cells strongly express it. The CSV is the same data as a table.

**What it's good for.** First go-to for "what cell type expresses receptor X?"

**Caveats.** Tabula has no brain, so canonical ADORA brain story (neurons) is absent. Cell-type labels can hide tissue/donor mixtures — see donor-stratified plot.

### `figures/tabula_sapiens_adora_dotplot_tissue.png` + `.csv`

> Same as above but by `tissue_in_publication` (28 broad tissues).

**How to read it.** Quick anatomical orientation. Tongue light up for ADORA2B → look at tongue cell-type breakdown (next entry).

**Caveats.** "Tissue" is a blunt instrument. The tongue ADORA2B signal turns out to be in basal cells, not taste cells. Always drill down.

### `figures/tabula_sapiens_tongue_adora_cell_type_breakdown.png` + `.csv`

> Drilldown: of Tabula's tongue cells, which cell types express ADORA2B?

**What it shows.** Dotplot restricted to tongue cells, rows = tongue cell types, columns = four receptors. Reveals that the ADORA2B "tongue" signal is mainly in **basal cells** (35% expressing) and stratified squamous epithelium (14%) — not taste receptor cells (6%, tiny n=32).

**What it's good for.** Textbook demonstration of why a tissue dotplot needs a cell-type drilldown. Used as the worked example in the [interpretation page](../../wiki/labs/001-adora-interpretation.md#tongue-adora2b-example).

### `figures/tabula_sapiens_adora_high_umap.png`

> UMAP of all 1.14M Tabula cells, colored by ADORA expression with a high-expression threshold.

**What it shows.** UMAP coordinates from `obsm/X_umap`. Cells with low/zero ADORA values are rendered very faintly; cells above the per-receptor threshold (top quartile among expressing cells) are rendered brightly.

**How to read it.** Look for spatial neighborhoods where the colored cells cluster. A bright island = a transcriptional state where ADORA is reliably high. Diffuse coloring = no strong pattern; the receptor is sprinkled across many cell states.

**What it's good for.** Visual sanity check. The dotplot is more interpretable for sparse genes, but the UMAP reveals when ADORA-positive cells form coherent populations vs scatter.

**Caveats.** UMAP geometry is not biology. Distance between far-apart islands is not meaningful. Cluster shape is not biological.

### `figures/tabula_sapiens_adora_high_umap_four_color.png`

> Same UMAP, but with one color per ADORA receptor (green/orange/purple/pink).

**What it shows.** ADORA1 high = teal, ADORA2A high = orange, ADORA2B high = purple, ADORA3 high = pink. Spatial niches for each receptor are visible.

**How to read it.** The four receptors occupy different parts of the UMAP — pink (ADORA3) clusters in the immune blob, purple (ADORA2B) sits in epithelial regions, teal/orange more diffuse. This is the "do the receptors hit different cell types?" question answered visually.

**Companion JSON:** `tabula_sapiens_adora_high_umap_four_color_summary.json` records the per-receptor thresholds (ADORA1 ≥ 0.684, ADORA2A ≥ 1.080, ADORA2B ≥ 0.779, ADORA3 ≥ 0.943) and cell counts above threshold.

### `figures/tabula_sapiens_embeddings_by_tissue.png` (and `_fine.png`)

> Six-panel comparison of the embeddings Tabula Sapiens ships: PCA, scVI latent, plain UMAP, scVI-corrected UMAP, two uncorrected UMAPs. All 1.14M cells colored by tissue.

**What it shows.** Same cells, six different 2D projections. The "uncorrected" UMAPs are what the data looks like without batch correction — tissues smear into each other because donor/assay effects masquerade as cell-state. The scVI-corrected UMAP cleanly separates tissues.

**How to read it.** A demonstration of why batch correction matters. Use the scVI-corrected `X_umap_scvi_full_donorassay` (third in the grid) as the "right" embedding for downstream coloring.

**`_fine.png` variant** uses 75 fine-grained tissue labels (`tissue`) instead of 28 broad ones (`tissue_in_publication`).

### `figures/tabula_sapiens_X_umap_interactive_cell_type.html`

53 MB interactive Plotly UMAP. Open in a browser; hover for per-cell metadata. Heavy file — not committed to git. Regenerate with `make_tabula_interactive_embedding.py` if needed.

---

## 4. The HBCA Brain Figures

Same shape as the Tabula figures but on HBCA non-neuronal. Useful for brain comparisons.

### `figures/hbca_adora_dotplot_cell_type.png` + `.csv`

> Dotplot of ADORA expression across the 11 HBCA non-neuronal cell types.

**How to read it.** Rows are HBCA's 11 broad labels (astrocyte, oligodendrocyte, OPC, microglial cell, endothelial cell, pericyte, ependymal cell, choroid plexus epithelial cell, fibroblast, vascular smooth muscle cell, Bergmann glial cell). Columns are the four ADORA receptors. Same dot encoding as Tabula plots.

**What it's good for.** Reading the brain glial ADORA story. ADORA1 in oligodendrocytes and OPCs jumps out. ADORA2B in Bergmann glia (cerebellum) and astrocytes. ADORA3 in microglia.

### `figures/hbca_adora_dotplot_supercluster_term.png` + `.csv`

> Same as above but using HBCA's `supercluster_term` (10 categories, slightly different grouping than `cell_type`).

**When to use which:** `cell_type` is more conventional; `supercluster_term` is HBCA's own coarser grouping. Compare both for robustness.

### `figures/hbca_adora_expression_summary.json` + `.csv`

Per-receptor cell counts and percent-expressing across all HBCA, plus the total cell count and gene indices used. Read this if you want the headline HBCA numbers in machine-readable form.

---

## 5. Tabula × HBCA Comparison Figures

Where the two atlases meet. Built by `compare_*.py` and `make_shared_umap_*.py`.

### `figures/compare_pseudobulk_adora_cell_type_overlap.png` + `.csv`

> Cell types present in both Tabula and HBCA, with their ADORA pseudobulk values side by side.

**What it shows.** Matched cell-type labels between the two atlases (e.g., "astrocyte" appears in both). Compares mean expression and percent expressing for each ADORA receptor across the two sources.

**How to read it.** If the two atlases agree on ADORA expression in a cell type, the signal is replicable. If they disagree wildly, one atlas has a confound (typically: HBCA's astrocytes are brain-specific; Tabula's astrocyte cells come from non-brain tissues and may not exist as a coherent group).

**Companion CSVs:**
- `compare_pseudobulk_adora_ranked.csv` — full ranked cross-atlas pseudobulk
- `compare_pseudobulk_adora_all_summaries.csv` — every cell type from both atlases
- `compare_pseudobulk_adora_hbca_roi_top.csv` — HBCA's `ROI` (anatomical region) ranking
- `compare_pseudobulk_adora_hbca_supercluster_top.csv` — HBCA supercluster ranking

### `figures/compare_umap_tabula_hbca_cell_type.png`

> Side-by-side UMAPs of Tabula (left) and HBCA (right), each in its own embedding space, colored by cell_type.

**How to read it.** Each atlas's UMAP is independent. The two panels share nothing structurally — same shape there isn't meaningful. Use this to compare the *kinds* of cell types each atlas covers.

### `figures/compare_umap_tabula_hbca_adora_high.png`

> Same two panels but colored by ADORA-high cells using the same per-receptor threshold scheme as the four-color UMAP.

**How to read it.** Look for whether brain-resident cell types (HBCA right panel) have ADORA-positive cells in the same spatial niches as their peripheral counterparts (Tabula left panel).

### `figures/compare_umap_tabula_hbca_summary.json`

The summary metadata for the comparison: per-atlas cell counts, label vocabularies, gene indices.

### `figures/compare_labels_summary.md` + `.json` + supporting CSVs

> Reconciling Tabula's and HBCA's label vocabularies — which labels exist in both, which only in one.

**How to read.** Open `compare_labels_summary.md` first. It tells you that of HBCA's 11 cell types, 4 have exact matches in Tabula and 0 have similar-but-not-exact matches. Useful before you compare any per-cell-type signal across the two atlases.

**Supporting CSVs:**
- `compare_labels_category_counts.csv` — cell counts per (atlas, column)
- `compare_labels_cell_type_overlap.csv` — which cell_type strings match across atlases
- `compare_labels_obs_columns.csv` — every obs column present in each atlas

### `figures/shared_umap_tabula_hbca/` (and `_harmony`, `_smoketest`, `_harmony_smoketest`)

> Joint UMAPs built on cells from both atlases at once, using only the genes shared between them.

**Four variants:**
- `shared_umap_tabula_hbca/` — no batch correction
- `shared_umap_tabula_hbca_harmony/` — Harmony batch-corrected
- `*_smoketest/` — small sample versions for quick iteration

**Files inside each directory:**
- `shared_umap_tabula_hbca_adora_high.png` — joint UMAP, ADORA-high coloring
- `shared_umap_tabula_hbca_dataset.png` — joint UMAP, colored by source (Tabula vs HBCA)
- `shared_umap_tabula_hbca_cells.csv` — per-cell metadata for the joint embedding
- `shared_umap_tabula_hbca_arrays.npz` — the raw 2D coordinates
- `shared_umap_tabula_hbca_genes.txt` — list of genes used (only genes present in both atlases)
- `shared_umap_tabula_hbca_summary.json` — run config and per-atlas counts

**How to read.** The Harmony version is what you want for "do brain glia from HBCA and the closest matches in Tabula land in the same neighborhood?" The uncorrected version shows what happens when you skip batch correction (atlases separate into two blobs).

---

## 6. Q1 Summary Cache Tables (Built By `make_q1_summary_artifacts.py`)

These feather tables are the structured machine-readable form of the Q1 result.

### `cache/pseudobulk_by_cell_type.feather`

> Mean expression, mean expression among expressing cells, and percent expressing per gene per cell type per source atlas.

**Columns:**
- `source` — "Tabula Sapiens" or "HBCA non-neuronal"
- `cell_type`
- `n_cells`
- `ADORA1_mean`, `ADORA1_mean_nonzero`, `ADORA1_pct_expressing` — and the same for ADORA2A, ADORA2B, ADORA3

**Load it:**
```python
import pandas as pd
pb = pd.read_feather('cache/pseudobulk_by_cell_type.feather')
pb.sort_values('ADORA3_mean_nonzero', ascending=False).head(20)
```

191 rows total (Tabula's 180 + HBCA's 11). The canonical Q1 result in table form.

### `cache/cross_receptor_overlap.feather`

> Per cell type, count of cells positive for each subset of the four ADORA receptors.

**Columns:**
- `source`, `cell_type`, `n_cells`
- `ADORA1_pos`, `ADORA2A_pos`, `ADORA2B_pos`, `ADORA3_pos` — singles
- `ADORA1+ADORA2A`, `ADORA1+ADORA2B`, … — pairs (6 columns)
- `ADORA1+ADORA2A+ADORA2B`, … — triples (4 columns)
- `ADORA1+ADORA2A+ADORA2B+ADORA3` — the quadruple

**How to use it:**
```python
import pandas as pd
ovl = pd.read_feather('cache/cross_receptor_overlap.feather')
# Which cell types co-express ADORA1 and ADORA3?
ovl.sort_values('ADORA1+ADORA3', ascending=False).head(10)
```

Microglia is the textbook multi-receptor population — see entry.md for the headline read.

### `cache/q1_summary.json`

Run metadata for the Q1 artifact pass: cell counts per atlas, gene list, row counts, the "what's missing" note about HBCA neurons.

---

## 7. Other Cache Files (Less Central)

### `cache/stratified_100_full_*`

A 100-cell coord-based proof-of-concept extraction done early in the project to verify that coord pushdown queries work fast (40 seconds for 100 cells × all 60k genes). Not part of Q1; kept as a reference of the technique. See [001 v3 stratified fetch](../../wiki/labs/001-v3-stratified-fetch.md) for context.

### `cache/_obs_human_primary_normal.parquet` + `_sample_metadata.parquet`

Outputs of the v3 SOMA fetch pipeline (Phase 1: global obs scan; Phase 2: stratified sample). Phase 3 (X read) was abandoned. The parquet caches are still valid and can feed a future SOMA-based fetch if we ever want to revisit. For now, the Tabula + HBCA H5AD path made these unnecessary.

### `cache/fetch.log` + `stats.jsonl`

Cumulative human and structured logs from every fetch script invocation. Append-only; oldest entries from the v1 SOMA experiments. Useful for forensics, not for analysis.

### `cache/download_tabula.log`

Tiny log of the Tabula Sapiens download timing.

---

## 8. What's Missing (And Where Each Belongs)

| Open item | Where it would land |
|---|---|
| HBCA neurons (30 GB H5AD) — closes brain ADORA1/2A | `cache/human_brain_cell_atlas/hbca_neurons_8e10f1c4.h5ad` |
| Heart Cell Atlas — closes cardiac atrial ADORA1 | `cache/heart_cell_atlas/` |
| Liver Cell Atlas — needed for Q5 CYP1A2 | `cache/liver_cell_atlas/` |
| Findley 2019 HUVEC + caffeine — Q10 | `cache/findley_2019_huvec/` |
| GTEx bulk pseudobulk cross-check | `cache/gtex_v8_bulk_adora.csv` |
| Joint Tabula + HBCA-neurons + Tabula brain regions UMAP | `figures/joint_full_brain_umap/` |

For wiki context on what each of those questions is and why we'd want them, see [research questions](../../wiki/concepts/research-questions.md).

---

Related: [entry.md](entry.md), [DATASET_PROCESSING.md](DATASET_PROCESSING.md), [wiki Lab 001 ADORA expression](../../wiki/labs/001-adora-expression.md), [wiki Lab 001 ADORA interpretation](../../wiki/labs/001-adora-interpretation.md), [wiki scRNA visualization and analysis](../../wiki/concepts/scrna-visualization-and-analysis.md).
