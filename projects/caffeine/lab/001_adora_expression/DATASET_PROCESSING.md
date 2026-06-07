# Processing Tabula Sapiens for Lab 001

This lab uses the cached Tabula Sapiens all-cells H5AD:

```text
projects/caffeine/lab/001_adora_expression/cache/tabula_sapiens_all_cells.h5ad
```

It is large, about 45 GB on disk, with 1,136,218 cells and 60,606 genes. Do not casually open it with `scanpy.read_h5ad()` or `anndata.read_h5ad(..., backed="r")`: even backed mode materializes the full `obs` dataframe, which can use several GB of RAM.

## Environment

Run from the repo root:

```bash
uv run python projects/caffeine/lab/001_adora_expression/<script>.py
```

Use `uv add <package>` if a dependency is missing. Do not use `pip install`.

## Useful Files

- `cache/tabula_sapiens_all_cells.h5ad` — full Tabula Sapiens H5AD.
- `cache/stratified_100_full_raw.h5ad` — small stratified test file for script development.
- `cache/stratified_100_full_X_raw.npz` — sparse raw matrix for the small test file.
- `cache/stratified_100_full_obs.parquet` and `cache/stratified_100_full_var.parquet` — small test metadata.
- `cache/tabula_sapiens_donor_summary.csv` — derived donor summary table.
- `plot_tabula_embeddings_by_tissue.py` — example HDF5-direct reader and plotting script.
- `make_tabula_interactive_embedding.py` — interactive Plotly/WebGL embedding browser.
- `summarize_tabula_donors.py` — donor-level summary generator.

## H5AD Structure

Important groups in the full file:

- `X` — main expression matrix, stored sparse.
- `raw/X` — raw counts.
- `layers/decontXcounts` — DecontX ambient-RNA-corrected counts.
- `layers/scale_data` — log-normalized and scaled expression.
- `obs` — cell metadata, mostly AnnData categorical columns.
- `var` — gene metadata.
- `obsm` — precomputed embeddings.
- `obsp/connectivities` and `obsp/distances` — neighbor graph.
- `uns` — scVI metadata, palettes, and plotting hints.

Precomputed embeddings under `obsm`:

- `X_pca` — 1,136,218 x 50
- `X_scvi` — 1,136,218 x 50
- `X_umap` — 1,136,218 x 2
- `X_umap_scvi_full_donorassay` — 1,136,218 x 2
- `X_uncorrected_alltissues_umap` — 1,136,218 x 2
- `X_uncorrected_umap` — 1,136,218 x 2

Useful categorical `obs` columns:

- `tissue_in_publication` — 28 broad tissue labels, best for readable overview plots.
- `tissue` — 75 finer anatomical labels.
- `cell_type` — 180 cell-type labels.
- `broad_cell_class` — 40 broader cell classes.
- `donor_id`, `donor_tissue`, `donor_tissue_assay`, `assay`, `compartment`.

## Safe Access Pattern

Prefer reading only the arrays/columns you need with `h5py`.

```python
from pathlib import Path

import h5py

h5ad = Path("projects/caffeine/lab/001_adora_expression/cache/tabula_sapiens_all_cells.h5ad")

with h5py.File(h5ad, "r") as f:
    umap = f["obsm/X_umap"][:]          # 1.1M x 2, okay in RAM
    codes = f["obs/tissue_in_publication/codes"][:]
    categories = [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in f["obs/tissue_in_publication/categories"][:]
    ]
```

AnnData categoricals are stored as:

```text
obs/<column>/codes
obs/<column>/categories
```

Codes are integer arrays. A code of `-1` means missing.

## Expression Access

For gene-expression questions, avoid loading all of `X` unless necessary. First find the gene indices in `var`, then read only those columns in chunks.

The ADORA genes used for Q1 are:

```text
ADORA1
ADORA2A
ADORA2B
ADORA3
```

If a quick whole-matrix prototype is needed, use the small stratified H5AD first. For full data, stream row chunks or selected columns from the sparse HDF5 representation.

Sparse CSR groups are usually shaped like:

```text
X/data
X/indices
X/indptr
```

Build only what you need:

```python
import h5py
import scipy.sparse as sp

shape = (1_136_218, 60_606)

with h5py.File(h5ad, "r") as f:
    X = sp.csr_matrix(
        (f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
        shape=shape,
    )
    adora = X[:, [11005, 5997, 12798, 55546]].toarray()
```

That pattern is acceptable for quick exploration on a roomy machine, but chunked aggregation is preferred for repeatable lab scripts.

## Current Visualization Outputs

Generate all embedding plots colored by broad tissue:

```bash
uv run python projects/caffeine/lab/001_adora_expression/plot_tabula_embeddings_by_tissue.py
```

Local output (ignored by git because it is large):

```text
projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_embeddings_by_tissue.png
```

Generate the finer 75-label tissue version:

```bash
uv run python projects/caffeine/lab/001_adora_expression/plot_tabula_embeddings_by_tissue.py \
  --color-by tissue \
  --out projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_embeddings_by_tissue_fine.png \
  --point-size 0.04 \
  --alpha 0.45
```

For fast drafts, add a subsample:

```bash
uv run python projects/caffeine/lab/001_adora_expression/plot_tabula_embeddings_by_tissue.py \
  --max-points 100000
```

Generate ADORA dotplots by cell type and by broad tissue:

```bash
uv run python projects/caffeine/lab/001_adora_expression/plot_adora_dotplot.py

uv run python projects/caffeine/lab/001_adora_expression/plot_adora_dotplot.py \
  --group-by tissue_in_publication \
  --min-cells 50 \
  --top-n 28 \
  --out projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_adora_dotplot_tissue.png \
  --table projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_adora_dotplot_tissue.csv
```

Generate a zoomable/clickable HTML explorer for `X_umap`, colored by cell type:

```bash
uv run python projects/caffeine/lab/001_adora_expression/make_tabula_interactive_embedding.py \
  --color-by cell_type \
  --out projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_X_umap_interactive_cell_type.html \
  --marker-size 2.0 \
  --max-marker-size 9.0
```

Local output (ignored by git because it is large):

```text
projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_X_umap_interactive_cell_type.html
```

The default interactive explorer samples 250,000 cells to keep the browser responsive. It includes cell index, donor, tissue, cell type, broad cell class, assay, sex, age/stage, ethnicity, disease, and compartment in the hover/click metadata. Marker size grows as you zoom in so local color structure stays visible; tune this with `--marker-size` and `--max-marker-size`. Use `--max-points 0` for all cells, but expect a much larger and slower HTML file.

Generate the same browser colored by broad tissue:

```bash
uv run python projects/caffeine/lab/001_adora_expression/make_tabula_interactive_embedding.py
```

Output:

```text
projects/caffeine/lab/001_adora_expression/figures/tabula_sapiens_X_umap_interactive.html
```

Generate the donor summary table:

```bash
uv run python projects/caffeine/lab/001_adora_expression/summarize_tabula_donors.py
```

Output:

```text
projects/caffeine/lab/001_adora_expression/cache/tabula_sapiens_donor_summary.csv
```

This file has 24 `donor_id` categories. A donor is one human tissue donor; many cells and often many tissues can come from the same donor.

## Agent Rules For This Dataset

- Keep all new scripts, notes, and figures under `projects/caffeine/lab/001_adora_expression/`.
- Put generated figures in `projects/caffeine/lab/001_adora_expression/figures/`.
- Prefer `h5py` and direct HDF5 reads for the full H5AD.
- Use the small stratified cache for development when possible.
- Do not rerun dimensionality reduction unless the question explicitly requires it; use the existing embeddings.
- Do not mix this work into `chapters/chapter2/`, which is the ESM2 protein-embedding workstream.
