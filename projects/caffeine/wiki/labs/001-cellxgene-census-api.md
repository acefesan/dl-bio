# Lab 001: CellxGene Census API

## Summary

Lab 001 uses the **CZ CELLxGENE Discover Census** through the `cellxgene_census` Python package. The Census is a versioned data object and API for querying standardized single-cell RNA data from CZ CELLxGENE Discover without manually downloading and harmonizing hundreds of source datasets.

In this lab, we use it to answer one narrow question:

> Which human cell types express [ADORA](../concepts/adenosine-receptors.md) receptor genes?

## Who Built This

CZ CELLxGENE is a single-cell data platform from the **Chan Zuckerberg Initiative**. CZI describes CZ CELLxGENE as a platform for scientists to access, analyze, and annotate high-dimensional single-cell data at scale.

Useful official sources:

- [CZ CELLxGENE Census documentation](https://chanzuckerberg.github.io/cellxgene-census/)
- [Census data and schema](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_schema.html)
- [Learning about the Census notebook](https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_census_info.html)
- [CZI Cell Science tools and funding page](https://chanzuckerberg.com/science/programs-resources/single-cell-biology/)

## What the API Is

The `cellxgene_census` Python package is a convenience layer around **TileDB-SOMA**. SOMA stands for **stack of matrices, annotated**. That name is literal: the data are organized as matrices plus cell and gene annotations. For how those matrices are physically partitioned on S3 and which query shapes that layout favors, see [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md).

The API gives us:

- metadata queries: "which cells/datasets/tissues exist?"
- gene metadata queries: "which gene symbols and Ensembl IDs are in the Census?"
- expression matrix slicing: "give me these genes for these cells"
- export into familiar Python objects such as `AnnData`

## Mental Model

```mermaid
flowchart TD
    Portal[CZ CELLxGENE Discover portal] --> Census[CZ CELLxGENE Census]
    Census --> SOMA[TileDB-SOMA data object]
    SOMA --> Obs[obs: cell metadata]
    SOMA --> Var[var: gene metadata]
    SOMA --> X[X matrix: cell x gene expression]
    Obs --> Query[Filter cells]
    Var --> Query2[Filter genes]
    Query --> Slice[Expression slice]
    Query2 --> Slice
    Slice --> AnnData[AnnData object in Python]
```

## What `open_soma()` Does

In the notebook:

```python
census = cellxgene_census.open_soma()
```

This opens the Census object. It does not mean "download the whole Census." It gives the notebook a handle to query remote, versioned data.

The script version uses:

```python
cellxgene_census.open_soma(context=ctx)
```

with a TileDB-SOMA context that increases S3 timeouts because the live query streams sparse matrix data and can fail on slow networks.

## The Important Top-Level Paths

| Path | Meaning |
|---|---|
| `census["census_info"]["summary"]` | high-level Census summary metadata |
| `census["census_info"]["datasets"]` | dataset-level metadata |
| `census["census_data"]["homo_sapiens"]` | human data object |
| `census["census_data"]["homo_sapiens"].obs` | cell metadata |
| `census["census_data"]["homo_sapiens"].ms["RNA"].var` | gene metadata for the RNA measurement |

## What `get_anndata()` Does

In the notebook/script, the key call is:

```python
cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    measurement_name="RNA",
    var_value_filter="feature_name in ['ADORA1','ADORA2A','ADORA2B','ADORA3']",
    obs_value_filter="tissue_general == 'brain' and disease == 'normal'",
    obs_column_names=["cell_type", "tissue", "tissue_general", "assay", "dataset_id", "donor_id"],
)
```

Read this as:

- organism: human cells only,
- measurement: RNA expression,
- gene filter: only the four ADORA genes,
- cell filter: only selected cells, such as normal brain cells,
- obs columns: bring along the metadata needed for grouping and sanity checks.

The returned object is an `AnnData`.

## What AnnData Is

`AnnData` is the standard Python object used by Scanpy and many single-cell workflows.

For this lab:

| AnnData slot | Meaning here |
|---|---|
| `adata.X` | sparse cell x gene expression matrix for the selected ADORA genes |
| `adata.obs` | cell metadata, one row per cell |
| `adata.var` | gene metadata, one row per gene |

So if `adata.shape == (200000, 4)`, that means 200,000 cells and four genes.

When the fetch script calls `adata.write_h5ad(...)`, that in-memory AnnData object is serialized as a local `.h5ad` file. The notebook can later reopen it with `anndata.read_h5ad(...)` without calling the live Census API. For the file layout and cache behavior, see [H5AD and AnnData cache](001-h5ad-anndata-cache.md).

## Why This API Is Useful Here

Without the Census, we would need to:

- identify many separate human tissue datasets,
- download each dataset,
- normalize metadata names,
- harmonize cell type labels,
- map gene IDs,
- concatenate matrices,
- track provenance.

The Census does much of that standardization upfront. Lab 001 can therefore focus on the biological question: where are ADORA receptors expressed?

Related pages: [Lab 001 overview](001-adora-expression.md), [Lab 001 data flow](001-data-flow.md), [H5AD and AnnData cache](001-h5ad-anndata-cache.md), [ADORA receptors](../concepts/adenosine-receptors.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [001 fetch stall post-mortem](001-fetch-stall-postmortem.md)
