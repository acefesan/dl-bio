# SOMA Axes And X

## Summary

In the Census, `X` is the expression matrix. It is not a dataframe. It is a sparse array whose coordinates point into two axis dataframes:

```text
obs.soma_joinid  -> X dim 0  -> cells
var.soma_joinid  -> X dim 1  -> genes/features
```

The axis meaning is established by the SOMA `Experiment` and `Measurement` structure:

```text
Experiment: homo_sapiens
├── obs                         cell axis dataframe
└── ms["RNA"]                   RNA measurement
    ├── var                     gene/feature axis dataframe
    └── X["raw"]                sparse expression array
```

So the alignment is not happening in pandas. It is part of the SOMA data model:

- `obs` describes the observation axis of the experiment,
- `var` describes the variable axis of a measurement,
- `X` stores values at coordinates from those two axes.

## The Three Pieces

For human RNA:

```python
human = census["census_data"]["homo_sapiens"]
rna = human.ms["RNA"]

obs = human.obs
var = rna.var
x_raw = rna.X["raw"]
```

These are different object types:

| Object | Type | Grain | Role |
|---|---|---|---|
| `human.obs` | SOMA `DataFrame` | one row per cell | cell metadata and cell-axis IDs |
| `rna.var` | SOMA `DataFrame` | one row per gene/feature | gene metadata and gene-axis IDs |
| `rna.X["raw"]` | SOMA `SparseNDArray` | one stored value per nonzero expression entry | expression values |

Observed schema for the 2025-11-08 human RNA objects:

```text
obs:
  soma_joinid: int64
  dataset_id, cell_type, tissue, disease, donor_id, ...

var:
  soma_joinid: int64
  feature_id, feature_name, feature_type, feature_length, ...

X["raw"]:
  soma_dim_0: int64
  soma_dim_1: int64
  soma_data: float
```

Read `X["raw"]` rows as sparse triples:

```text
(soma_dim_0, soma_dim_1, soma_data)
```

For this measurement:

```text
soma_dim_0 = a cell's obs.soma_joinid
soma_dim_1 = a gene's var.soma_joinid
soma_data  = expression value for that cell/gene pair
```

## A Tiny Toy Example

Imagine `obs` has three cells:

| obs.soma_joinid | cell_type | tissue |
|---:|---|---|
| 10 | astrocyte | brain |
| 11 | neuron | brain |
| 12 | T cell | blood |

And `var` has two genes:

| var.soma_joinid | feature_name |
|---:|---|
| 100 | ADORA1 |
| 101 | ADORA2A |

Then `X["raw"]` might store only nonzero values:

| soma_dim_0 | soma_dim_1 | soma_data |
|---:|---:|---:|
| 10 | 100 | 2.3 |
| 11 | 101 | 0.7 |
| 12 | 101 | 1.1 |

Interpreted through the axes, this means:

| cell | gene | expression |
|---|---|---:|
| astrocyte / brain | ADORA1 | 2.3 |
| neuron / brain | ADORA2A | 0.7 |
| T cell / blood | ADORA2A | 1.1 |

The sparse array does not itself store `cell_type` or `feature_name`. It stores coordinates. The axis dataframes give those coordinates meaning.

## What Layer Does The Alignment?

The alignment is defined at the SOMA `Experiment` / `Measurement` layer.

```text
Experiment
  owns obs, the observation axis

Measurement
  owns var, the variable axis for that modality
  owns X layers whose dimensions use obs and var coordinates
```

The dataframe layer provides axis metadata. The sparse array layer provides values. The experiment/measurement layer says which axis dataframe belongs to which array dimension.

This is why a standalone `SparseNDArray` schema only says:

```text
soma_dim_0
soma_dim_1
soma_data
```

By itself, that schema does not say "cell" or "gene". Inside `human.ms["RNA"].X["raw"]`, SOMA convention and object placement say:

```text
dim 0 -> human.obs
dim 1 -> human.ms["RNA"].var
```

## How Filters Become Coordinates

When we write:

```python
obs_value_filter = "tissue_general == 'brain' and cell_type == 'astrocyte'"
var_value_filter = "feature_name in ['ADORA1', 'ADORA2A']"
```

the conceptual flow is:

1. Filter `obs` to find matching cell `soma_joinid`s.
2. Filter `var` to find matching gene `soma_joinid`s.
3. Read `X` entries whose `soma_dim_0` is in the selected cell IDs and whose `soma_dim_1` is in the selected gene IDs.
4. Return the result as sparse triples, an Arrow table, or an `AnnData`, depending on the API path.

In other words, filters on human-readable metadata become coordinate selections on `X`.

## Why This Matters For Lab 001

Lab 001 asks for a tiny gene set across many cells:

```text
many obs.soma_joinid values
few var.soma_joinid values
```

Logically, that is a small result: many cells x four genes. Physically, the Census `X` array is cell-major on S3. That means selecting a few genes does not necessarily avoid reading cell chunks that contain all genes for those cells.

This is the bridge between two wiki pages:

- this page explains the logical axis alignment,
- [TileDB-SOMA storage](tiledb-soma-storage.md) explains the physical fragment layout and why some logical slices are slow.

## Relationship To AnnData

When `get_anndata()` returns an `AnnData`, it materializes the same logical relationship in a more familiar object:

```text
SOMA human.obs              -> adata.obs
SOMA rna.var                -> adata.var
SOMA rna.X["raw"] slice     -> adata.X
```

In AnnData, row order and column order are already chosen for the returned slice. You usually inspect `adata.obs` and `adata.var` directly rather than manually joining by `soma_joinid`.

But the source structure still matters: the AnnData was assembled from axis-coordinate reads against SOMA.

## Practical Introspection Cells

These cells are safe structure checks; they do not read the expression matrix values:

```python
human = census["census_data"]["homo_sapiens"]
rna = human.ms["RNA"]

print(human.obs.schema)
print(rna.var.schema)
print(rna.X["raw"].schema)
```

To inspect a few ADORA gene axis rows:

```python
adora_var = (
    rna.var.read(
        value_filter="feature_name in ['ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3']",
        column_names=["soma_joinid", "feature_id", "feature_name"],
    )
    .concat()
    .to_pandas()
)
```

To inspect a small cell metadata sample:

```python
brain_obs = (
    human.obs.read(
        value_filter="tissue_general == 'brain' and disease == 'normal'",
        column_names=["soma_joinid", "dataset_id", "cell_type", "tissue_general"],
    )
    .concat()
    .to_pandas()
)
```

Be careful with broad `obs` reads: metadata is cheaper than expression, but full-organism `obs` is still large.

## Common Confusions

`X` is not a dataframe.
: `X` is a sparse array. Its rows are stored coordinate/value entries, not cell metadata rows.

`obs.soma_joinid` is not the same thing as pandas index.
: It is the cell-axis coordinate used by SOMA. A pandas dataframe materialized from `obs` may also have a pandas index, but that is separate.

`var.soma_joinid` is not the gene symbol.
: It is the gene-axis coordinate. Use `feature_name` for symbols such as `ADORA1`; use `feature_id` for Ensembl-style IDs.

`dataset_id` does not index `X`.
: `dataset_id` is provenance metadata in `obs`. It can filter cells, but `X` coordinates still use `obs.soma_joinid`.

Related pages: [Census core objects](census-core-objects.md), [TileDB-SOMA storage](tiledb-soma-storage.md), [001 CellxGene Census API](../labs/001-cellxgene-census-api.md), [001 H5AD and AnnData cache](../labs/001-h5ad-anndata-cache.md)
