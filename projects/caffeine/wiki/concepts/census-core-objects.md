# Census Core Objects

## Summary

The CZ CELLxGENE Census is easiest to understand as a remote object tree. The top-level object is a TileDB-SOMA `Collection`; inside it are more containers, table-like metadata objects, and matrix-like expression objects.

The confusing part is that several things look "dataframe-ish" at different stages:

| Word | In this project | Where it lives |
|---|---|---|
| SOMA `DataFrame` | remote, table-like TileDB-SOMA object | S3 / TileDB |
| `TableReadIter` | streaming reader returned by `.read()` | Python handle over remote read |
| Arrow `Table` | in-memory columnar table returned by `.concat()` | local Python memory |
| pandas `DataFrame` | familiar local dataframe returned by `.to_pandas()` | local Python memory |

So this chain:

```python
datasets_df = (
    census["census_info"]["datasets"]
    .read()
    .concat()
    .to_pandas()
)
```

means:

```text
remote SOMA DataFrame -> read iterator -> Arrow Table -> pandas DataFrame
```

## Object Tree

For the 2025-11-08 Census release, the top-level shape observed in Lab 001 is:

```text
census                                      Collection
├── census_info                             Collection
│   ├── datasets                            DataFrame
│   ├── organisms                           DataFrame
│   ├── summary                             DataFrame
│   └── summary_cell_counts                 DataFrame
├── census_data                             Collection
│   ├── homo_sapiens                        Experiment
│   │   ├── obs                             DataFrame
│   │   └── ms["RNA"]                       Measurement
│   │       ├── var                         DataFrame
│   │       └── X["raw"]                    SparseNDArray
│   ├── mus_musculus                        Experiment
│   └── other organisms...
└── census_spatial_sequencing               Collection
```

## Collection

A `Collection` is a named container. It is close to a directory or group: it holds child SOMA objects under string keys.

Examples:

```python
census.keys()
census["census_info"].keys()
census["census_data"].keys()
```

Collections can contain dataframes, experiments, arrays, or other collections. They do not themselves represent rows of biological data; they organize the object tree.

## DataFrame

A SOMA `DataFrame` is a remote table-like object. It has rows, columns, a schema, and a URI, but it is not a pandas dataframe until materialized.

Examples:

```python
census["census_info"]["datasets"]
census["census_data"]["homo_sapiens"].obs
census["census_data"]["homo_sapiens"].ms["RNA"].var
```

The most important Census dataframes for Lab 001 are:

| SOMA DataFrame | Grain | Purpose |
|---|---|---|
| `census_info["datasets"]` | one row per source dataset | source provenance: title, collection, DOI, H5AD path, total cells |
| `homo_sapiens.obs` | one row per human cell | cell metadata: tissue, cell type, disease, donor, dataset ID |
| `homo_sapiens.ms["RNA"].var` | one row per gene | gene metadata: feature ID, feature name, gene length |

## TableReadIter

Calling `.read()` on a SOMA DataFrame does not immediately return all the rows. It returns a `TableReadIter`, a streamable iterator over Arrow tables.

```python
datasets_iter = census["census_info"]["datasets"].read()
datasets_iter
```

That object means "the read has been described and can be consumed." This design matters because `obs` can contain hundreds of millions of cells. Loading the whole thing just because `.read()` was called would be dangerous.

For small metadata tables, it is normal to immediately collect the iterator:

```python
datasets_arrow = datasets_iter.concat()
datasets_df = datasets_arrow.to_pandas()
```

For large reads, especially `obs` or `X`, keeping the iterator chunked is often the safer pattern.

## Arrow Table

An Arrow `Table` is a local, in-memory, columnar table. It is what `.concat()` returns after consuming a `TableReadIter`.

This is the "table" meant in many SOMA examples:

```python
table = census["census_info"]["datasets"].read().concat()
```

It is not the remote storage object anymore. It is a local Arrow representation that can be converted to pandas:

```python
df = table.to_pandas()
```

Arrow tables are useful because they preserve column types efficiently and can represent dictionary-encoded columns from the Census without immediately turning everything into Python objects.

## pandas DataFrame

A pandas `DataFrame` is the local analysis object most Python users expect.

```python
datasets = (
    census["census_info"]["datasets"]
    .read()
    .concat()
    .to_pandas()
)
```

Use pandas for small-to-medium metadata, summaries, joins, display, and plotting. Do not casually materialize full organism `obs` or expression `X` into pandas.

## Experiment

A SOMA `Experiment` is the core single-cell bundle for one organism. It holds:

- `obs`: rows are cells,
- `ms`: measurements, such as RNA,
- each measurement's `var`: rows are features/genes,
- each measurement's `X`: matrix layers.

For human RNA:

```python
human = census["census_data"]["homo_sapiens"]
human.obs
human.ms["RNA"].var
human.ms["RNA"].X["raw"]
```

This mirrors the AnnData mental model:

```text
AnnData.obs  ~ SOMA Experiment.obs
AnnData.var  ~ SOMA Measurement.var
AnnData.X    ~ SOMA Measurement.X["raw"]
```

## Measurement

A `Measurement` is one modality inside an experiment. For Lab 001, the relevant measurement is `RNA`.

```python
rna = human.ms["RNA"]
```

The measurement owns the gene axis (`var`) and the expression layers (`X`). A Census can also contain other modalities or spatial sequencing objects, but Lab 001 only needs RNA expression.

## SparseNDArray

`human.ms["RNA"].X["raw"]` is a SOMA sparse array. Logically, it is a cell x gene expression matrix. Physically, TileDB stores it as sparse coordinates and values rather than a dense rectangle.

At the logical level, each nonzero value is like:

```text
(cell_soma_joinid, gene_soma_joinid, expression_value)
```

This is why `obs.soma_joinid` and `var.soma_joinid` matter: they are the coordinates that align metadata rows to matrix axes.

For a fuller explanation of how `obs`, `var`, and `X` align, see [SOMA axes and X](soma-axes-and-x.md). For more on the physical S3 layout and fragments, see [TileDB-SOMA storage](tiledb-soma-storage.md).

## `soma_joinid`

`soma_joinid` is an internal coordinate used to join SOMA axes to arrays. It is not a universal biological ID.

The meaning depends on the object:

| Location | What `soma_joinid` indexes |
|---|---|
| `homo_sapiens.obs.soma_joinid` | the cell axis of human expression matrices |
| `homo_sapiens.ms["RNA"].var.soma_joinid` | the gene/feature axis of RNA expression matrices |
| `census_info["datasets"].soma_joinid` | the row coordinate inside the dataset catalog table |
| `census_info["summary"].soma_joinid` | the row coordinate inside the summary table |

Do not join unrelated tables by `soma_joinid`. Join datasets to cells with `dataset_id`, not `soma_joinid`.

Correct:

```python
obs_with_dataset_titles = obs.merge(
    datasets[["dataset_id", "collection_name", "dataset_title"]],
    on="dataset_id",
    how="left",
)
```

Wrong:

```python
obs.merge(datasets, on="soma_joinid")
```

The only time `soma_joinid` is the right join key is when aligning an axis dataframe with its own matrix axis. For example, `obs.soma_joinid` aligns cells to rows of `X`, and `var.soma_joinid` aligns genes to columns of `X`.

## Dataset IDs Versus Collections

`dataset_id` and `collection_id` are provenance identifiers from CELLxGENE Discover:

| ID | Meaning |
|---|---|
| `dataset_id` | one source dataset / H5AD-style unit |
| `collection_id` | a Discover collection containing one or more datasets |

One collection can contain multiple datasets. Many cells in `obs` can share the same `dataset_id`, because they came from the same source dataset.

This is different from a SOMA `Collection`. Unfortunately, the same English word appears in two places:

- CELLxGENE **collection**: a biological/provenance grouping on the Discover portal.
- SOMA `Collection`: a storage/API container object.

Context decides which one is meant.

## Practical Rules For Lab 001

1. Use `census["census_info"]["datasets"]` to understand source datasets and provenance.
2. Use `human.obs` to ask which cells exist and to filter by tissue, disease, cell type, donor, and `dataset_id`.
3. Use `human.ms["RNA"].var` to map gene names like `ADORA1` to Census feature rows.
4. Use `human.ms["RNA"].X["raw"]` or `get_anndata()` to retrieve expression values.
5. Treat `soma_joinid` as an axis coordinate, not as a universal row ID.
6. Treat `.read()` as the start of a streaming read; `.concat().to_pandas()` is the moment data become a local dataframe.

Related pages: [SOMA axes and X](soma-axes-and-x.md), [001 CellxGene Census API](../labs/001-cellxgene-census-api.md), [TileDB-SOMA storage](tiledb-soma-storage.md), [001 H5AD and AnnData cache](../labs/001-h5ad-anndata-cache.md), [001 data flow](../labs/001-data-flow.md)
