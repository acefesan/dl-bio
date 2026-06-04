# Census Experiment Tree

## Summary

The human Census object is a SOMA `Experiment`: a container for one organism's cell metadata, measurement metadata, and expression arrays. In Lab 001, the relevant experiment is:

```text
census["census_data"]["homo_sapiens"]
```

Observed for Census release `2025-11-08`:

```text
homo_sapiens                           Experiment
├── obs                                DataFrame
└── ms                                 Collection
    └── RNA                            Measurement
        ├── var                        DataFrame
        ├── feature_dataset_presence_matrix
        └── X                          Collection
            ├── raw                    SparseNDArray
            └── normalized             SparseNDArray
```

The short version:

| Object | Question it answers |
|---|---|
| `obs` | What cells exist, and what metadata describe them? |
| `ms["RNA"].var` | What genes/features exist on the RNA measurement axis? |
| `ms["RNA"].X["raw"]` | What raw expression values were measured for cell-gene pairs? |
| `ms["RNA"].X["normalized"]` | What normalized expression values are available for cell-gene pairs? |
| `feature_dataset_presence_matrix` | Which genes/features are present in which source datasets? |

For the object vocabulary behind `Experiment`, `Measurement`, `DataFrame`, and `SparseNDArray`, see [Census core objects](census-core-objects.md). For how the axes line up with `X`, see [SOMA axes and X](soma-axes-and-x.md).

## `obs`

`obs` is the observation dataframe. In single-cell data, an observation is usually one cell or nucleus.

In the human Census:

```text
human.obs
```

has one row per human cell. Its columns describe where the cell came from, what biological label it carries, how it was assayed, and basic quality-summary fields.

Examples:

| Column | Meaning |
|---|---|
| `soma_joinid` | cell-axis coordinate used by `X` |
| `cell_type` | harmonized cell type label |
| `tissue_general` | broad tissue label |
| `dataset_id` | source dataset provenance |
| `is_primary_data` | whether this is the primary Census copy of the cell |

The full column glossary is in [Census obs columns](census-obs-columns.md).

## `ms`

`ms` means measurements. It is a collection of modalities measured on the experiment's observation axis.

For this project, the active measurement is:

```text
human.ms["RNA"]
```

That means the measurement describes RNA expression for the cells in `human.obs`.

The Census can contain other measurement types in other contexts, but Lab 001 uses RNA expression only.

## `var`

`var` is the variable dataframe for a measurement. In RNA, the variables are genes or gene-like features.

```text
human.ms["RNA"].var
```

has one row per RNA feature. Its columns describe gene identifiers, gene symbols, feature type, and feature-level summary counts.

Examples:

| Column | Meaning |
|---|---|
| `soma_joinid` | gene-axis coordinate used by `X` |
| `feature_id` | stable feature identifier, usually Ensembl-style |
| `feature_name` | human-readable gene symbol such as `ADORA1` |
| `nnz` | number of nonzero expression entries for that feature |

The full column glossary is in [Census var columns](census-var-columns.md).

## `X`

`X` is the expression matrix collection for a measurement.

```text
human.ms["RNA"].X["raw"]
human.ms["RNA"].X["normalized"]
```

Each layer is a sparse cell x gene matrix. Stored entries are triples:

```text
(cell_soma_joinid, gene_soma_joinid, expression_value)
```

The sparse array does not store `cell_type` or `feature_name` directly. It stores coordinates. The coordinates point back to `obs.soma_joinid` and `var.soma_joinid`.

For the layer glossary, see [Census X layers and feature presence](census-x-layers-and-feature-presence.md).

## `feature_dataset_presence_matrix`

The presence matrix records which RNA features are present in which source datasets. This matters because not every source dataset measured or retained exactly the same gene set.

Use it when asking:

- Is a missing expression value biologically zero, or was the gene not measured in that dataset?
- Which datasets can support a given gene panel?
- Are all ADORA receptor genes available in the cells I am about to compare?

For Lab 001, this is a useful sanity-check object before interpreting absence of receptor expression as biology.

Related pages: [Census core objects](census-core-objects.md), [Census obs columns](census-obs-columns.md), [Census var columns](census-var-columns.md), [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [SOMA axes and X](soma-axes-and-x.md)
