# Census X Layers And Feature Presence

## Summary

Inside the human RNA measurement, `X` contains expression matrices:

```text
human.ms["RNA"].X["raw"]
human.ms["RNA"].X["normalized"]
```

Each layer is a sparse cell x gene array. Rows are cell coordinates from `obs.soma_joinid`; columns are feature coordinates from `var.soma_joinid`.

The RNA measurement also contains:

```text
human.ms["RNA"].feature_dataset_presence_matrix
```

That object records which features are present in which source datasets.

For the full object tree, see [Census experiment tree](census-experiment-tree.md). For physical storage behavior, see [TileDB-SOMA storage](tiledb-soma-storage.md).

## Sparse Expression Triples

The `X` layers are not dataframes. They are sparse arrays.

At the table level, a stored nonzero looks like:

```text
soma_dim_0  soma_dim_1  soma_data
cell_id     gene_id     expression_value
```

Where:

| Sparse column | Meaning |
|---|---|
| `soma_dim_0` | cell coordinate, matching `obs.soma_joinid` |
| `soma_dim_1` | gene/feature coordinate, matching `var.soma_joinid` |
| `soma_data` | expression value in that X layer |

The readable metadata live outside `X`:

```text
obs.soma_joinid -> cell_type, tissue, disease, dataset_id, ...
var.soma_joinid -> feature_name, feature_id, feature_type, ...
```

## `raw`

`X["raw"]` is the raw expression layer exposed by the Census.

Use it when the analysis wants the least-transformed expression values available through the Census object. For Lab 001, raw expression is useful for sparse receptor detection summaries:

| Question | Raw layer use |
|---|---|
| Does this cell have any observed ADORA count? | Check nonzero raw expression. |
| Which cell types have detectable receptor expression? | Summarize nonzero fraction and raw signal by cell type. |
| Is expression extremely sparse? | Count nonzero entries per receptor. |

Raw counts are still shaped by source dataset processing, assay chemistry, depth, and Census harmonization. They are not magic ground truth.

## `normalized`

`X["normalized"]` is a normalized expression layer exposed by the Census.

Use it when comparing expression magnitude across cells or groups after a standard Census normalization. This is generally more convenient for plots and broad contrasts than raw counts, but it still inherits single-cell dropout and dataset composition effects.

For Lab 001, a sensible pattern is:

| Output | Layer |
|---|---|
| Detection fraction: percent of cells with expression above zero | `raw` |
| Mean expression among cells or groups | `normalized`, with raw as a sanity check |
| Robust dotplot | detection fraction from `raw`, color/intensity from `normalized` |

## `feature_dataset_presence_matrix`

The presence matrix answers a different question from expression:

```text
Was this feature present in this source dataset's measured feature set?
```

This matters because a zero-like absence in `X` can mean either:

| Situation | Interpretation |
|---|---|
| Feature present in dataset, no expression stored for a cell. | The cell likely had zero observed expression for that feature. |
| Feature absent from dataset. | The value is not biologically interpretable as zero. |

For ADORA work, the presence matrix can prevent a quiet but serious mistake: calling a cell type "ADORA-negative" when many of its source datasets did not include the receptor gene in the measured feature set.

## Coordinate Filters

The efficient expression read uses coordinates on both axes:

```python
cell_coords = sampled_obs["soma_joinid"].astype(int).tolist()
gene_coords = adora_var["soma_joinid"].astype(int).tolist()

cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    measurement_name="RNA",
    obs_coords=cell_coords,
    var_coords=gene_coords,
    obs_column_names=[...],
)
```

This differs from:

```python
var_value_filter="feature_name in ['ADORA1', ...]"
```

The value-filter form is easy to write, but the coordinate form matches TileDB's dimension indexing. See the pushdown section in [TileDB-SOMA storage](tiledb-soma-storage.md).

## Practical Rules

1. Treat `X` as sparse coordinate data, not as a dataframe with metadata columns.
2. Use `obs` and `var` to interpret `X` coordinates.
3. Prefer `var_coords` after resolving gene symbols to `var.soma_joinid`.
4. Use the presence matrix before interpreting missing receptor expression.
5. Keep both raw-detection and normalized-magnitude summaries when making cell-type plots.

Related pages: [Census experiment tree](census-experiment-tree.md), [Census obs columns](census-obs-columns.md), [Census var columns](census-var-columns.md), [SOMA axes and X](soma-axes-and-x.md), [TileDB-SOMA storage](tiledb-soma-storage.md), [001 v3 stratified fetch](../labs/001-v3-stratified-fetch.md)
