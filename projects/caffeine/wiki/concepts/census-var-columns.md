# Census Var Columns

## Summary

`human.ms["RNA"].var` is the feature metadata table for the human RNA measurement. In RNA expression, a variable is usually a gene.

Observed path:

```text
census["census_data"]["homo_sapiens"].ms["RNA"].var
```

Each row describes one RNA feature. The corresponding expression values live in `X`, where the feature's `soma_joinid` is used as the gene-axis coordinate.

For the object tree, see [Census experiment tree](census-experiment-tree.md). For axis alignment, see [SOMA axes and X](soma-axes-and-x.md).

## Column Glossary

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `soma_joinid` | Internal feature-axis coordinate. This is the coordinate used by `X.soma_dim_1`. | Use for `var_coords` when doing dimension-pushed gene reads. |
| `feature_id` | Stable feature identifier, commonly Ensembl-style for genes. | Good machine key for genes across tools. |
| `feature_name` | Human-readable feature name, commonly gene symbol. | Use for readable gene panels such as `ADORA1`, `ADORA2A`, `ADORA2B`, `ADORA3`. |
| `feature_type` | Feature category, such as gene-like RNA feature. | Sanity-check that selected features are the expected kind. |
| `feature_length` | Feature length as recorded in the Census feature metadata. | Usually not needed for simple receptor expression summaries. |
| `nnz` | Number of nonzero expression entries for this feature across the Census matrix. | Rough global prevalence of detectable expression. |
| `n_measured_obs` | Number of observations/cells for which this feature was in the measured feature set. | Helps separate "not expressed" from "not measured". |

## `feature_name` Versus `feature_id`

For human-facing work, `feature_name` is easier:

```text
ADORA1
ADORA2A
ADORA2B
ADORA3
```

For exact joins and durable references, `feature_id` is safer because gene symbols can change, collide, or carry aliases across resources.

For storage-layer coordinate reads, neither `feature_name` nor `feature_id` is the coordinate. The coordinate is `var.soma_joinid`.

The efficient pattern is:

```python
adora_var = (
    rna.var.read(
        value_filter="feature_name in ['ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3']",
        column_names=["soma_joinid", "feature_id", "feature_name"],
    )
    .concat()
    .to_pandas()
)

gene_coords = adora_var["soma_joinid"].astype(int).tolist()
```

Then pass `gene_coords` as `var_coords` for the expression read when possible.

## Why `n_measured_obs` Matters

Single-cell source datasets do not always have identical feature sets. If a gene has no stored expression for a cell, there are two different possibilities:

| Case | Meaning |
|---|---|
| The gene was measured and has zero expression. | A biologically interpretable zero, within assay limits. |
| The gene was not in that source dataset's measured feature set. | Missing-by-design, not evidence of no expression. |

`n_measured_obs` gives a feature-level summary of how broadly that feature was measured. For dataset-level checks, use `feature_dataset_presence_matrix`; see [Census X layers and feature presence](census-x-layers-and-feature-presence.md).

## Practical Rules

1. Use `feature_name` to find familiar genes.
2. Keep `feature_id` in outputs for provenance.
3. Convert selected genes to `var.soma_joinid` when doing expression reads.
4. Check feature presence before interpreting missing expression too strongly.
5. Remember that `var.soma_joinid` only indexes the RNA feature axis; it is not comparable to `obs.soma_joinid`.

Related pages: [Census experiment tree](census-experiment-tree.md), [Census obs columns](census-obs-columns.md), [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [SOMA axes and X](soma-axes-and-x.md), [TileDB-SOMA storage](tiledb-soma-storage.md)
