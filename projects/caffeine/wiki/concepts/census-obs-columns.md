# Census Obs Columns

## Summary

`human.obs` is the cell metadata table for the human Census experiment. Each row is one cell or nucleus. The columns answer three kinds of questions:

| Column family | Question |
|---|---|
| Identity | Which cell is this on the SOMA axis? |
| Provenance | Which dataset, donor, assay, and tissue did it come from? |
| Biology | What cell type, disease state, sex, ethnicity, and developmental stage are annotated? |
| QC summaries | How much RNA signal and how many genes were measured? |

Observed path:

```text
census["census_data"]["homo_sapiens"].obs
```

For where `obs` sits in the object tree, see [Census experiment tree](census-experiment-tree.md).

## Identity And Provenance

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `soma_joinid` | Internal cell-axis coordinate. This is the coordinate used by `X.soma_dim_0`. | Use for `obs_coords` and for joining returned expression rows back to `obs`. |
| `dataset_id` | UUID for the source CELLxGENE dataset/H5AD that contributed the cell. | Use to trace cells back to source studies and donor/dataset composition. |
| `observation_joinid` | Source-side observation identifier carried through Census ingest. | Useful if tracing a cell back to its source object; not usually needed for grouping. |
| `donor_id` | Donor identifier from the source dataset after Census harmonization. | Needed for donor-aware summaries and avoiding pseudoreplication. |
| `is_primary_data` | Boolean indicating the primary Census copy of a cell. | Filter to `True` for most analyses to avoid duplicate representation across datasets. |

`dataset_id` is the bridge to `census["census_info"]["datasets"]`. Do not join datasets to cells by `soma_joinid`; those coordinates are local to each SOMA dataframe or axis.

## Assay Columns

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `assay` | Human-readable assay label, such as a single-cell or single-nucleus RNA-seq protocol. | Stratify or sanity-check whether receptor detection differs by technology. |
| `assay_ontology_term_id` | Ontology identifier for the assay label. Usually an EFO-style term. | Stable machine-readable assay key. Better than string labels for joins. |
| `suspension_type` | Whether the source material is cells, nuclei, or another suspension category. | Important because single-nucleus RNA can under-detect some transcripts relative to whole-cell RNA. |

An ontology term ID is the controlled-vocabulary version of a label. The readable label is nice for plots; the ontology ID is safer for exact matching.

## Cell Type Columns

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `cell_type` | Harmonized readable cell type label. | Primary grouping variable for ADORA expression. |
| `cell_type_ontology_term_id` | Ontology ID for the cell type, usually from the Cell Ontology. | Stable key for grouping labels across releases or datasets. |

Example:

```text
cell_type = "mast cell"
cell_type_ontology_term_id = controlled vocabulary ID for mast cell
```

Use `cell_type` for human-facing tables and plots. Use the ontology ID when exact identity matters across datasets, releases, or synonym changes.

## Tissue Columns

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `tissue` | Specific tissue or anatomical structure label. | Fine-grained tissue comparisons. |
| `tissue_ontology_term_id` | Ontology ID for `tissue`. | Stable exact tissue key. |
| `tissue_general` | Broader tissue grouping. | Coarse stratification across the body. |
| `tissue_general_ontology_term_id` | Ontology ID for the broad tissue grouping. | Stable exact broad tissue key. |
| `tissue_type` | Category describing tissue context, such as tissue versus organoid/cell culture-like material where present. | Sanity-check whether cells come from primary tissue-like contexts. |

In Lab 001, `tissue_general` is useful for broad body-wide maps, while `tissue` is useful once a broad category needs to be split into anatomical subregions.

## Disease And Development Columns

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `disease` | Disease or health-state label. | Filter to normal controls or compare disease contexts. |
| `disease_ontology_term_id` | Ontology ID for disease. | Stable exact disease key. |
| `development_stage` | Age or developmental stage label. | Separate fetal, pediatric, adult, or other stage-specific patterns. |
| `development_stage_ontology_term_id` | Ontology ID for developmental stage. | Stable exact stage key. |

The v3 fetch uses:

```python
value_filter="is_primary_data == True and disease == 'normal'"
```

That keeps primary, normal cells before drawing the stratified sample.

## Donor Demographic Columns

| Column | Plain meaning | Lab 001 use |
|---|---|---|
| `sex` | Donor sex label as harmonized by the Census. | Possible covariate or QC split. |
| `sex_ontology_term_id` | Ontology ID for sex. | Stable exact key. |
| `self_reported_ethnicity` | Donor self-reported ethnicity label, where available. | Usually a covariate/provenance field, not a primary grouping variable here. |
| `self_reported_ethnicity_ontology_term_id` | Ontology ID for ethnicity label. | Stable exact key. |

These fields are often missing or unevenly distributed across studies. Treat them as useful provenance, not as uniformly sampled population metadata.

## Expression Summary Columns

These columns summarize the expression profile for each cell over measured genes.

| Column | Plain meaning | Interpretation |
|---|---|---|
| `raw_sum` | Sum of raw counts for the cell over measured variables. | Library size / total captured RNA signal. |
| `nnz` | Number of genes/features with nonzero expression in the cell. | Detected feature count. |
| `raw_mean_nnz` | Mean raw expression among nonzero entries. | Average positive signal, excluding zeros. |
| `raw_variance_nnz` | Variance of raw expression among nonzero entries. | Spread of positive signal, excluding zeros. |
| `n_measured_vars` | Number of features measured for that cell's source dataset. | Helps distinguish true zeros from genes absent from the measured feature set. |

These are not ADORA-specific. They describe each cell's overall expression capture and are useful for QC, assay comparisons, and detecting odd source datasets.

## Practical Rules

1. Use `soma_joinid` for coordinate reads, not for biological interpretation.
2. Use `dataset_id` for provenance and joining to dataset metadata.
3. Use readable labels for plots, ontology IDs for stable exact grouping.
4. Filter `is_primary_data == True` unless there is a specific reason to keep duplicate Census representations.
5. Interpret absence of expression only after checking whether the feature was measured in the source dataset.

Related pages: [Census experiment tree](census-experiment-tree.md), [Census core objects](census-core-objects.md), [SOMA axes and X](soma-axes-and-x.md), [Census X layers and feature presence](census-x-layers-and-feature-presence.md), [001 v3 stratified fetch](../labs/001-v3-stratified-fetch.md)
