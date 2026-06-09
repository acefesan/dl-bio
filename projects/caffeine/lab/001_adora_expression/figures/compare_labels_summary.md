# Label comparison: Tabula Sapiens vs HBCA non-neuronal

- Tabula `cell_type` categories: 180
- HBCA non-neuronal `cell_type` categories: 11
- Exact normalized HBCA-to-Tabula cell-type matches: 4 / 11
- Similar non-exact HBCA-to-Tabula cell-type matches: 0 / 11

## Best labels for ADORA comparisons

### Cell Identity

- Tabula: `cell_type`
- HBCA: `cell_type or supercluster_term`
- Why: Tabula has fine whole-body cell types, while HBCA non-neuronal has 11 broad cell_type labels and 10 supercluster_term labels. For ADORA summaries, use HBCA cell_type/supercluster_term inside the brain dataset and compare cautiously to Tabula cell_type.

### Broad Class

- Tabula: `broad_cell_class`
- HBCA: `supercluster_term`
- Why: These are the closest broad biological groupings across datasets, useful for comparing immune, endothelial, epithelial-like, and glial/non-neuronal compartments.

### Anatomical Context

- Tabula: `tissue_in_publication or tissue`
- HBCA: `ROIGroupCoarse, ROIGroupFine, roi, tissue`
- Why: Tabula tissue labels describe organs/tissues; HBCA labels describe brain regions. They should be presented side by side rather than normalized into one shared tissue axis.

## Notes

- Existing UMAP coordinates are dataset-specific and should not be merged as if they share axes.
- For ADORA pattern comparison, dotplots/pseudobulk tables by label are more defensible than comparing raw UMAP geometry.
- HBCA non-neuronal is brain-specific; Tabula Sapiens is whole-body and lacks brain tissue in this cache.
