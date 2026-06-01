# Public Data Landscape

## Summary

Most datasets available for this project are baseline reference atlases, not direct caffeine perturbations. The value comes from joining them with the few caffeine-specific datasets.

## Key Resources

| Resource | What it contributes |
|---|---|
| ENCODE | ATAC-seq, ChIP-seq, TF binding, histone marks, Hi-C in tissues and cell lines |
| Roadmap Epigenomics | 127 reference epigenomes and ChromHMM states |
| GEO | Direct caffeine datasets, EWAS datasets, treatment expression studies |
| Human Cell Atlas / CellxGene Census | single-cell expression and cell-type labels across tissues |
| GTEx v8 | bulk tissue expression and eQTLs |
| CistromeDB / ChIP-Atlas | TF binding evidence near genes like ADORA2A or CYP1A2 |
| JASPAR / HOCOMOCO | transcription factor motifs |
| FANTOM5 | promoter and enhancer usage |
| 4D Nucleome | 3D chromatin context |
| LINCS L1000 / CMap | compound expression signatures |
| PharmGKB / DrugBank | pharmacogenomics and drug-target metadata |

## Practical Split

Use CellxGene/GTEx first for expression mapping. Use Roadmap/ENCODE/HCA scATAC second for chromatin context. Use GEO direct perturbation data as the anchor for what a real caffeine response looks like. Use GWAS/EWAS for population-level constraints.

Related pages: [direct caffeine epigenomics](direct-caffeine-epigenomics.md), [gwas ewas pharmacogenomics](gwas-ewas-pharmacogenomics.md), [cell type response model](cell-type-response-model.md), [001 adora expression](../labs/001-adora-expression.md)
