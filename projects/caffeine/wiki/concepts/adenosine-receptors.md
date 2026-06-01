# Adenosine Receptors

## Summary

Adenosine receptors are the core caffeine-response machinery for normal exposure. Mapping where ADORA1, ADORA2A, ADORA2B, and ADORA3 are expressed is the first filter for predicting which cell types can directly respond to caffeine.

**ADORA** is the human gene-symbol prefix for **ADenosine Receptor** genes. The genes are named ADORA1, ADORA2A, ADORA2B, and ADORA3; the receptor proteins are often called A1, A2A, A2B, and A3.

For the receptor-by-receptor signal flow from adenosine binding to cAMP changes, see [cAMP signaling and ADORA cascades](camp-signaling.md). For the `Gi/o`, `Gs`, `Golf`, and `Gq` decoder, see [G-protein coupling](g-protein-coupling.md).

## Receptor Cheat Sheet

| Gene symbol | Protein shorthand | Coupling | Caffeine effect | Expected enriched contexts |
|---|---|---|---|---|
| ADORA1 | A1 | [Gi/o](g-protein-coupling.md) | Blocks adenosine-mediated inhibition of adenylyl cyclase; can raise cAMP relative to adenosine-bound state | brain, adipose, heart atria, kidney, testis |
| ADORA2A | A2A | [Gs/Golf](g-protein-coupling.md) | Blocks adenosine-mediated cAMP increase through A2A | striatal medium spiny neurons, T cells, macrophages, NK cells, endothelial cells, platelets |
| ADORA2B | A2B | [Gs](g-protein-coupling.md) | Blocks low-affinity stress/pathology-associated adenosine signaling | intestinal epithelium, endothelial cells, cardiac fibroblasts, astrocytes, dendritic cells |
| ADORA3 | A3 | [Gi/o](g-protein-coupling.md) | Blocks inhibitory Gi/o signaling | mast cells, neutrophils, macrophages, eosinophils, synovial tissue, lung |

## Modeling Implication

The same caffeine dose can push different cell types in different directions because receptor composition differs. A cell dominated by ADORA1 does not have the same cAMP response as a cell dominated by ADORA2A.

## First Project Use

[001 adora expression](../labs/001-adora-expression.md) maps ADORA receptor expression across human cell types using CellxGene Census. That output should seed later chromatin and perturbation analyses.

New to the terms here (Gi/o, Gs, cAMP, coupling)? Start with [pharmacology vocabulary](pharmacology-vocabulary.md), then read [G-protein coupling](g-protein-coupling.md) and [cAMP signaling and ADORA cascades](camp-signaling.md).

Related pages: [G-protein coupling](g-protein-coupling.md), [cAMP signaling](camp-signaling.md), [pharmacology vocabulary](pharmacology-vocabulary.md), [caffeine molecular targets](caffeine-molecular-targets.md), [cell type response model](cell-type-response-model.md), [public data landscape](public-data-landscape.md)
