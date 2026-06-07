# Learning Path

This is the shortest ramp from "I know basic biology" to "I can contribute to this project."

## 1. Widen the Frame First

Read [hank caffeine video](raw/hank-caffeine-video.md), [purine alkaloids](concepts/purine-alkaloids.md), [xanthines](concepts/xanthines.md), [caffeine in plants](concepts/caffeine-in-plants.md), and [caffeine cultural history](concepts/caffeine-cultural-history.md).

Goal: understand why caffeine is not just "the coffee molecule." It is a convergently evolved plant metabolite, an ecological defense chemical, a pollinator-behavior molecule, and a culturally diverse exposure.

## 2. Human Mechanism

Read [pharmacology vocabulary](concepts/pharmacology-vocabulary.md), [caffeine molecular targets](concepts/caffeine-molecular-targets.md), [adenosine receptors](concepts/adenosine-receptors.md), [G-protein coupling](concepts/g-protein-coupling.md), [cAMP signaling and ADORA cascades](concepts/camp-signaling.md), and [secondary caffeine targets](concepts/secondary-caffeine-targets.md).

Goal: understand why normal caffeine exposure is mostly receptor antagonism, while many in vitro effects involve higher concentrations and need careful interpretation.

## 3. Drink Chemistry Is Not One Thing

Read [methylxanthines and tea effects](concepts/methylxanthines-and-tea-effects.md).

Goal: keep pure caffeine, coffee, tea, cacao, mate, guarana, yaupon, and synthetic caffeine analytically separate.

## 4. Connect Signaling to Gene Regulation

Read [signaling to transcription](concepts/signaling-to-transcription.md) and [epigenomics vocabulary](concepts/epigenomics-vocabulary.md).

Goal: be able to explain how cAMP, calcium, NF-kB, Nrf2, mTOR, CREB, NFAT, and HDAC export can show up as accessibility, histone, methylation, or expression signals.

## 5. Learn the Data Landscape

Read [public data landscape](concepts/public-data-landscape.md) and [direct caffeine epigenomics](concepts/direct-caffeine-epigenomics.md).

Goal: know which datasets are reference atlases versus true caffeine perturbations, and why the project is mostly an integration problem.

## 6. Understand the Predictive Model

Read [cell type response model](concepts/cell-type-response-model.md).

Goal: understand the four-layer prediction stack: receptor expression, receptor-locus accessibility, TF motif activity, and GRN perturbation.

## 7. Start With the Active Lab

Read [single-cell RNA-seq measurement](concepts/single-cell-rna-seq-measurement.md), [epithelial cell types](concepts/epithelial-cell-types.md), [immune cell types](concepts/immune-cell-types.md), [001 CellxGene Census API](labs/001-cellxgene-census-api.md), [001 data flow](labs/001-data-flow.md), [001 H5AD and AnnData cache](labs/001-h5ad-anndata-cache.md), [001 notebook guide](labs/001-notebook-guide.md), [001 adora expression](labs/001-adora-expression.md), and [001 ADORA interpretation](labs/001-adora-interpretation.md), then open `../lab/001_adora_expression/explore_adora_expression.ipynb`.

Goal: produce the first cell-type-by-receptor map and use it to prioritize later chromatin and GRN analyses.

## 8. Pick the Next Research Question

Read [caffeine sensitivity genetics](concepts/caffeine-sensitivity-genetics.md) and [research questions](concepts/research-questions.md).

Recommended next steps after Q1:

- Q9: chromatin accessibility at [ADORA](concepts/adenosine-receptors.md) loci.
- Q10: NFAT response across vascular cell types.
- Q13: CellOracle perturbation simulation.
