# Cell-Type Response Model

## Summary

The project predicts caffeine responsiveness with four layers. Each layer filters or explains which cell types should respond even without direct caffeine perturbation data.

## Layer 1: Receptor Expression

Cells expressing [ADORA](adenosine-receptors.md) genes have direct molecular machinery for caffeine response. ADORA means the adenosine receptor gene family: ADORA1, ADORA2A, ADORA2B, and ADORA3. This is the purpose of [001 adora expression](../labs/001-adora-expression.md).

Output: cell type x receptor expression matrix.

## Layer 2: Receptor-Locus Accessibility

Cells with accessible [ADORA](adenosine-receptors.md) promoters or enhancers may be primed for receptor expression, even if mRNA is low in a static atlas.

Output: cell type x [ADORA](adenosine-receptors.md) regulatory accessibility matrix.

## Layer 3: TF Activity

Cells with accessible motifs for CREB, AP-1, NF-kB, Nrf2, AHR, HNF4A, MEF2, or NFAT may be positioned to execute caffeine-linked transcriptional programs.

Output: cell type x TF motif activity matrix.

## Layer 4: GRN Perturbation

CellOracle or related GRN tools can simulate receptor blockade or downstream TF perturbation to predict gene-expression response.

Output: predicted cell-type-specific caffeine response signatures.

## Scoring Sketch

A first caffeine-responsiveness score could combine:

- [ADORA](adenosine-receptors.md) expression percentile,
- [ADORA](adenosine-receptors.md) promoter/enhancer accessibility,
- caffeine-relevant TF motif activity,
- overlap with direct HUVEC caffeine-response genes,
- overlap with caffeine GWAS/EWAS regulatory annotations.

Related pages: [adenosine receptors](adenosine-receptors.md), [signaling to transcription](signaling-to-transcription.md), [computational pipeline](computational-pipeline.md), [research questions](research-questions.md)
