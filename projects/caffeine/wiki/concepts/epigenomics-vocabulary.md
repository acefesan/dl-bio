# Epigenomics Vocabulary

For the basic pharmacology terms (receptor, ligand, gene, transcription factor), see [pharmacology vocabulary](pharmacology-vocabulary.md). This page owns the genomics, regulatory, and assay terms.

## Genomics and Regulatory Basics

| Term | Plain meaning |
|---|---|
| Chromatin | DNA plus the proteins (histones) it wraps around. How tightly it is packed controls whether genes can be read |
| Chromatin accessibility | How "open" a stretch of chromatin is. Open regions can be bound by regulators and tend to be active; closed regions are silenced |
| Promoter | The stretch of DNA right at a gene's start where transcription is switched on — the gene's "on ramp" |
| Enhancer | A regulatory stretch of DNA, often far from the gene, that boosts how strongly the gene is transcribed |
| Motif | A short, recurring DNA sequence pattern that a transcription factor recognizes and binds. "CREB motif" = the pattern CREB looks for |
| DNA methylation | A small chemical tag (a methyl group) added to DNA, usually at CpG sites; often dials gene activity down. What WGBS and EWAS measure |
| Locus | A specific location in the genome (plural: loci). "ADORA loci" = where the ADORA genes sit |
| Pseudobulk | Adding up single-cell measurements within a cell type to imitate a bulk sample — steadier numbers per cell type |
| TPM | Transcripts Per Million: a normalized unit for expression so samples can be compared fairly |
| Dotplot | A grid plot where dot color = average expression and dot size = fraction of cells expressing the gene. Used to show ADORA across cell types |

## Core Assays

| Term | Meaning | Project use |
|---|---|---|
| ATAC-seq | Measures open chromatin using Tn5 transposase insertion | Identify accessible [ADORA](adenosine-receptors.md) regulatory elements and caffeine-responsive regions |
| DNase-seq | Measures open chromatin through DNase I cutting | Baseline accessibility in ENCODE/Roadmap |
| ChIP-seq | Maps protein or histone-mark binding to DNA | Histone marks and TF binding at caffeine-relevant loci |
| WGBS | Whole-genome bisulfite methylation profiling | Genome-wide methylation if available |
| EWAS | Association scan between methylation sites and traits | Coffee/caffeine methylation associations |
| eQTL | Variant associated with expression variation | Connect caffeine GWAS variants to gene regulation |
| ChromHMM | Segments genome into chromatin states from histone marks | Interpret enhancer/promoter states across cell types |
| GRN | Gene regulatory network | Predict downstream response to receptor perturbation |

## Histone Marks

- H3K4me3: active promoters.
- H3K27ac: active enhancers and promoters.
- H3K4me1: poised or active enhancers.
- H3K27me3: Polycomb-associated repression.
- H3K9me3: constitutive heterochromatin.

## Project Translation

The project is not only asking whether genes change expression. It asks whether caffeine-relevant receptors, enhancers, TF motifs, methylation sites, and chromatin states identify cell types that are primed to respond.

Related pages: [signaling to transcription](signaling-to-transcription.md), [public data landscape](public-data-landscape.md), [computational pipeline](computational-pipeline.md)
