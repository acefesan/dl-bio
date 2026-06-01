# Proposal Source

Source file: `../../PROPOSAL.md`

## Summary

The proposal argues that caffeine likely rewires gene regulation across many tissues, but direct caffeine-treatment epigenomic experiments are rare. The project should therefore predict caffeine's regulatory effects by integrating:

- baseline epigenomic atlases,
- caffeine-responsive gene lists,
- adenosine receptor expression,
- GWAS and EWAS signals,
- signaling pathway knowledge,
- single-cell RNA and ATAC atlases,
- gene regulatory network inference.

## Key Claims

- Physiological caffeine primarily acts as a non-selective antagonist of four adenosine receptors: ADORA1, ADORA2A, ADORA2B, and ADORA3.
- Higher caffeine concentrations also affect [PDEs, ryanodine receptors, ATM/ATR, and mTOR](../concepts/secondary-caffeine-targets.md), but many of those effects are supraphysiological relative to normal coffee intake.
- Caffeine perturbation reaches chromatin through transcription-factor pathways such as CREB, AP-1, NF-kB, Nrf2, AHR, HNF4A, MEF2, and NFAT.
- Findley et al. 2019 and Boye et al. 2024 are the central direct endothelial caffeine epigenomics datasets.
- No systematic caffeine ChIP-seq, Hi-C, single-cell multiome, or atlas-scale GRN perturbation project exists yet.

## Most Important Gaps

- No genome-wide histone ChIP-seq after caffeine treatment.
- No caffeine-treatment 3D genome data.
- No single-cell multi-omics under caffeine exposure.
- No systematic ATAC-seq time course.
- EWAS evidence is largely blood-only.
- NFAT/calcium signaling from endothelial work has not been mapped across other cell types.

## Links

Related pages: [caffeine molecular targets](../concepts/caffeine-molecular-targets.md), [direct caffeine epigenomics](../concepts/direct-caffeine-epigenomics.md), [computational pipeline](../concepts/computational-pipeline.md), [research questions](../concepts/research-questions.md)
