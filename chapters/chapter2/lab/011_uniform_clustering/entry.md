# 011 — Uniform-Truncation Multi-Model Clustering

**Date:** 2026-03-15
**Models:** ESM2 150M, 650M, 3B
**Status:** complete

## Hypothesis
Entry 006 compared models but each used different UMAP coordinates (computed from
separately-embedded protein sets with different max-seq-length filters). The 650M
set had 902 extra long-sequence proteins. This entry re-does the comparison on
a **uniform protein set** (141,091 common proteins) and drops HOG-based metrics
(class sizes too small and unbalanced). Instead we measure KMeans silhouette in
both the original high-D embedding space and in UMAP 2D space.

## Setup
- **Script:** `chapters/chapter2/05_model_size_clustering.py`
- **Common protein set:** 141,091 proteins present in all 3 embedding files
  (650M filtered from 141,993 → 141,091 by intersection)
- **Sampling:** 500 proteins × 22 taxa = 11,000 per run
- **KMeans:** k=20, n_init=10, cosine silhouette in original space
- **UMAP:** n_neighbors=15, min_dist=0.1, cosine metric, euclidean silhouette in 2D
- **Seeds:** 42, 123, 456, 789 (4 runs per model = 12 total)

## Results

| Model | Dim  | Avg Sil (orig) | Std   | Avg Sil (UMAP) | Std   |
|-------|------|-----------------|-------|-----------------|-------|
| 150M  | 640  | **0.0518**      | 0.004 | **-0.067**      | 0.028 |
| 650M  | 1280 | 0.0370          | 0.003 | -0.091          | 0.007 |
| 3B    | 2560 | 0.0001          | 0.003 | -0.174          | 0.009 |

Per-seed breakdown:

| Model | Seed 42 | Seed 123 | Seed 456 | Seed 789 |
|-------|---------|----------|----------|----------|
| 150M orig | 0.0495 | 0.0471 | 0.0520 | 0.0584 |
| 650M orig | 0.0410 | 0.0337 | 0.0384 | 0.0350 |
| 3B orig   | 0.0045 | -0.0009 | -0.0010 | -0.0022 |
| 150M UMAP | -0.068 | -0.048 | -0.112 | -0.038 |
| 650M UMAP | -0.099 | -0.079 | -0.089 | -0.095 |
| 3B UMAP   | -0.161 | -0.172 | -0.178 | -0.184 |

### Per-taxon silhouette breakdown (seed=42, original space, cosine)

| Taxon              | 150M    | 650M    | 3B      |
|--------------------|---------|---------|---------|
| P. falciparum      | 0.0926  | **0.2130** | 0.0180  |
| T. brucei          | 0.0828  | 0.0944  | 0.0135  |
| X. laevis          | 0.0758  | 0.0366  | 0.0106  |
| G. gallus          | 0.0672  | 0.0289  | 0.0003  |
| S. pombe           | 0.0627  | 0.0449  | 0.0119  |
| O. sativa          | 0.0604  | 0.0414  | 0.0105  |
| M. musculus        | 0.0588  | 0.0197  | -0.0039 |
| H. sapiens         | 0.0568  | 0.0398  | 0.0124  |
| S. cerevisiae      | 0.0554  | 0.0527  | 0.0072  |
| D. melanogaster    | 0.0522  | 0.0309  | -0.0086 |
| D. rerio           | 0.0506  | 0.0184  | -0.0012 |
| C. elegans         | 0.0483  | 0.0354  | 0.0025  |
| D. discoideum      | 0.0289  | 0.0885  | 0.0100  |
| A. thaliana        | 0.0335  | -0.0007 | -0.0107 |
| E. coli            | 0.0180  | -0.0031 | -0.0006 |
| E. nidulans        | 0.0299  | -0.0051 | -0.0003 |
| M. tuberculosis    | 0.0293  | 0.0218  | -0.0077 |
| P. aeruginosa      | 0.0203  | 0.0189  | 0.0048  |
| **All 22 taxa positive?** | **yes** | no (3 neg) | no (5 neg) |

*P. falciparum* (malaria parasite) is the strongest signal in all models —
phylogenetically distant from everything else. The 650M score (0.213) is 4× its
overall mean, suggesting the model concentrates discriminative power on the most
divergent organisms. The 3B model barely distinguishes any taxon (max 0.018).

## Interpretation
**Inverse scaling confirmed on uniform protein set.** Filtering to common proteins
does not change the story — 150M > 650M >> 3B consistently across all 4 seeds.

New finding: **UMAP silhouette is negative for all models.** KMeans clusters found
in high-D don't correspond to compact regions in UMAP 2D. The 3B model is worst
(-0.174), suggesting its KMeans "clusters" are especially incoherent geometrically.

The 150M UMAP silhouette has the highest variance (σ=0.028 vs 0.007–0.009),
indicating its UMAP structure is more seed-dependent — the clusters that form are
real but fragile under different initializations.

## Next steps
- Visualize UMAP plots colored by KMeans cluster for each model×seed to see the geometry
- Try PCA whitening before clustering to test the anisotropy hypothesis (open question #1)
- Run KMeans directly in UMAP space (not just evaluate high-D labels there)
