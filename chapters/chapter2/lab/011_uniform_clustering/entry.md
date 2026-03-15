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
