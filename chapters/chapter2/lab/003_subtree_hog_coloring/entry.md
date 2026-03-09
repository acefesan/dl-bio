# 003 — Subtree HOG Coloring (Level-1 Sub-HOGs)

**Date:** 2026-03-06
**Model:** ESM2 650M (reusing UMAP from run 001)
**Status:** complete

## Hypothesis
The root HOG coloring in entry 001 was too coarse and showed most proteins as gray
(no root HOG). Going one level deeper in the HOG tree — coloring by level-1 sub-HOG
within the largest root HOGs — should reveal whether proteins in the same evolutionary
subfamily cluster together in the embedding space.

## Setup
- UMAP coordinates from: `runs/default_20260216_213052/clustering/umap_coordinates.csv`
  (11k sampled proteins, ESM2 650M embeddings, seed=42)
- HOG IDs from: `runs/esm2_3b_20260304_015127/dataset/cafa3_annotations.feather`
- Top 4 root HOGs by sample count: 801468 (203 proteins), 136254 (52), 792940 (41), 801484 (26)
- Level-1 sub-HOG coloring; sub-HOGs with <5 proteins grouped as "Other"
- Script: `chapters/chapter2/04_subtree_hog_analysis.py`

## Results

| Root HOG | Proteins in sample | % of sample | # Level-1 sub-HOGs |
|----------|-------------------|-------------|---------------------|
| 801468   | 203               | 1.8%        | many (diverse)      |
| 136254   | 52                | 0.5%        | fewer               |
| 792940   | 41                | 0.4%        | fewer               |
| 801484   | 26                | 0.2%        | few                 |

Key observation: even the largest root HOG (801468) contains only ~2% of the 11k sample,
meaning the HOG-annotated proteins are sparse in the UMAP.

## Figures
- `figures/subtree_overview.png` — all 4 root HOGs highlighted on single UMAP
- `figures/subtree_hog_801468_level1.png` — largest HOG, colored by L1 sub-HOG
- `figures/subtree_hog_136254_level1.png`
- `figures/subtree_hog_792940_level1.png`
- `figures/subtree_hog_801484_level1.png`

## Interpretation
TODO: examine figures to draw conclusions about whether sub-HOG members cluster.

## Next steps
- Try level-2 or level-3 sub-HOGs for finer resolution
- Compute quantitative silhouette scores at sub-HOG level (not just root HOG)
- Consider plotting only proteins WITH a HOG assignment (drop the gray background majority)
  to see the HOG structure more clearly on its own UMAP
- Run 150M model and compare sub-HOG clustering
