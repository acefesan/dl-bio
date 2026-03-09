# 006 — Sub-HOG Depth Sweep (Levels 1–3, Multi-Seed)

**Date:** 2026-03-07
**Models:** ESM2 150M, 650M, 3B
**Status:** complete

## Hypothesis
Entry 003 ran level-1 sub-HOG coloring with a single seed and found the top root HOGs
only cover ~2% of the sample. Going deeper (levels 2, 3) could reveal finer evolutionary
structure within those families. Running multiple seeds tests whether any apparent
clustering is stable or a UMAP initialization artifact.

This addresses open question #3: "does level-2/3 sub-HOG coloring reveal finer structure?"

## Setup
- **Script:** `chapters/chapter2/04_subtree_hog_analysis.py`
- **HOG levels tested:** 1, 2, 3
- **Top HOGs:** 4 (HOGs 801468, 136254, 792940, 801484)
- **min_subhog_size:** 3 (smaller groups than entry 003's threshold of 5)
- **Seeds:** 42, 123, 456, 789 (for 150M and 3B); 42 only for 650M
- **UMAP sources:**
  - 650M: `runs/default_20260216_213052/clustering/umap_coordinates.csv`
  - 150M: `runs/esm2_150m_20260306_225615/clustering/seed_N/umap_coordinates.csv`
  - 3B: `runs/esm2_3b_20260304_015127/clustering/seed_N/umap_coordinates.csv`
- **Annotations:** `runs/esm2_3b_20260304_015127/dataset/cafa3_annotations.feather`

Total: 27 runs (3 models × variable seeds × 3 levels), producing 135 figures + 27 metadata files.

## Results

### Figure organization

Every figure directory contains its own `metadata.json` with the exact command used,
the git commit, seed, and parameters. Structure:

```
figures/
├── 150m/seed_{42,123,456,789}/level_{1,2,3}/
├── 650m/seed_42/level_{1,2,3}/
└── 3b/seed_{42,123,456,789}/level_{1,2,3}/
```

Each directory contains:
- `subtree_hog_NNNNNN_levelN.png` — one per top root HOG
- `subtree_overview.png` — all 4 root HOGs highlighted
- `subtree_summary.json` — protein counts per HOG
- `metadata.json` — exact reproduction command

### Observations

Examine figures to assess:
1. Do sub-HOG members cluster together at level 2/3 (finer groups)?
2. Is the pattern consistent across seeds for the same model?
3. Does any model show sub-HOG clustering that the root-HOG analysis missed?

## Figures
- `figures/{model}/seed_{N}/level_{L}/subtree_hog_801468_level{L}.png` — largest HOG
- `figures/{model}/seed_{N}/level_{L}/subtree_overview.png` — all 4 HOGs highlighted
- 135 total PNGs across all combinations

## Interpretation
TODO: user should examine the figures across seeds and levels to determine:
- Whether deeper HOG levels reveal clustering invisible at root level
- Whether patterns are seed-stable (reliable) or artifacts

## Next steps
- If sub-HOG clustering appears: compute silhouette scores at sub-HOG level
- If no clustering: the embedding space may genuinely not organize by evolutionary
  lineage — try GO-term functional labels instead (open question #4)
- Consider UMAP of only HOG-annotated proteins (drop the gray majority) for clearer view
