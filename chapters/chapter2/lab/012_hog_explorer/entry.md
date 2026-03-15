# 012 — HOG Hierarchy Explorer

**Date:** 2026-03-15
**Model:** N/A (exploration notebook)
**Status:** complete

## Hypothesis
Previous entries used root-HOG labels for clustering evaluation, but those classes
are extremely small and unbalanced (entry 011 dropped them for this reason). Can we
find a sub-HOG level within a single large root HOG that provides a well-balanced
set of classes suitable for clustering evaluation?

## Setup
- **Notebook:** `chapters/chapter2/lab/012_hog_explorer/hog_explorer.ipynb`
- **Annotations:** `cafa3_annotations.feather` (92,160 proteins with hog_id)
- **Focus HOG:** 801468 (2,236 proteins, max depth 14, by far the largest)

## Results

### HOG 801468 balance by level (all proteins, min_class=5)

| Level | Classes | Coverage    | Min | Max | Mean | CV   |
|-------|---------|-------------|-----|-----|------|------|
| 1     | 113     | 1070/2236   | 5   | 65  | 9.5  | 0.80 |
| 2     | 84      | 662/2236    | 5   | 25  | 7.9  | 0.47 |
| 3     | 53      | 426/2236    | 5   | 19  | 8.0  | 0.46 |
| 4     | 43      | 304/2236    | 5   | 18  | 7.1  | 0.39 |
| 5     | 33      | 217/2236    | 5   | 18  | 6.6  | 0.36 |
| 6     | 30      | 182/2236    | 5   | 9   | 6.1  | 0.22 |

**Level 6 has the best balance** (CV=0.22, sizes 5-9) but only 8% coverage.
**Level 3 is the best trade-off**: 53 classes, 19% coverage, CV=0.46.

### Key structural finding
Most proteins (541/2236 = 24%) terminate at depth 1 — they belong to a level-1
sub-HOG with no further subdivision. Only 553 (25%) reach depth >= 5. The HOG
hierarchy is very top-heavy: most evolutionary diversification captured by OMA
happens in the first 1-2 splits.

### No other HOG comes close
HOG 801468 has 553 depth>=5 proteins; the next largest (792940) has only 109.
For deep-level analysis, 801468 is the only viable candidate.

## Notebook contents
- Balance report function: query any HOG at any level
- Depth distribution plots
- UMAP visualization colored by sub-HOG (configurable model/level)
- Per-class silhouette scores (original + UMAP space)
- Multi-model comparison (150M/650M/3B) on sub-HOG labels
- Species-per-sub-HOG heatmap

## Interpretation
The HOG hierarchy is too shallow for most proteins to reach level 5+. **Level 3
is the practical sweet spot** for HOG 801468: enough classes (53) with reasonable
balance (CV=0.46) and ~20% coverage. For finer analysis, filtering to depth >= 5
gives 553 proteins with 20 classes at level 5 (CV=0.36).

## Next steps
- Run the notebook to generate UMAP and silhouette results at level 3
- Compare sub-HOG clustering across 150M/650M/3B at level 3
- Investigate whether sub-HOG classes track species boundaries (heatmap)
