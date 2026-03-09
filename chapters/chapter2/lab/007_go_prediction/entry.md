# 007 — GO Function Prediction from ESM2 Embeddings

**Date:** 2026-03-07
**Models:** ESM2 150M, 650M, 3B
**Status:** complete

## Hypothesis
Entries 001–004 showed inverse scaling for unsupervised clustering (larger ESM2 = worse
silhouette). Does this pattern hold for supervised GO term prediction? If the 3B
embeddings are genuinely less informative, a simple MLP trained on them should also
underperform — ruling out the "geometry artifact" explanation (open question #5).

## Setup
- **Script:** `chapters/chapter2/05_go_prediction.py`
- **Architecture:** Dense(512) → GELU → Dense(256) → GELU → Dense(303), same as `dlfb.proteins.model.Model`
- **Loss:** sigmoid binary cross-entropy (multi-label, 303 GO targets)
- **Train/valid/test:** 3574 / 1191 / 1192 proteins (from preprocessed 150M dataset)
- **Embeddings:** swapped in from each model's feather file (640D, 1280D, 2560D)
- **Training:** 300 steps, batch_size=64, lr=1e-3, Adam, seed=42
- **Metrics:** per-target accuracy, recall, precision, auPRC, auROC (ch2 metrics)

## Results

| Model | Params | Loss | Accuracy | Recall | Precision | auPRC | auROC | Time |
|-------|--------|------|----------|--------|-----------|-------|-------|------|
| 150M | 537K | 0.0742 | 0.9790 | 0.0557 | 0.1381 | 0.2135 | 0.8212 | 11s |
| 650M | 865K | **0.0718** | **0.9796** | **0.0749** | **0.1647** | **0.2389** | **0.8376** | 12s |
| 3B | 1.5M | 0.0934 | 0.9759 | 0.0241 | 0.0720 | 0.0905 | 0.6279 | 10s |

650M best on every metric. 3B dramatically worse — auPRC 0.09 vs 0.24, auROC 0.63 vs 0.84.
3B training curve shows slower loss decay and validation auPRC plateaus at ~0.08.

## Figures
- `figures/go_prediction_comparison.png` — training curves (loss, val loss, val auPRC)
- `figures/go_prediction_test_bars.png` — bar chart of all 6 test metrics by model

## Interpretation
The inverse scaling pattern from unsupervised clustering (entry 004) extends to
supervised prediction. This is strong evidence that the issue is **not** a geometry
artifact — the 3B embeddings are genuinely less informative for this dataset.

Possible explanations:
1. **Anisotropy** (open question #1): 3B embeddings may be concentrated in a narrow
   cone, making linear probes ineffective without whitening
2. **Feature granularity**: 3B may encode features at a scale too fine for 5K proteins
3. **Dimensionality curse**: 2560D with 3574 training samples → ~0.7 samples/dim

The 650M > 150M ordering (unlike clustering where 150M won) suggests that supervised
methods can extract more from richer representations, while unsupervised methods suffer
from the curse of dimensionality.

## Next steps
- PCA/whitening on 3B embeddings before training to test anisotropy hypothesis
- Increase training steps (currently 300) to see if 3B eventually catches up
- Try larger hidden_dim for the 3B model to give it more capacity
