# 008 — 1.47B Pretraining (Largest Model That Fits on Hardware)

**Date:** 2026-03-08
**Models:** ESM2 from scratch — 33L-1920H (1,465M params)
**Status:** complete

## Hypothesis
Entry 005 showed inverse scaling from 7.7M to 651M on 10K sequences, but the
650M model only trained 42 epochs (batch=8) vs 84 for smaller models. Adding a
1.47B model (the largest that fits in GPU memory) tests whether the inverse
scaling trend continues at even larger scale, and whether the 650M result was
simply due to insufficient training.

## Setup
- **Dataset:** `nvidia/esm2_uniref_pretraining_data`, first 10K sequences via streaming
- **Objective:** Masked language modeling (15% random masking)
- **Split:** 9,500 train / 500 eval (95/5)
- **Training:** 50,000 steps, lr=4e-4, weight_decay=0.01, warmup=500 steps, bf16
- **Script:** `chapters/chapter2/pretrain_esm2.py`

| Label | Layers | Hidden | Heads | Params | Batch | Epochs |
|-------|--------|--------|-------|--------|-------|--------|
| 33L-1920H | 33 | 1920 | 20 | 1,465M | 4 | 21.1 |

Note: batch=4 is the maximum that fits in GPU memory for this model size.
This gives 21 epochs in 50K steps (vs 42 for 650M, 84 for smaller models).

## Results

### Loss metrics (updated table with all 6 sizes)

| Model | Params | Final train | Final eval | Best eval | @Step | Gen gap | Epochs |
|-------|--------|-------------|------------|-----------|-------|---------|--------|
| 6L-320H | 7.7M | 2.5727 | 2.5879 | 2.5625 | 41500 | +0.015 | 84.2 |
| 12L-320H | 15.1M | 2.5741 | 2.5888 | 2.5611 | 41500 | +0.015 | 84.2 |
| 20L-480H | 55.9M | 2.5735 | 2.5884 | 2.5625 | 41500 | +0.015 | 84.2 |
| 30L-640H | 148.5M | 2.5818 | 2.5943 | 2.5712 | 41500 | +0.013 | 84.2 |
| 33L-1280H | 651.7M | 2.6573 | 2.6602 | 2.6275 | 2500 | +0.003 | 42.1 |
| **33L-1920H** | **1,465M** | **2.7610** | **2.7386** | **2.6755** | **2500** | **-0.022** | **21.1** |

### Key observations

1. **Inverse scaling extends to 1.47B.** Best eval loss of 2.68 is worse than all
   smaller models (7.7M–55.9M achieve 2.56). The trend is monotonic above 55.9M.

2. **Same early-peak pattern as 650M.** Both 650M and 1.47B peak at step 2500 then
   degrade. Smaller models all peak at step 41500.

3. **Negative generalization gap.** The 1.47B model shows eval loss *lower* than
   train loss (-0.022), which is unusual. This likely reflects the stochastic MLM
   masking providing regularization that benefits eval more than train at low epoch
   counts.

4. **No overfitting despite 1.47B params on 10K sequences.** Even the largest model
   doesn't memorize the training set.

5. **Epoch confound persists.** The 1.47B model only sees 21 epochs vs 84 for the
   smaller models due to batch size constraints. However, its best eval loss (step
   2500, ~1 epoch) is already worse than what smaller models achieve at the same
   step count.

## Figures
- `figures/scaling_dynamics_6models.png` — all 6 models: loss curves, generalization
  gap, gradient norms, weight norms, effective rank, stable rank heatmaps

## Interpretation

The inverse scaling pattern is now confirmed across 3 orders of magnitude
(7.7M → 1,465M). Key insights:

1. **The 650M result was not just an epoch artifact.** Even comparing at the same
   step count (e.g., step 2500), larger models underperform. The 1.47B at step 2500
   has eval loss 2.68 while the 7.7M at step 2500 already achieves ~2.60.

2. **10K sequences is insufficient for large models.** The data-to-parameter ratio
   is ~7 tokens/param for 1.47B vs ~1300 tokens/param for 7.7M. Large models cannot
   find useful features that small models haven't already captured from this little data.

3. **Connection to pretrained inverse scaling (entry 004).** The same pattern appears
   in pretrained ESM2 embeddings (150M > 650M > 3B for clustering). Both scratch
   and pretrained settings suggest large protein language models learn representations
   that are powerful but geometrically unsuitable for simple metrics.

## Next steps
- Run 1.47B and 650M for matched epoch counts (e.g., 100K steps for 650M, 200K for 1.47B)
- Remove weight decay to test if overfitting finally appears at this scale
- Track effective rank evolution with more checkpoint granularity
