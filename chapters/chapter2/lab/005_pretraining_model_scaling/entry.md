# 005 — Pretraining Model Scaling (7.7M → 650M on 10K sequences)

**Date:** 2026-03-07
**Models:** ESM2 from scratch — 5 sizes (7.7M, 15.1M, 55.9M, 148.5M, 651.7M)
**Status:** complete

## Hypothesis
Training progressively larger ESM2 models from scratch on a fixed 10K-sequence
dataset should eventually cause overfitting — the model has enough capacity to
memorize the training set, causing eval loss to diverge from train loss. We also
track weight spectral dynamics (effective rank, weight norms) to see how
implicit regularization changes with scale.

## Setup
- **Dataset:** `nvidia/esm2_uniref_pretraining_data`, first 10K sequences via streaming
- **Objective:** Masked language modeling (15% random masking)
- **Split:** 9,500 train / 500 eval (95/5)
- **Tokenizer:** `facebook/esm2_t6_8M_UR50D` (shared ESM2 vocab, 33 tokens)
- **Training:** 50,000 steps, lr=4e-4, weight_decay=0.01, warmup=500 steps, bf16
- **Script:** `chapters/chapter2/pretrain_esm2.py`
- **Analysis:** `chapters/chapter2/analyze_training_dynamics.py`

| Label | Layers | Hidden | Heads | Params | Batch | Epochs |
|-------|--------|--------|-------|--------|-------|--------|
| 6L-320H | 6 | 320 | 20 | 7.7M | 16 | 84.2 |
| 12L-320H | 12 | 320 | 20 | 15.1M | 16 | 84.2 |
| 20L-480H | 20 | 480 | 20 | 55.9M | 16 | 84.2 |
| 30L-640H | 30 | 640 | 20 | 148.5M | 16 | 84.2 |
| 33L-1280H | 33 | 1280 | 20 | 651.7M | 8 | 42.1 |

Note: 650M model used batch=8 (GPU memory), so it saw half the epochs in 50K steps.

## Results

### Loss metrics

| Model | Params | Final train | Final eval | Best eval | @Step | Gen gap |
|-------|--------|-------------|------------|-----------|-------|---------|
| 6L-320H | 7.7M | 2.5727 | 2.5879 | 2.5625 | 41500 | +0.015 |
| 12L-320H | 15.1M | 2.5741 | 2.5888 | 2.5611 | 41500 | +0.015 |
| 20L-480H | 55.9M | 2.5735 | 2.5884 | 2.5625 | 41500 | +0.015 |
| 30L-640H | 148.5M | 2.5818 | 2.5943 | 2.5712 | 41500 | +0.013 |
| 33L-1280H | 651.7M | 2.6573 | 2.6602 | 2.6275 | 2500 | +0.003 |

### Key observations

1. **No overfitting at any scale.** Generalization gap stays at +0.01–0.02 for all
   models. Even 650M params on 10K sequences doesn't overfit with MLM + weight decay.

2. **Inverse scaling for loss.** The 650M model has the *worst* eval loss (2.63 vs
   2.56 for smaller models). The 148M also slightly underperforms (2.57). Models
   7.7M–55.9M are essentially identical (2.56).

3. **650M best eval at step 2500 then degrades.** Unlike smaller models (all best at
   step 41500), the 650M peaks very early and gets worse — but note it only trained
   42 epochs (batch=8) vs 84 for others.

4. **All models best at same step.** The 7.7M through 148M models all peak at step
   41500. The learning dynamics are synchronized despite 20x parameter difference.

### Weight spectral dynamics

5. **Weight norm grows with model size.** Total Frobenius norm at end of training:
   ~110 (6L) → ~220 (12L) → ~250 (20L) → ~350 (30L) → ~430 (33L).

6. **Effective rank decline is steeper in deeper models.** Normalized effective rank
   of attention layers: 6L drops to ~0.63, 12L to ~0.60, 20L to ~0.58, 30L to ~0.52,
   33L to ~0.48. Confirms Huh et al.'s low-rank simplicity bias prediction.

7. **Stable rank heatmap shows layer specialization.** In all models, early layers
   stay low-rank while the final layer maintains highest stable rank. This gradient
   is more pronounced in deeper models.

## Figures
- `figures/scaling_dynamics_comparison.png` — all 5 models: loss, gap, grad norm,
  weight norms, effective rank, stable rank heatmaps
- `figures/6L_vs_12L_dynamics.png` — detailed 6L vs 12L comparison (earlier analysis)

## Interpretation

**Why no overfitting?** Three factors conspire:
1. Stochastic MLM masking generates different targets each epoch, so the effective
   dataset size is much larger than 10K
2. Weight decay (0.01) provides explicit regularization
3. SGD + weight decay implicitly minimizes rank (Galanti et al.), compressing excess
   capacity rather than memorizing

**Why does the 650M model perform worse?** Two likely factors:
1. With batch=8 it only saw 42 epochs vs 84 — half the data exposure. A fair
   comparison would need ~100K steps.
2. Larger models need more data to reach their potential. 10K sequences is not enough
   signal for 650M parameters — the model can't find useful features beyond what
   7.7M already captures.

**Connection to entry 004 (inverse scaling of pretrained ESM2 embeddings):**
Entry 004 showed pretrained ESM2 embeddings get *worse* for clustering as model
size increases (150M > 650M > 3B). Here we see a similar pattern from scratch:
bigger models don't benefit from small data. The mechanisms may be related —
larger models may learn representations that are powerful but geometrically
unfriendly to simple distance metrics.

## Next steps
- Run 650M for 100K steps (batch=8) to give it fair epoch count
- Remove weight decay to see if overfitting finally appears
- Try even smaller datasets (1K, 500 sequences) to force memorization
- Dataset cartography: track per-sequence confidence across checkpoints
- Compare spectral dynamics to pretrained checkpoints (facebook/esm2_t6_8M_UR50D)
