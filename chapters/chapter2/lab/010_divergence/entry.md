# 010 — 200-Epoch Sweep: Divergence in Large Models

**Date:** 2026-03-09
**Models:** ESM2 from scratch — 6 sizes (7.7M to 1,465M)
**Status:** partial (4/6 models complete, 650M diverged, 1.47B killed)

## Hypothesis
Entries 005/008 had an epoch confound: smaller models trained 84 epochs (batch=16)
while 650M trained 42 (batch=8) and 1.47B trained 21 (batch=4). Running all models
for exactly 200 epochs eliminates this confound.

## Setup
- Resumed from checkpoint-50000 of each original run
- Target: 200 epochs (steps adjusted per batch size)
- eval_steps=1000, save_steps=20000, save_total_limit=2
- **Bug:** `--resume` loaded weights but reset the step counter and lr schedule

| Model | Batch | Target steps | Resumed from | Actual behavior |
|-------|-------|-------------|--------------|-----------------|
| 6L-320H | 16 | 118,800 | step 50000 | 200 fresh epochs from warm weights |
| 12L-320H | 16 | 118,800 | step 50000 | same |
| 20L-480H | 16 | 118,800 | step 50000 | same |
| 30L-640H | 16 | 118,800 | step 50000 | same |
| 33L-1280H | 8 | 237,600 | step 50000 | **DIVERGED at epoch 44** |
| 33L-1920H | 4 | 475,000 | step 50000 | killed (expected same divergence) |

## Results

### Stable models (batch=16, effective LR = 2.5e-5)

| Model | Params | Final train | Final eval | Status |
|-------|--------|-------------|------------|--------|
| 6L-320H | 7.7M | 1.4921 | 2.5862 | converged |
| 12L-320H | 15.1M | 1.4923 | 2.5832 | converged |
| 20L-480H | 55.9M | 1.4909 | 2.5853 | converged |
| 30L-640H | 148.5M | 1.4949 | 2.5879 | converged |

### Divergent model (batch=8, effective LR = 5.0e-5)

| Model | Params | Eval @ ep 0 | Eval @ ep 44 | Eval @ ep 200 |
|-------|--------|-------------|--------------|---------------|
| 33L-1280H | 651.7M | 2.6819 | 2.8760 | **4.3819** |

Divergence began at epoch ~44 (step 52000). Loss went from 2.68 to >4.3 — worse
than random initialization (log(33) = 3.50).

## Root Cause

1. **Resume bug:** `--resume` loaded checkpoint weights but created a new output dir,
   resetting the step counter to 0. The lr schedule restarted (warmup + peak lr).
2. **Double effective LR:** batch=8 gives eta_eff = 4e-4/8 = 5e-5, vs 2.5e-5 for
   batch=16 models. The 650M model received 2x the effective learning rate.
3. **Sharper loss landscape:** the 650M model has lower effective rank (0.48 vs 0.63),
   meaning sharper curvature. Peak lr=4e-4 exceeded the stability threshold.

## Figures
- `figures/divergence_figure.png` — stable vs divergent loss curves + lr schedule
- `figures/divergence_conceptual.png` — loss landscape + gradient noise diagrams
- `divergence_article.pdf` — detailed explainer article

## Interpretation
The divergence is a **learning rate / batch size / model size interaction**. Smaller
batches produce noisier gradients, and larger models have sharper loss landscapes.
The combination makes large models with small batches inherently unstable at the
same nominal learning rate.

## Next steps
- Rerun all 6 models from scratch (no resume) for 200 epochs
- Scale lr with batch size: lr = 4e-4 * (batch/16) → 2e-4 for 650M, 1e-4 for 1.47B
- Or use a single lr that's safe for all: 1e-4 with gradient clipping
