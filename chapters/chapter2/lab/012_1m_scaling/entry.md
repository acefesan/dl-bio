# 012 — 1M-Sequence Scaling Experiment

**Date:** 2026-03-15
**Models:** ESM2 from scratch — 5 sizes (7.7M to 651.7M)
**Status:** in-progress

## Hypothesis
Entries 005/008/010 trained on 10K UniRef sequences and found all models converge to
~2.56 eval loss (perplexity ~12.9) regardless of size. The bottleneck was data, not
capacity. With 1M sequences (100x more data), larger models should separate from
smaller ones, approaching the ESM2 paper's perplexities (8M: 10.45, 150M: 7.75,
650M: 6.95).

## Setup

**Data:** `nvidia/esm2_uniref_pretraining_data` full split, limited to first 1M
sequences via `--num-sequences 1000000`. 95/5 train/eval split (950K/50K).

**Models:** All trained from random initialization with identical hyperparameters
(except batch size for GPU memory).

| Name | Layers | Hidden | Heads | Params | Batch | LR | Est. VRAM |
|------|--------|--------|-------|--------|-------|----|-----------|
| 6L-320H | 6 | 320 | 20 | 7.7M | 16 | 4e-4 | 0.2 GB |
| 12L-320H | 12 | 320 | 20 | 15.1M | 16 | 4e-4 | 0.3 GB |
| 20L-480H | 20 | 480 | 20 | 55.9M | 16 | 4e-4 | 1.1 GB |
| 30L-640H | 30 | 640 | 20 | 148.5M | 16 | 4e-4 | 3.0 GB |
| 33L-1280H | 33 | 1280 | 20 | 651.7M | 8 | 2e-4 | 13.0 GB |

**Key fixes from entry 010:**
- Gradient clipping (`max_grad_norm=1.0`) to prevent divergence
- 650M model gets halved LR (2e-4) to compensate for batch=8 (higher gradient noise)
- All runs from scratch (no resume bug)

**Training:** 50,000 steps each, eval every 2500, save every 10000.
- 1M seqs at batch 16 → 62,500 steps/epoch → 50K steps ≈ 0.8 epochs
- This means each sequence is seen less than once on average — no memorization risk

**Hardware:** RTX 5090 32GB, bf16

## Runtime Estimates

| Model | Steps/sec (est) | Total time (est) |
|-------|----------------|-----------------|
| 6L-320H | ~50 | ~17 min |
| 12L-320H | ~35 | ~24 min |
| 20L-480H | ~15 | ~55 min |
| 30L-640H | ~8 | ~1.7 hrs |
| 33L-1280H | ~3 | ~4.6 hrs |
| **Total** | | **~7.5 hrs** |

## Commands

```bash
# Run from project root with venv activated
cd /home/acefsan/src/dl_bio

# 6L-320H (7.7M)
.venv/bin/python chapters/chapter2/04_pretrain_esm2.py \
    --num-sequences 1000000 --num-layers 6 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --lr 4e-4 --num-steps 50000 --eval-steps 2500 --save-steps 10000

# 12L-320H (15.1M)
.venv/bin/python chapters/chapter2/04_pretrain_esm2.py \
    --num-sequences 1000000 --num-layers 12 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --lr 4e-4 --num-steps 50000 --eval-steps 2500 --save-steps 10000

# 20L-480H (55.9M)
.venv/bin/python chapters/chapter2/04_pretrain_esm2.py \
    --num-sequences 1000000 --num-layers 20 --hidden-size 480 --num-heads 20 \
    --batch-size 16 --lr 4e-4 --num-steps 50000 --eval-steps 2500 --save-steps 10000

# 30L-640H (148.5M)
.venv/bin/python chapters/chapter2/04_pretrain_esm2.py \
    --num-sequences 1000000 --num-layers 30 --hidden-size 640 --num-heads 20 \
    --batch-size 16 --lr 4e-4 --num-steps 50000 --eval-steps 2500 --save-steps 10000

# 33L-1280H (651.7M) — halved LR for smaller batch
.venv/bin/python chapters/chapter2/04_pretrain_esm2.py \
    --num-sequences 1000000 --num-layers 33 --hidden-size 1280 --num-heads 20 \
    --batch-size 8 --lr 2e-4 --num-steps 50000 --eval-steps 2500 --save-steps 10000
```

## Expected Outcomes

1. **Larger models should achieve lower loss** with 1M sequences (unlike 10K where all converged)
2. **Gap should narrow at 50K steps** — larger models learn slower per step but extract more per sample
3. **No divergence** with gradient clipping + scaled LR
4. Comparison to ESM2 paper (trained on 187M seqs, 500K steps):
   - Our 8M model: expect perplexity ~11-12 (paper: 10.45)
   - Our 150M model: expect perplexity ~9-10 (paper: 7.75)
   - Our 650M model: expect perplexity ~8-10 (paper: 6.95)

## Results

*(to be filled after runs complete)*

## Next steps

- If models separate: extend to 200K+ steps or full dataset
- Compare learning curves (loss vs step) across scales
- Spectral analysis on checkpoints to compare rank dynamics with more data
