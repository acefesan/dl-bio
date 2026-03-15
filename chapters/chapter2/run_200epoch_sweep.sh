#!/bin/bash
# Entry 010: 200-epoch pretraining sweep — all 6 model sizes
# Each model trains for exactly 200 epochs on 10K sequences (9500 train).
# Resumes from checkpoint-50000 of the original runs where available.
#
# Steps per epoch = ceil(9500 / batch_size)
# Total steps = 200 * steps_per_epoch

set -euo pipefail
cd /home/acefsan/src/dl_bio

PYTHON=".venv/bin/python"
SCRIPT="chapters/chapter2/pretrain_esm2.py"
COMMON="--dataset-slice --lr 4e-4 --eval-steps 1000 --save-steps 20000 --save-total-limit 2"

echo "=== 200-epoch pretraining sweep ==="
echo "Start: $(date)"
echo ""

# 1. 6L-320H (7.7M) — batch=16, 200*594=118800 steps
# Resume from entry 005 checkpoint-50000 (already at 50000 steps = 84 epochs)
echo "[1/6] 6L-320H (7.7M) — 118800 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 6 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --num-steps 118800 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260304_133526/checkpoint-50000
echo "[1/6] Done: $(date)"
echo ""

# 2. 12L-320H (15.1M) — batch=16, 118800 steps
echo "[2/6] 12L-320H (15.1M) — 118800 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 12 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --num-steps 118800 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260304_192041/checkpoint-50000
echo "[2/6] Done: $(date)"
echo ""

# 3. 20L-480H (55.9M) — batch=16, 118800 steps
echo "[3/6] 20L-480H (55.9M) — 118800 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 20 --hidden-size 480 --num-heads 20 \
    --batch-size 16 --num-steps 118800 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260307_102803/checkpoint-50000
echo "[3/6] Done: $(date)"
echo ""

# 4. 30L-640H (148.5M) — batch=16, 118800 steps
echo "[4/6] 30L-640H (148.5M) — 118800 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 30 --hidden-size 640 --num-heads 20 \
    --batch-size 16 --num-steps 118800 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260307_112149/checkpoint-50000
echo "[4/6] Done: $(date)"
echo ""

# 5. 33L-1280H (651.7M) — batch=8, 200*1188=237600 steps
echo "[5/6] 33L-1280H (651.7M) — 237600 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 33 --hidden-size 1280 --num-heads 20 \
    --batch-size 8 --num-steps 237600 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260307_130445/checkpoint-50000
echo "[5/6] Done: $(date)"
echo ""

# 6. 33L-1920H (1465M) — batch=4, 200*2375=475000 steps
echo "[6/6] 33L-1920H (1465M) — 475000 steps, resume from 50000"
$PYTHON $SCRIPT $COMMON \
    --num-layers 33 --hidden-size 1920 --num-heads 20 \
    --batch-size 4 --num-steps 475000 \
    --resume chapters/chapter2/runs/pretrain_esm2_20260308_095902/checkpoint-50000
echo "[6/6] Done: $(date)"
echo ""

echo "=== All runs complete: $(date) ==="
