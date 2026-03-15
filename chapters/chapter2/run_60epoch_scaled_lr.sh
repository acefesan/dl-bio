#!/bin/bash
# Entry 011: 60-epoch pretraining with scaled learning rates
# Scale lr proportionally with batch size to prevent divergence
# No resume — train from scratch for fair comparison
#
# LR scaling: lr_new = lr_base * (batch / batch_base)
# where batch_base = 16, lr_base = 4e-4

set -euo pipefail
cd /home/acefsan/src/dl_bio

PYTHON=".venv/bin/python"
SCRIPT="chapters/chapter2/pretrain_esm2.py"
COMMON="--dataset-slice --eval-steps 500 --save-steps 10000 --save-total-limit 2"

echo "=== 60-epoch pretraining with scaled learning rates ==="
echo "Start: $(date)"
echo ""

# 1. 6L-320H — batch=16, 60*594=35640 steps, lr=4e-4
echo "[1/6] 6L-320H (7.7M) — 35640 steps, lr=4e-4"
$PYTHON $SCRIPT $COMMON \
    --num-layers 6 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --num-steps 35640 --lr 4e-4
echo "[1/6] Done: $(date)"
echo ""

# 2. 12L-320H — batch=16, 35640 steps, lr=4e-4
echo "[2/6] 12L-320H (15.1M) — 35640 steps, lr=4e-4"
$PYTHON $SCRIPT $COMMON \
    --num-layers 12 --hidden-size 320 --num-heads 20 \
    --batch-size 16 --num-steps 35640 --lr 4e-4
echo "[2/6] Done: $(date)"
echo ""

# 3. 20L-480H — batch=16, 35640 steps, lr=4e-4
echo "[3/6] 20L-480H (55.9M) — 35640 steps, lr=4e-4"
$PYTHON $SCRIPT $COMMON \
    --num-layers 20 --hidden-size 480 --num-heads 20 \
    --batch-size 16 --num-steps 35640 --lr 4e-4
echo "[3/6] Done: $(date)"
echo ""

# 4. 30L-640H — batch=16, 35640 steps, lr=4e-4
echo "[4/6] 30L-640H (148.5M) — 35640 steps, lr=4e-4"
$PYTHON $SCRIPT $COMMON \
    --num-layers 30 --hidden-size 640 --num-heads 20 \
    --batch-size 16 --num-steps 35640 --lr 4e-4
echo "[4/6] Done: $(date)"
echo ""

# 5. 33L-1280H (650M) — batch=8, 60*1188=71280 steps, lr=2e-4 (scaled by 8/16)
echo "[5/6] 33L-1280H (651.7M) — 71280 steps, lr=2e-4 (scaled)"
$PYTHON $SCRIPT $COMMON \
    --num-layers 33 --hidden-size 1280 --num-heads 20 \
    --batch-size 8 --num-steps 71280 --lr 2e-4
echo "[5/6] Done: $(date)"
echo ""

# 6. 33L-1920H (1.47B) — batch=4, 60*2375=142500 steps, lr=1e-4 (scaled by 4/16)
echo "[6/6] 33L-1920H (1465M) — 142500 steps, lr=1e-4 (scaled)"
$PYTHON $SCRIPT $COMMON \
    --num-layers 33 --hidden-size 1920 --num-heads 20 \
    --batch-size 4 --num-steps 142500 --lr 1e-4
echo "[6/6] Done: $(date)"
echo ""

echo "=== All runs complete: $(date) ==="
