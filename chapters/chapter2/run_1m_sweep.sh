#!/bin/bash
set -e
cd /home/acefsan/src/dl_bio
PYTHON=".venv/bin/python"
TRAIN="chapters/chapter2/04_pretrain_esm2.py"
COMMON="--num-sequences 1000000 --num-steps 50000 --eval-steps 2500 --save-steps 10000"

echo "=== 1M Scaling Sweep — $(date) ==="

echo -e "\n>>> 6L-320H (7.7M) — $(date)"
$PYTHON $TRAIN $COMMON --num-layers 6 --hidden-size 320 --num-heads 20 --batch-size 16 --lr 4e-4

echo -e "\n>>> 12L-320H (15.1M) — $(date)"
$PYTHON $TRAIN $COMMON --num-layers 12 --hidden-size 320 --num-heads 20 --batch-size 16 --lr 4e-4

echo -e "\n>>> 20L-480H (55.9M) — $(date)"
$PYTHON $TRAIN $COMMON --num-layers 20 --hidden-size 480 --num-heads 20 --batch-size 16 --lr 4e-4

echo -e "\n>>> 30L-640H (148.5M) — $(date)"
$PYTHON $TRAIN $COMMON --num-layers 30 --hidden-size 640 --num-heads 20 --batch-size 16 --lr 4e-4

echo -e "\n>>> 33L-1280H (651.7M) — $(date)"
$PYTHON $TRAIN $COMMON --num-layers 33 --hidden-size 1280 --num-heads 20 --batch-size 8 --lr 2e-4

echo -e "\n=== Sweep complete — $(date) ==="
