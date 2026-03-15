#!/usr/bin/env python2
"""Pretrain ESM2 from scratch using masked language modeling on UniRef data.

Trains a randomly initialized ESM2 model (8M param variant by default)
using the standard MLM objective: mask 15% of tokens and predict them.

Uses NVIDIA's ESM2 UniRef pretraining data from HuggingFace Hub.

Usage:
    # Quick smoke test (~2 min on GPU)
    python chapters/chapter2/04_pretrain_esm2.py --dataset-slice --num-steps 50 --batch-size 4

    # Toy run with defaults (10k sequences, 5000 steps)
    python chapters/chapter2/04_pretrain_esm2.py --dataset-slice

    # Full dataset (requires significant compute)
    python chapters/chapter2/04_pretrain_esm2.py
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

import torch
from datasets import load_dataset
from transformers import (
    DataCollatorForLanguageModeling,
    EsmConfig,
    EsmForMaskedLM,
    EsmTokenizer,
    Trainer,
    TrainingArguments,
    set_seed,
)

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHAPTER_DIR = Path(__file__).parent

# Tokenizer from the 8M param ESM2 variant (shared vocabulary across all ESM2 sizes)
TOKENIZER_CHECKPOINT = "facebook/esm2_t6_8M_UR50D"

# NVIDIA's pretraining data on HuggingFace Hub
DATASET_NAME = "nvidia/esm2_uniref_pretraining_data"
DATASET_SLICE = "UniRef50_10K_cluster"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pretrain ESM2 from scratch with masked language modeling"
    )

    # Model architecture
    parser.add_argument("--num-layers", type=int, default=6,
                        help="Number of transformer layers (default: 6)")
    parser.add_argument("--hidden-size", type=int, default=320,
                        help="Hidden dimension size (default: 320)")
    parser.add_argument("--num-heads", type=int, default=20,
                        help="Number of attention heads (default: 20)")

    # Data
    parser.add_argument("--dataset-slice", action="store_true",
                        help="Use 10K subset instead of full UniRef dataset")
    parser.add_argument("--max-length", type=int, default=512,
                        help="Max sequence length for tokenizer (default: 512)")

    # Training
    parser.add_argument("--batch-size", type=int, default=16,
                        help="Per-device training batch size (default: 16)")
    parser.add_argument("--num-steps", type=int, default=5000,
                        help="Total training steps (default: 5000)")
    parser.add_argument("--lr", type=float, default=4e-4,
                        help="Learning rate (default: 4e-4)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")

    # Checkpointing / logging
    parser.add_argument("--eval-steps", type=int, default=500,
                        help="Evaluate every N steps (default: 500)")
    parser.add_argument("--save-steps", type=int, default=1000,
                        help="Save checkpoint every N steps (default: 1000)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: runs/pretrain_esm2_{timestamp})")
    parser.add_argument("--resume", type=str, default=None,
                        help="Resume training from checkpoint directory")

    return parser.parse_args()


def tokenize_sequences(examples, tokenizer, max_length):
    """Tokenize protein sequences, truncating to max_length."""
    return tokenizer(
        examples["sequence"],
        truncation=True,
        max_length=max_length,
        padding=False,
        return_special_tokens_mask=True,
    )


def main():
    args = parse_args()
    set_seed(args.seed)

    # Output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = CHAPTER_DIR / "runs" / f"pretrain_esm2_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"Device: {'cuda' if torch.cuda.is_available() else 'cpu'}")

    # --- Tokenizer ---
    print(f"\nLoading tokenizer from {TOKENIZER_CHECKPOINT}...")
    tokenizer = EsmTokenizer.from_pretrained(TOKENIZER_CHECKPOINT)

    # --- Model (random weights) ---
    print(f"\nInitializing ESM2 from scratch: {args.num_layers} layers, "
          f"{args.hidden_size} hidden, {args.num_heads} heads")
    config = EsmConfig(
        vocab_size=tokenizer.vocab_size,
        num_hidden_layers=args.num_layers,
        hidden_size=args.hidden_size,
        num_attention_heads=args.num_heads,
        intermediate_size=args.hidden_size * 4,
        max_position_embeddings=args.max_length + 2,  # +2 for special tokens
    )
    model = EsmForMaskedLM(config)
    num_params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {num_params:,} ({num_params / 1e6:.1f}M)")

    # --- Dataset ---
    subset = DATASET_SLICE if args.dataset_slice else None
    print(f"\nLoading dataset: {DATASET_NAME}" +
          (f" (slice: {subset})" if subset else " (full)"))
    dataset = load_dataset(DATASET_NAME, subset, trust_remote_code=True)

    # The dataset may have a single split — use whatever is available
    if "train" in dataset:
        raw_dataset = dataset["train"]
    else:
        split_name = list(dataset.keys())[0]
        uaw_dataset = dataset[split_name]
    print(f"Sequences loaded: {len(raw_dataset):,}")

    # Tokenize
    print("Tokenizing sequences...")
    tokenized = raw_dataset.map(
        lambda ex: tokenize_sequences(ex, tokenizer, args.max_length),
        batched=True,
        remove_columns=raw_dataset.column_names,
        desc="Tokenizing",
    )

    # Split into train/eval (95/5)
    split = tokenized.train_test_split(test_size=0.05, seed=args.seed)
    train_dataset = split["train"]
    eval_dataset = split["test"]
    print(f"Train: {len(train_dataset):,}  |  Eval: {len(eval_dataset):,}")

    # --- MLM data collator ---
    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer,
        mlm=True,
        mlm_probability=0.15,
    )

    # --- Training ---
    training_args = TrainingArguments(
        output_dir=str(output_dir),
        max_steps=args.num_steps,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        learning_rate=args.lr,
        weight_decay=0.01,
        warmup_steps=min(500, args.num_steps // 10),
        eval_strategy="steps",
        eval_steps=args.eval_steps,
        save_strategy="steps",
        save_steps=args.save_steps,
        logging_steps=50,
        seed=args.seed,
        bf16=torch.cuda.is_available(),
        dataloader_num_workers=2,
        report_to="none",
        save_total_limit=3,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset,
        data_collator=data_collator,
    )

    # Train
    print(f"\nStarting training for {args.num_steps} steps...")
    checkpoint = args.resume if args.resume else None
    train_result = trainer.train(resume_from_checkpoint=checkpoint)

    # Save final model + metrics
    trainer.save_model(str(output_dir / "final_model"))
    metrics = train_result.metrics
    trainer.log_metrics("train", metrics)
    trainer.save_metrics("train", metrics)

    # Evaluate
    print("\nRunning final evaluation...")
    eval_metrics = trainer.evaluate()
    trainer.log_metrics("eval", eval_metrics)
    trainer.save_metrics("eval", eval_metrics)

    print(f"\nTraining complete!")
    print(f"  Final train loss: {metrics['train_loss']:.4f}")
    print(f"  Final eval loss:  {eval_metrics['eval_loss']:.4f}")
    print(f"  Model saved to:   {output_dir / 'final_model'}")


if __name__ == "__main__":
