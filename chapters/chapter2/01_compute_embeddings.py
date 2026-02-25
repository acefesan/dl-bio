#!/usr/bin/env python3
"""Compute protein embeddings for all species using ESM2.

Reads protein sequences from a FASTA file, computes mean-pooled embeddings
using the ESM2 protein language model, and saves the results as a Feather file.
Supports batch checkpointing and resume for long-running computations.

Usage:
    python 01_compute_embeddings.py --limit 100  # Test run with 100 sequences
    python 01_compute_embeddings.py              # Full run (all sequences)
    python 01_compute_embeddings.py --resume     # Resume from partial results
    python 01_compute_embeddings.py --max-seq-length 4000  # Skip very long sequences
"""

import argparse
import time
from dataclasses import dataclass
from math import ceil
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from transformers import AutoTokenizer, EsmModel

# Paths relative to this script's location (chapters/chapter2/)
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATASETS_PATH = PROJECT_ROOT / "assets/proteins/datasets"

FASTA_FILE = DATASETS_PATH / "train_sequences.fasta"
DEFAULT_OUTPUT = DATASETS_PATH / "all_species_embeddings.feather"
DEFAULT_CHECKPOINT_DIR = PROJECT_ROOT / "checkpoints"

MODEL_CHECKPOINT = "facebook/esm2_t33_650M_UR50D"


@dataclass
class TimingStats:
    """Track timing statistics for embedding computation."""
    sequence_lengths: list[int]
    batch_times: list[float]
    batch_sizes: list[int]

    def __init__(self):
        self.sequence_lengths = []
        self.batch_times = []
        self.batch_sizes = []

    def record_batch(self, lengths: list[int], elapsed: float):
        self.sequence_lengths.extend(lengths)
        self.batch_times.append(elapsed)
        self.batch_sizes.append(len(lengths))

    def summary(self) -> dict:
        total_seqs = len(self.sequence_lengths)
        total_time = sum(self.batch_times)
        avg_per_seq = total_time / total_seqs if total_seqs > 0 else 0

        return {
            "total_sequences": total_seqs,
            "total_time_sec": total_time,
            "avg_time_per_seq_ms": avg_per_seq * 1000,
            "avg_seq_length": np.mean(self.sequence_lengths) if self.sequence_lengths else 0,
            "min_seq_length": min(self.sequence_lengths) if self.sequence_lengths else 0,
            "max_seq_length": max(self.sequence_lengths) if self.sequence_lengths else 0,
            "sequences_per_second": total_seqs / total_time if total_time > 0 else 0,
        }

    def print_summary(self):
        stats = self.summary()
        print("\n" + "=" * 60)
        print("TIMING SUMMARY")
        print("=" * 60)
        print(f"Total sequences processed: {stats['total_sequences']}")
        print(f"Total time: {stats['total_time_sec']:.2f} sec")
        print(f"Avg time per sequence: {stats['avg_time_per_seq_ms']:.2f} ms")
        print(f"Throughput: {stats['sequences_per_second']:.2f} seq/sec")
        print(f"Sequence lengths: min={stats['min_seq_length']}, "
              f"avg={stats['avg_seq_length']:.1f}, max={stats['max_seq_length']}")
        print("=" * 60)


def get_device() -> torch.device:
    """Return CUDA device if available, otherwise CPU."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def get_mean_embeddings(
    sequences: list[str],
    tokenizer,
    model,
    device: torch.device,
) -> np.ndarray:
    """Compute mean embedding for each sequence using a protein LM."""
    model_inputs = tokenizer(sequences, padding=True, return_tensors="pt")
    model_inputs = {k: v.to(device) for k, v in model_inputs.items()}

    model = model.to(device)
    model.eval()

    with torch.no_grad():
        outputs = model(**model_inputs)
        mean_embeddings = outputs.last_hidden_state.mean(dim=1)

    return mean_embeddings.detach().cpu().numpy()


def load_sequences(fasta_path: Path, limit: int | None = None) -> pd.DataFrame:
    """Load protein sequences from FASTA file."""
    print(f"Loading sequences from {fasta_path}...")

    data = []
    fasta_sequences = SeqIO.parse(open(fasta_path), "fasta")

    for i, fasta in enumerate(fasta_sequences):
        if limit is not None and i >= limit:
            break
        data.append({
            "EntryID": fasta.id,
            "Sequence": str(fasta.seq),
            "Length": len(fasta.seq),
        })

    df = pd.DataFrame(data)
    print(f"Loaded {len(df)} sequences")
    return df


def compute_embeddings_with_timing(
    sequence_df: pd.DataFrame,
    tokenizer,
    model,
    batch_size: int = 32,
    checkpoint_dir: Path | None = None,
    max_seq_length: int | None = None,
) -> tuple[np.ndarray, TimingStats, np.ndarray, pd.DataFrame | None]:
    """Compute embeddings with incremental saving.

    Returns embeddings, stats, original_indices, and skipped_df.
    """
    device = get_device()
    print(f"Using device: {device}")

    if device.type == "cuda":
        print(f"GPU: {torch.cuda.get_device_name()}")
        print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")

    # Sort by length to minimize padding overhead within batches
    original_indices = sequence_df.index.to_numpy()
    sequence_df = sequence_df.sort_values("Length").reset_index(drop=True)
    print(f"Sorted sequences by length (min={sequence_df['Length'].min()}, max={sequence_df['Length'].max()})")

    # Filter out sequences that are too long
    skipped_df = None
    if max_seq_length is not None:
        mask = sequence_df["Length"] <= max_seq_length
        skipped_df = sequence_df[~mask].copy()
        sequence_df = sequence_df[mask].reset_index(drop=True)
        if len(skipped_df) > 0:
            print(f"Skipping {len(skipped_df)} sequences longer than {max_seq_length} aa")
            original_indices = original_indices[mask.values if hasattr(mask, 'values') else mask]

    # Save sequence metadata for resume
    if checkpoint_dir is not None:
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        meta_path = checkpoint_dir / "metadata.npz"
        np.savez(meta_path,
                 entry_ids=sequence_df["EntryID"].values,
                 lengths=sequence_df["Length"].values,
                 original_indices=original_indices)

    model = model.to(device)
    model.eval()

    n_batches = ceil(len(sequence_df) / batch_size)
    stats = TimingStats()
    all_embeddings = []

    # Check which batches already exist (for resume)
    start_batch = 0
    if checkpoint_dir is not None:
        existing = sorted(checkpoint_dir.glob("batch_*.npy"))
        if existing:
            for f in existing:
                all_embeddings.append(np.load(f))
            start_batch = len(existing)
            n_existing_seqs = sum(e.shape[0] for e in all_embeddings)
            print(f"Found {start_batch} existing batches ({n_existing_seqs} sequences), resuming...")

    print(f"\nProcessing {len(sequence_df)} sequences in {n_batches} batches (batch_size={batch_size})...")

    for i in range(start_batch, n_batches):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, len(sequence_df))

        batch_df = sequence_df.iloc[start_idx:end_idx]
        batch_seqs = list(batch_df["Sequence"])
        batch_lengths = list(batch_df["Length"])

        if device.type == "cuda":
            torch.cuda.synchronize()

        start_time = time.perf_counter()

        embeddings = get_mean_embeddings(batch_seqs, tokenizer, model, device)

        if device.type == "cuda":
            torch.cuda.synchronize()

        elapsed = time.perf_counter() - start_time

        stats.record_batch(batch_lengths, elapsed)
        all_embeddings.append(embeddings)

        # Save batch immediately
        if checkpoint_dir is not None:
            batch_path = checkpoint_dir / f"batch_{i:05d}.npy"
            np.save(batch_path, embeddings)

        # Progress update every 10 batches
        if (i + 1) % 10 == 0 or i == n_batches - 1:
            current_stats = stats.summary()
            max_len_in_batch = max(batch_lengths)
            print(f"  Batch {i+1}/{n_batches} | "
                  f"MaxLen: {max_len_in_batch} | "
                  f"Elapsed: {current_stats['total_time_sec']:.1f}s | "
                  f"Rate: {current_stats['sequences_per_second']:.1f} seq/s")

    return np.vstack(all_embeddings), stats, original_indices, skipped_df


def main():
    parser = argparse.ArgumentParser(description="Compute protein embeddings using ESM2")
    parser.add_argument("--limit", type=int, default=None,
                        help="Limit number of sequences to process (for testing)")
    parser.add_argument("--batch-size", type=int, default=32,
                        help="Batch size for embedding computation")
    parser.add_argument("--output", type=str, default=None,
                        help="Output file path (default: assets/proteins/datasets/all_species_embeddings.feather)")
    parser.add_argument("--checkpoint-dir", type=str, default=None,
                        help="Directory to save batch checkpoints")
    parser.add_argument("--max-seq-length", type=int, default=None,
                        help="Skip sequences longer than this (to avoid OOM)")
    parser.add_argument("--model", type=str, default=MODEL_CHECKPOINT,
                        help=f"HuggingFace model checkpoint (default: {MODEL_CHECKPOINT})")
    args = parser.parse_args()

    fasta_path = FASTA_FILE
    output_path = Path(args.output) if args.output else DEFAULT_OUTPUT
    checkpoint_dir = Path(args.checkpoint_dir) if args.checkpoint_dir else DEFAULT_CHECKPOINT_DIR
    model_checkpoint = args.model

    print(f"Model: {model_checkpoint}")

    # Load model and tokenizer
    print("Loading tokenizer and model...")
    tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)
    model = EsmModel.from_pretrained(model_checkpoint)

    # Load sequences
    sequence_df = load_sequences(fasta_path, limit=args.limit)

    # Compute embeddings
    embeddings, stats, original_indices, skipped_df = compute_embeddings_with_timing(
        sequence_df.copy(),
        tokenizer,
        model,
        batch_size=args.batch_size,
        checkpoint_dir=checkpoint_dir,
        max_seq_length=args.max_seq_length,
    )

    stats.print_summary()

    # Estimate full run time if this was a limited run
    if args.limit is not None:
        total_sequences = sum(1 for _ in SeqIO.parse(open(fasta_path), "fasta"))
        summary = stats.summary()
        estimated_total_time = total_sequences / summary["sequences_per_second"]
        print(f"\nFULL RUN ESTIMATE:")
        print(f"  Total sequences in dataset: {total_sequences}")
        print(f"  Estimated total time: {estimated_total_time / 60:.1f} minutes ({estimated_total_time / 3600:.2f} hours)")

    # Save embeddings (in original order)
    print(f"\nSaving embeddings to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    embedding_cols = [f"ME:{i+1}" for i in range(embeddings.shape[1])]
    embeddings_df = pd.DataFrame(embeddings, columns=embedding_cols)

    sorted_sequence_df = sequence_df.sort_values("Length").reset_index(drop=True)
    if args.max_seq_length is not None:
        sorted_sequence_df = sorted_sequence_df[sorted_sequence_df["Length"] <= args.max_seq_length].reset_index(drop=True)

    result_df = pd.concat([sorted_sequence_df, embeddings_df], axis=1)

    # Restore original order
    restore_order = np.argsort(original_indices)
    result_df = result_df.iloc[restore_order].reset_index(drop=True)

    result_df.to_feather(output_path)

    print(f"Saved {len(result_df)} embeddings with {embeddings.shape[1]} dimensions")

    # Report skipped sequences
    if skipped_df is not None and len(skipped_df) > 0:
        skipped_path = output_path.with_suffix(".skipped.csv")
        skipped_df.to_csv(skipped_path, index=False)
        print(f"Saved {len(skipped_df)} skipped sequences to {skipped_path}")

    print("Done!")


if __name__ == "__main__":
    main()
