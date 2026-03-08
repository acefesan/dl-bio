#!/usr/bin/env python3
"""Analyze training dynamics from ESM2 pretraining runs.

Loads checkpoint data and computes metrics beyond simple train/eval loss:
  1. Loss curves + generalization gap
  2. Weight norm evolution per layer
  3. Effective rank of weight matrices (spectral analysis)
  4. Gradient norm trends
  5. Per-layer spectral distribution snapshots

Usage:
    python chapters/chapter2/analyze_training_dynamics.py <run_dir> [<run_dir2> ...]
    python chapters/chapter2/analyze_training_dynamics.py runs/pretrain_esm2_* --compare
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch

CHAPTER_DIR = Path(__file__).parent


def load_training_log(run_dir: Path) -> dict:
    """Extract train/eval metrics from the last checkpoint's trainer_state.json."""
    # Find the last checkpoint
    checkpoints = sorted(
        [d for d in run_dir.iterdir() if d.name.startswith("checkpoint-")],
        key=lambda d: int(d.name.split("-")[1]),
    )
    if not checkpoints:
        print(f"No checkpoints found in {run_dir}")
        sys.exit(1)

    state_file = checkpoints[-1] / "trainer_state.json"
    with open(state_file) as f:
        state = json.load(f)

    train_logs = []
    eval_logs = []
    for entry in state["log_history"]:
        if "eval_loss" in entry:
            eval_logs.append(entry)
        elif "loss" in entry and "train_loss" not in entry:
            train_logs.append(entry)

    return {
        "train": train_logs,
        "eval": eval_logs,
        "checkpoints": checkpoints,
    }


def plot_loss_curves(logs: dict, run_name: str, ax: plt.Axes):
    """Plot train loss, eval loss, and generalization gap."""
    train_steps = [e["step"] for e in logs["train"]]
    train_loss = [e["loss"] for e in logs["train"]]
    eval_steps = [e["step"] for e in logs["eval"]]
    eval_loss = [e["eval_loss"] for e in logs["eval"]]

    ax.plot(train_steps, train_loss, alpha=0.5, label=f"{run_name} train", linewidth=0.8)
    ax.plot(eval_steps, eval_loss, label=f"{run_name} eval", linewidth=1.5)
    ax.set_xlabel("Step")
    ax.set_ylabel("Loss")
    ax.legend(fontsize=8)
    ax.set_title("Train vs Eval Loss")
    ax.grid(True, alpha=0.3)


def plot_generalization_gap(logs: dict, run_name: str, ax: plt.Axes):
    """Plot the gap between train and eval loss over time."""
    # Interpolate train loss at eval steps
    train_steps = np.array([e["step"] for e in logs["train"]])
    train_loss = np.array([e["loss"] for e in logs["train"]])
    eval_steps = np.array([e["step"] for e in logs["eval"]])
    eval_loss = np.array([e["eval_loss"] for e in logs["eval"]])

    train_at_eval = np.interp(eval_steps, train_steps, train_loss)
    gap = eval_loss - train_at_eval

    ax.plot(eval_steps, gap, label=run_name, linewidth=1.5)
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Step")
    ax.set_ylabel("Eval Loss - Train Loss")
    ax.set_title("Generalization Gap")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_grad_norm(logs: dict, run_name: str, ax: plt.Axes):
    """Plot gradient norm evolution."""
    steps = [e["step"] for e in logs["train"] if "grad_norm" in e]
    grad_norms = [e["grad_norm"] for e in logs["train"] if "grad_norm" in e]

    ax.plot(steps, grad_norms, alpha=0.5, label=run_name, linewidth=0.8)
    ax.set_xlabel("Step")
    ax.set_ylabel("Gradient Norm")
    ax.set_title("Gradient Norm")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


def compute_weight_metrics(checkpoint_dir: Path) -> dict:
    """Compute weight norm and effective rank for each layer from a checkpoint."""
    model_path = checkpoint_dir / "model.safetensors"
    if not model_path.exists():
        return None

    from safetensors.torch import load_file
    state_dict = load_file(str(model_path))

    metrics = {}
    total_norm = 0.0
    for name, param in state_dict.items():
        if param.ndim < 2:
            continue
        # Reshape to 2D for SVD
        w = param.float().reshape(param.shape[0], -1)
        norm = torch.norm(w).item()
        total_norm += norm ** 2

        # Effective rank via Shannon entropy of normalized singular values
        try:
            s = torch.linalg.svdvals(w)
            s_norm = s / s.sum()
            s_norm = s_norm[s_norm > 1e-10]  # filter zeros
            entropy = -(s_norm * torch.log(s_norm)).sum().item()
            eff_rank = np.exp(entropy)
        except Exception:
            eff_rank = float("nan")

        # Stable rank: ||W||_F^2 / ||W||_2^2
        stable_rank = (norm ** 2) / (s[0].item() ** 2) if s[0].item() > 0 else float("nan")

        metrics[name] = {
            "frobenius_norm": norm,
            "spectral_norm": s[0].item(),
            "effective_rank": eff_rank,
            "stable_rank": stable_rank,
            "max_rank": min(w.shape),
        }

    metrics["_total_weight_norm"] = np.sqrt(total_norm)
    return metrics


def analyze_weight_dynamics(checkpoints: list[Path], sample_every: int = 10):
    """Analyze weight metrics across sampled checkpoints."""
    sampled = checkpoints[::sample_every]
    if checkpoints[-1] not in sampled:
        sampled.append(checkpoints[-1])

    results = []
    for ckpt in sampled:
        step = int(ckpt.name.split("-")[1])
        print(f"  Loading checkpoint {step}...")
        metrics = compute_weight_metrics(ckpt)
        if metrics:
            results.append((step, metrics))

    return results


def plot_weight_norms(weight_results: list, run_name: str, ax: plt.Axes):
    """Plot total weight norm over training."""
    steps = [r[0] for r in weight_results]
    norms = [r[1]["_total_weight_norm"] for r in weight_results]
    ax.plot(steps, norms, "o-", label=run_name, markersize=3)
    ax.set_xlabel("Step")
    ax.set_ylabel("Total Weight Norm (Frobenius)")
    ax.set_title("Weight Norm Growth")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_effective_rank(weight_results: list, run_name: str, ax: plt.Axes):
    """Plot mean effective rank across attention layers."""
    steps = []
    mean_ranks = []
    for step, metrics in weight_results:
        ranks = []
        for name, m in metrics.items():
            if name.startswith("_"):
                continue
            if "attention" in name and "weight" in name:
                ranks.append(m["effective_rank"] / m["max_rank"])  # normalized
        if ranks:
            steps.append(step)
            mean_ranks.append(np.mean(ranks))

    ax.plot(steps, mean_ranks, "o-", label=run_name, markersize=3)
    ax.set_xlabel("Step")
    ax.set_ylabel("Normalized Effective Rank")
    ax.set_title("Attention Layer Effective Rank")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_layer_rank_heatmap(weight_results: list, run_name: str, ax: plt.Axes):
    """Heatmap of effective rank per layer over training."""
    # Collect attention query weight effective ranks per layer
    layer_names = []
    for name in weight_results[0][1]:
        if name.startswith("_"):
            continue
        if "attention.self.query.weight" in name:
            layer_names.append(name)
    layer_names.sort()

    if not layer_names:
        ax.text(0.5, 0.5, "No attention layers found", ha="center", va="center")
        return

    steps = [r[0] for r in weight_results]
    matrix = np.zeros((len(layer_names), len(steps)))
    for j, (step, metrics) in enumerate(weight_results):
        for i, name in enumerate(layer_names):
            if name in metrics:
                matrix[i, j] = metrics[name]["stable_rank"]

    im = ax.imshow(matrix, aspect="auto", cmap="viridis",
                   extent=[steps[0], steps[-1], len(layer_names) - 0.5, -0.5])
    ax.set_xlabel("Step")
    ax.set_ylabel("Layer")
    ax.set_yticks(range(len(layer_names)))
    ax.set_yticklabels([f"L{i}" for i in range(len(layer_names))], fontsize=7)
    ax.set_title(f"Stable Rank per Layer ({run_name})")
    plt.colorbar(im, ax=ax, label="Stable Rank")


def main():
    parser = argparse.ArgumentParser(description="Analyze ESM2 training dynamics")
    parser.add_argument("run_dirs", nargs="+", type=str, help="Run directories to analyze")
    parser.add_argument("--sample-every", type=int, default=10,
                        help="Sample every N checkpoints for weight analysis (default: 10)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output path for figure (default: <first_run>/training_dynamics.png)")
    args = parser.parse_args()

    run_dirs = [Path(d) for d in args.run_dirs]
    for d in run_dirs:
        if not d.exists():
            print(f"ERROR: {d} not found")
            sys.exit(1)

    output_path = Path(args.output) if args.output else run_dirs[0] / "training_dynamics.png"
    n_runs = len(run_dirs)

    # --- Load logs ---
    all_logs = {}
    for run_dir in run_dirs:
        name = run_dir.name
        print(f"Loading logs from {name}...")
        all_logs[name] = load_training_log(run_dir)

    # --- Create figure ---
    n_cols = 3
    n_rows = 2 + n_runs  # loss/gap/grad row + weight rows per run
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows))
    fig.suptitle("Training Dynamics Analysis", fontsize=14, y=0.98)

    # Row 1: Loss curves, generalization gap, gradient norms
    for name, logs in all_logs.items():
        plot_loss_curves(logs, name, axes[0, 0])
        plot_generalization_gap(logs, name, axes[0, 1])
        plot_grad_norm(logs, name, axes[0, 2])

    # Row 2+: Weight analysis per run
    for i, (name, logs) in enumerate(all_logs.items()):
        print(f"\nAnalyzing weights for {name}...")
        weight_results = analyze_weight_dynamics(
            logs["checkpoints"], sample_every=args.sample_every
        )

        row = 1 + i
        plot_weight_norms(weight_results, name, axes[row, 0])
        plot_effective_rank(weight_results, name, axes[row, 1])
        plot_layer_rank_heatmap(weight_results, name, axes[row, 2])

    # Last row: comparison overlay if multiple runs
    if n_runs > 1:
        row = 1 + n_runs
        all_weight_results = {}
        for name, logs in all_logs.items():
            all_weight_results[name] = analyze_weight_dynamics(
                logs["checkpoints"], sample_every=args.sample_every
            )
        for name, wr in all_weight_results.items():
            plot_weight_norms(wr, name, axes[row, 0])
            plot_effective_rank(wr, name, axes[row, 1])
        axes[row, 0].set_title("Weight Norm Comparison")
        axes[row, 1].set_title("Effective Rank Comparison")
        axes[row, 2].axis("off")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {output_path}")

    # --- Print summary stats ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, logs in all_logs.items():
        train = logs["train"]
        evals = logs["eval"]
        print(f"\n{name}:")
        print(f"  Steps:           {train[-1]['step']}")
        print(f"  Epochs:          {train[-1]['epoch']:.1f}")
        print(f"  Final train loss: {train[-1]['loss']:.4f}")
        print(f"  Final eval loss:  {evals[-1]['eval_loss']:.4f}")
        print(f"  Gen. gap:         {evals[-1]['eval_loss'] - train[-1]['loss']:+.4f}")
        print(f"  Best eval loss:   {min(e['eval_loss'] for e in evals):.4f} "
              f"(step {min(evals, key=lambda e: e['eval_loss'])['step']})")


if __name__ == "__main__":
    main()
