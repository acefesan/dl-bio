#!/usr/bin/env python3
"""GO function prediction from ESM2 embeddings using a Flax/JAX MLP.

Trains a simple sequential model (Dense → GELU → Dense → GELU → Dense) to predict
GO term annotations from protein embeddings. Compares training performance across
different ESM2 embedding sizes (150M/640D, 650M/1280D, 3B/2560D).

Uses the same train/valid/test split from the preprocessed 150M dataset, swapping
in embeddings from each model.

Usage:
    python 05_go_prediction.py --output-dir results/go_prediction
    python 05_go_prediction.py --models 150m,650m    # subset of models
    python 05_go_prediction.py --num-steps 500        # more training steps
"""

import argparse
import json
import time
from pathlib import Path

import flax.linen as nn
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
import optax
import pandas as pd
from flax.training.train_state import TrainState
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    precision_score,
    recall_score,
    roc_auc_score,
)

PROJECT_ROOT = Path(__file__).parent.parent.parent

EMBEDDING_PATHS = {
    "150m": PROJECT_ROOT / "assets/proteins/datasets/esm2_150m_embeddings.feather",
    "650m": PROJECT_ROOT / "assets/proteins/datasets/all_species_embeddings.feather",
    "3b": PROJECT_ROOT / "assets/proteins/datasets/esm2_3b_embeddings.feather",
}

DATASET_DIR = PROJECT_ROOT / "assets/proteins/datasets"

EMBEDDING_DIMS = {"150m": 640, "650m": 1280, "3b": 2560}


# =============================================================================
# Model (same architecture as dlfb.proteins.model.Model)
# =============================================================================

class GOPredictor(nn.Module):
    """Simple MLP for GO function prediction."""
    num_targets: int
    hidden_dim: int = 256

    @nn.compact
    def __call__(self, x):
        x = nn.Sequential([
            nn.Dense(self.hidden_dim * 2),
            jax.nn.gelu,
            nn.Dense(self.hidden_dim),
            jax.nn.gelu,
            nn.Dense(self.num_targets),
        ])(x)
        return x


# =============================================================================
# Data loading
# =============================================================================

def load_split(split: str, model_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Load GO labels from preprocessed 150M dataset, swap in embeddings from model_name."""
    # Load the preprocessed dataset (has GO: and ME: columns for 150M)
    df = pd.read_feather(
        DATASET_DIR / f"protein_dataset_{split}_facebook_esm2_t30_150M_UR50D.feather"
    )

    # Extract GO targets
    go_cols = sorted([c for c in df.columns if c.startswith("GO:")])
    targets = df[go_cols].values.astype(np.float32)

    # Load embeddings for the requested model
    emb_df = pd.read_feather(EMBEDDING_PATHS[model_name])
    me_cols = sorted([c for c in emb_df.columns if c.startswith("ME:")])
    emb_df = emb_df[["EntryID"] + me_cols]

    # Join on EntryID
    merged = df[["EntryID"]].merge(emb_df, on="EntryID", how="inner")
    embeddings = merged[me_cols].values.astype(np.float32)

    return embeddings, targets


def make_batches(embeddings, targets, batch_size, shuffle=False, rng=None):
    """Yield batches from numpy arrays."""
    n = len(embeddings)
    indices = np.arange(n)
    if shuffle and rng is not None:
        rng.shuffle(indices)
    for start in range(0, n, batch_size):
        idx = indices[start:start + batch_size]
        if len(idx) < batch_size:
            continue
        yield {"embedding": embeddings[idx], "target": targets[idx]}


# =============================================================================
# Training
# =============================================================================

@jax.jit
def train_step(state, batch):
    """Single training step."""
    def loss_fn(params):
        logits = state.apply_fn({"params": params}, batch["embedding"])
        loss = optax.sigmoid_binary_cross_entropy(logits, batch["target"]).mean()
        return loss

    loss, grads = jax.value_and_grad(loss_fn)(state.params)
    state = state.apply_gradients(grads=grads)
    return state, loss


def compute_per_target_metrics(targets: np.ndarray, probs: np.ndarray, thresh: float = 0.5) -> dict:
    """Compute per-target metrics (same as dlfb chapter 2), then average."""
    metric_names = ["accuracy", "recall", "precision", "auprc", "auroc"]
    per_target = {m: [] for m in metric_names}

    for i in range(targets.shape[1]):
        t, p = targets[:, i], probs[:, i]
        if t.sum() == 0:
            continue
        preds = (p >= thresh).astype(int)
        per_target["accuracy"].append(accuracy_score(t, preds))
        per_target["recall"].append(recall_score(t, preds, zero_division=0.0))
        per_target["precision"].append(precision_score(t, preds, zero_division=0.0))
        per_target["auprc"].append(average_precision_score(t, p))
        try:
            per_target["auroc"].append(roc_auc_score(t, p))
        except ValueError:
            pass

    return {m: float(np.mean(v)) if v else 0.0 for m, v in per_target.items()}


def evaluate(state, embeddings, targets, batch_size):
    """Evaluate on a full split, return loss and per-target metrics dict."""
    losses = []
    all_probs = []
    all_targets = []

    for batch in make_batches(embeddings, targets, batch_size):
        logits = state.apply_fn({"params": state.params}, batch["embedding"])
        loss = optax.sigmoid_binary_cross_entropy(logits, batch["target"]).mean()
        losses.append(float(loss))
        all_probs.append(jax.nn.sigmoid(logits))
        all_targets.append(batch["target"])

    mean_loss = np.mean(losses)
    probs = np.vstack(all_probs)
    tgts = np.vstack(all_targets)
    metrics = compute_per_target_metrics(tgts, probs)
    metrics["loss"] = float(mean_loss)

    return metrics


def train_model(
    model_name: str,
    num_steps: int = 300,
    batch_size: int = 64,
    lr: float = 1e-3,
    eval_every: int = 30,
    seed: int = 42,
    hidden_dim: int = 256,
) -> dict:
    """Train GO predictor with given embedding model, return metrics."""
    print(f"\n{'=' * 60}")
    print(f"Training with {model_name.upper()} embeddings ({EMBEDDING_DIMS[model_name]}D)")
    print(f"{'=' * 60}")

    # Load data
    train_emb, train_tgt = load_split("train", model_name)
    valid_emb, valid_tgt = load_split("valid", model_name)
    test_emb, test_tgt = load_split("test", model_name)

    n_targets = train_tgt.shape[1]
    emb_dim = train_emb.shape[1]
    print(f"  Train: {len(train_emb)} proteins, {emb_dim}D embeddings, {n_targets} GO targets")
    print(f"  Valid: {len(valid_emb)} proteins")
    print(f"  Test:  {len(test_emb)} proteins")

    # Init model
    model = GOPredictor(num_targets=n_targets, hidden_dim=hidden_dim)
    rng = jax.random.PRNGKey(seed)
    dummy = jnp.ones((1, emb_dim))
    variables = model.init(rng, dummy)
    tx = optax.adam(lr)
    state = TrainState.create(apply_fn=model.apply, params=variables["params"], tx=tx)

    n_params = sum(x.size for x in jax.tree.leaves(state.params))
    print(f"  Model params: {n_params:,}")

    # Training loop
    train_metrics = []
    valid_metrics = []
    np_rng = np.random.default_rng(seed)
    t0 = time.time()

    for step in range(num_steps):
        # Get one epoch of shuffled batches, cycle through
        for batch in make_batches(train_emb, train_tgt, batch_size, shuffle=True, rng=np_rng):
            state, loss = train_step(state, batch)
            train_metrics.append({"step": step, "loss": float(loss)})
            break  # one batch per step

        if step % eval_every == 0 or step == num_steps - 1:
            val_metrics = evaluate(state, valid_emb, valid_tgt, batch_size)
            valid_metrics.append({"step": step, **val_metrics})
            print(f"  Step {step:4d} | train_loss: {float(loss):.4f} | val_loss: {val_metrics['loss']:.4f} | val_auPRC: {val_metrics['auprc']:.4f}")

    elapsed = time.time() - t0

    # Final test evaluation
    test_metrics = evaluate(state, test_emb, test_tgt, batch_size)
    print(f"\n  Test loss: {test_metrics['loss']:.4f} | Test auPRC: {test_metrics['auprc']:.4f}")
    print(f"  Training time: {elapsed:.1f}s")

    return {
        "model_name": model_name,
        "embedding_dim": emb_dim,
        "n_targets": n_targets,
        "n_params": int(n_params),
        "num_steps": num_steps,
        "hidden_dim": hidden_dim,
        "lr": lr,
        "batch_size": batch_size,
        "seed": seed,
        "train_metrics": train_metrics,
        "valid_metrics": valid_metrics,
        "test_metrics": test_metrics,
        "train_time_s": round(elapsed, 1),
        "best_val_loss": min(m["loss"] for m in valid_metrics),
        "best_val_auprc": max(m["auprc"] for m in valid_metrics),
    }


# =============================================================================
# Plotting
# =============================================================================

def plot_comparison(all_results: dict, output_dir: Path):
    """Training curves for all models on one figure."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = {"150m": "#2196F3", "650m": "#4CAF50", "3b": "#F44336"}

    for model_name, results in all_results.items():
        c = colors.get(model_name, "gray")
        label = f"{model_name.upper()} ({results['embedding_dim']}D)"

        # Train loss
        train_df = pd.DataFrame(results["train_metrics"])
        axes[0].plot(train_df["step"], train_df["loss"], alpha=0.4, c=c, linewidth=0.5)
        # Smoothed
        smooth = train_df["loss"].rolling(20, min_periods=1).mean()
        axes[0].plot(train_df["step"], smooth, c=c, label=label, linewidth=2)

        # Valid loss
        valid_df = pd.DataFrame(results["valid_metrics"])
        axes[1].plot(valid_df["step"], valid_df["loss"], c=c, label=label, linewidth=2, marker="o", markersize=3)

        # Valid auPRC
        axes[2].plot(valid_df["step"], valid_df["auprc"], c=c, label=label, linewidth=2, marker="o", markersize=3)

    axes[0].set_title("Train Loss", fontsize=13)
    axes[0].set_xlabel("Step")
    axes[0].set_ylabel("Loss")
    axes[0].legend(fontsize=9)

    axes[1].set_title("Validation Loss", fontsize=13)
    axes[1].set_xlabel("Step")
    axes[1].set_ylabel("Loss")
    axes[1].legend(fontsize=9)

    axes[2].set_title("Validation auPRC", fontsize=13)
    axes[2].set_xlabel("Step")
    axes[2].set_ylabel("Mean auPRC")
    axes[2].legend(fontsize=9)

    plt.suptitle("GO Function Prediction: ESM2 Embedding Comparison", fontsize=14, y=1.02)
    plt.tight_layout()
    out_path = output_dir / "go_prediction_comparison.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {out_path}")
    return out_path


def plot_final_bar(all_results: dict, output_dir: Path):
    """Bar chart of final test metrics."""
    models = list(all_results.keys())
    metric_keys = ["loss", "accuracy", "recall", "precision", "auprc", "auroc"]
    colors = ["#2196F3", "#4CAF50", "#F44336"]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    labels = [f"{m.upper()}\n({all_results[m]['embedding_dim']}D)" for m in models]

    for ax_idx, metric in enumerate(metric_keys):
        values = [all_results[m]["test_metrics"][metric] for m in models]
        axes[ax_idx].bar(labels, values, color=colors[:len(models)])
        axes[ax_idx].set_title(f"Test {metric}", fontsize=13)
        axes[ax_idx].set_ylabel(metric)
        for i, v in enumerate(values):
            axes[ax_idx].text(i, v + 0.002, f"{v:.4f}", ha="center", fontsize=10)

    plt.suptitle("GO Prediction: Final Test Performance", fontsize=14)
    plt.tight_layout()
    out_path = output_dir / "go_prediction_test_bars.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")
    return out_path


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="GO prediction from ESM2 embeddings")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--models", type=str, default="150m,650m,3b",
                        help="Comma-separated model names (default: 150m,650m,3b)")
    parser.add_argument("--num-steps", type=int, default=300, help="Training steps (default: 300)")
    parser.add_argument("--batch-size", type=int, default=64, help="Batch size (default: 64)")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate (default: 1e-3)")
    parser.add_argument("--hidden-dim", type=int, default=256, help="MLP hidden dim (default: 256)")
    parser.add_argument("--eval-every", type=int, default=30, help="Eval frequency (default: 30)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_names = [m.strip() for m in args.models.split(",")]

    all_results = {}
    for model_name in model_names:
        results = train_model(
            model_name=model_name,
            num_steps=args.num_steps,
            batch_size=args.batch_size,
            lr=args.lr,
            eval_every=args.eval_every,
            seed=args.seed,
            hidden_dim=args.hidden_dim,
        )
        all_results[model_name] = results

        # Save per-model results
        with open(output_dir / f"results_{model_name}.json", "w") as f:
            json.dump(results, f, indent=2)

    # Comparison plots
    plot_comparison(all_results, output_dir)
    plot_final_bar(all_results, output_dir)

    # Summary table
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    print(f"{'Model':<8} {'Dim':>6} {'Params':>10} {'Loss':>8} {'Acc':>8} {'Rec':>8} {'Prec':>8} {'auPRC':>8} {'auROC':>8} {'Time':>6}")
    print("-" * 80)
    for m, r in all_results.items():
        t = r["test_metrics"]
        print(f"{m.upper():<8} {r['embedding_dim']:>6} {r['n_params']:>10,} {t['loss']:>8.4f} {t['accuracy']:>8.4f} {t['recall']:>8.4f} {t['precision']:>8.4f} {t['auprc']:>8.4f} {t['auroc']:>8.4f} {r['train_time_s']:>5.0f}s")

    # Save combined summary
    summary = {m: {k: v for k, v in r.items() if k not in ("train_metrics", "valid_metrics")}
               for m, r in all_results.items()}
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
