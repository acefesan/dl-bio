#!/usr/bin/env python3
"""Grid search over classifier architecture for GO prediction.

Tests whether 3B's poor performance (entry 007) is due to an undersized classifier.
Sweeps hidden_dim and num_layers across all 3 embedding sizes.

Usage:
    python 06_go_grid_search.py --output-dir results/go_grid_search
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
# Model — variable depth
# =============================================================================

class GOPredictor(nn.Module):
    """MLP with configurable depth for GO function prediction."""
    num_targets: int
    hidden_dim: int = 256
    num_layers: int = 2  # number of hidden layers

    @nn.compact
    def __call__(self, x):
        for _ in range(self.num_layers):
            x = nn.Dense(self.hidden_dim)(x)
            x = jax.nn.gelu(x)
        x = nn.Dense(self.num_targets)(x)
        return x


# =============================================================================
# Data loading (reuse from 05_go_prediction.py)
# =============================================================================

_data_cache = {}

def load_split(split: str, model_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Load GO labels from preprocessed 150M dataset, swap in embeddings."""
    key = (split, model_name)
    if key in _data_cache:
        return _data_cache[key]

    df = pd.read_feather(
        DATASET_DIR / f"protein_dataset_{split}_facebook_esm2_t30_150M_UR50D.feather"
    )
    go_cols = sorted([c for c in df.columns if c.startswith("GO:")])
    targets = df[go_cols].values.astype(np.float32)

    emb_df = pd.read_feather(EMBEDDING_PATHS[model_name])
    me_cols = sorted([c for c in emb_df.columns if c.startswith("ME:")])
    emb_df = emb_df[["EntryID"] + me_cols]

    merged = df[["EntryID"]].merge(emb_df, on="EntryID", how="inner")
    embeddings = merged[me_cols].values.astype(np.float32)

    _data_cache[key] = (embeddings, targets)
    return embeddings, targets


def make_batches(embeddings, targets, batch_size, shuffle=False, rng=None):
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
    def loss_fn(params):
        logits = state.apply_fn({"params": params}, batch["embedding"])
        loss = optax.sigmoid_binary_cross_entropy(logits, batch["target"]).mean()
        return loss
    loss, grads = jax.value_and_grad(loss_fn)(state.params)
    state = state.apply_gradients(grads=grads)
    return state, loss


def compute_per_target_metrics(targets: np.ndarray, probs: np.ndarray, thresh: float = 0.5) -> dict:
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
    losses, all_probs, all_targets = [], [], []
    for batch in make_batches(embeddings, targets, batch_size):
        logits = state.apply_fn({"params": state.params}, batch["embedding"])
        loss = optax.sigmoid_binary_cross_entropy(logits, batch["target"]).mean()
        losses.append(float(loss))
        all_probs.append(jax.nn.sigmoid(logits))
        all_targets.append(batch["target"])
    probs = np.vstack(all_probs)
    tgts = np.vstack(all_targets)
    metrics = compute_per_target_metrics(tgts, probs)
    metrics["loss"] = float(np.mean(losses))
    return metrics


def train_single(
    model_name: str,
    hidden_dim: int,
    num_layers: int,
    num_steps: int = 300,
    batch_size: int = 64,
    lr: float = 1e-3,
    seed: int = 42,
) -> dict:
    """Train one configuration, return test metrics."""
    train_emb, train_tgt = load_split("train", model_name)
    valid_emb, valid_tgt = load_split("valid", model_name)
    test_emb, test_tgt = load_split("test", model_name)

    n_targets = train_tgt.shape[1]
    emb_dim = train_emb.shape[1]

    model = GOPredictor(num_targets=n_targets, hidden_dim=hidden_dim, num_layers=num_layers)
    rng = jax.random.PRNGKey(seed)
    variables = model.init(rng, jnp.ones((1, emb_dim)))
    tx = optax.adam(lr)
    state = TrainState.create(apply_fn=model.apply, params=variables["params"], tx=tx)

    n_params = sum(x.size for x in jax.tree.leaves(state.params))
    np_rng = np.random.default_rng(seed)
    t0 = time.time()

    best_val_auprc = 0.0
    for step in range(num_steps):
        for batch in make_batches(train_emb, train_tgt, batch_size, shuffle=True, rng=np_rng):
            state, loss = train_step(state, batch)
            break

        if step % 50 == 0 or step == num_steps - 1:
            val = evaluate(state, valid_emb, valid_tgt, batch_size)
            best_val_auprc = max(best_val_auprc, val["auprc"])

    elapsed = time.time() - t0
    test = evaluate(state, test_emb, test_tgt, batch_size)

    return {
        "model_name": model_name,
        "embedding_dim": emb_dim,
        "hidden_dim": hidden_dim,
        "num_layers": num_layers,
        "n_params": int(n_params),
        "test_metrics": test,
        "best_val_auprc": best_val_auprc,
        "train_time_s": round(elapsed, 1),
    }


# =============================================================================
# Plotting
# =============================================================================

def plot_grid(results: list[dict], output_dir: Path):
    """Heatmaps of test auPRC: hidden_dim × num_layers, one per embedding."""
    model_names = sorted(set(r["model_name"] for r in results))
    hidden_dims = sorted(set(r["hidden_dim"] for r in results))
    num_layers_list = sorted(set(r["num_layers"] for r in results))
    colors = {"150m": "Blues", "650m": "Greens", "3b": "Reds"}

    fig, axes = plt.subplots(1, len(model_names), figsize=(6 * len(model_names), 5))
    if len(model_names) == 1:
        axes = [axes]

    for ax, mn in zip(axes, model_names):
        grid = np.zeros((len(hidden_dims), len(num_layers_list)))
        for r in results:
            if r["model_name"] != mn:
                continue
            hi = hidden_dims.index(r["hidden_dim"])
            li = num_layers_list.index(r["num_layers"])
            grid[hi, li] = r["test_metrics"]["auprc"]

        im = ax.imshow(grid, cmap=colors.get(mn, "viridis"), aspect="auto",
                       vmin=0, vmax=max(r["test_metrics"]["auprc"] for r in results) * 1.05)
        ax.set_xticks(range(len(num_layers_list)))
        ax.set_xticklabels([str(n) for n in num_layers_list])
        ax.set_yticks(range(len(hidden_dims)))
        ax.set_yticklabels([str(h) for h in hidden_dims])
        ax.set_xlabel("Num hidden layers")
        ax.set_ylabel("Hidden dim")
        ax.set_title(f"{mn.upper()} ({EMBEDDING_DIMS[mn]}D)", fontsize=13)

        for i in range(len(hidden_dims)):
            for j in range(len(num_layers_list)):
                ax.text(j, i, f"{grid[i, j]:.3f}", ha="center", va="center",
                        fontsize=10, fontweight="bold",
                        color="white" if grid[i, j] > grid.max() * 0.6 else "black")

        plt.colorbar(im, ax=ax, label="Test auPRC")

    plt.suptitle("GO Prediction: Classifier Architecture Grid Search", fontsize=14, y=1.02)
    plt.tight_layout()
    out_path = output_dir / "go_grid_search_heatmap.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")
    return out_path


def plot_scaling(results: list[dict], output_dir: Path):
    """auPRC vs classifier params, colored by embedding model."""
    colors = {"150m": "#2196F3", "650m": "#4CAF50", "3b": "#F44336"}
    markers = {1: "o", 2: "s", 3: "^"}

    fig, ax = plt.subplots(figsize=(10, 6))
    for r in results:
        ax.scatter(r["n_params"], r["test_metrics"]["auprc"],
                   c=colors[r["model_name"]], marker=markers[r["num_layers"]],
                   s=80, alpha=0.8, edgecolors="black", linewidth=0.5)

    # Legend for embedding models
    for mn, c in colors.items():
        ax.scatter([], [], c=c, s=80, label=f"{mn.upper()} ({EMBEDDING_DIMS[mn]}D)")
    # Legend for layers
    for nl, m in markers.items():
        ax.scatter([], [], c="gray", marker=m, s=80, label=f"{nl} hidden layer{'s' if nl > 1 else ''}")

    ax.set_xlabel("Classifier params", fontsize=12)
    ax.set_ylabel("Test auPRC", fontsize=12)
    ax.set_xscale("log")
    ax.set_title("GO Prediction: auPRC vs Classifier Size", fontsize=14)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = output_dir / "go_grid_search_scaling.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")
    return out_path


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Grid search GO prediction classifier")
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--models", type=str, default="150m,650m,3b")
    parser.add_argument("--hidden-dims", type=str, default="128,256,512,1024")
    parser.add_argument("--num-layers", type=str, default="1,2,3")
    parser.add_argument("--num-steps", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_names = [m.strip() for m in args.models.split(",")]
    hidden_dims = [int(h) for h in args.hidden_dims.split(",")]
    num_layers_list = [int(n) for n in args.num_layers.split(",")]

    total = len(model_names) * len(hidden_dims) * len(num_layers_list)
    print(f"Grid search: {len(model_names)} models × {len(hidden_dims)} hidden_dims × {len(num_layers_list)} depths = {total} runs")

    all_results = []
    for i, model_name in enumerate(model_names):
        for hidden_dim in hidden_dims:
            for num_layers in num_layers_list:
                run_idx = len(all_results) + 1
                print(f"\n[{run_idx}/{total}] {model_name.upper()} | hidden={hidden_dim} | layers={num_layers}")
                result = train_single(
                    model_name=model_name,
                    hidden_dim=hidden_dim,
                    num_layers=num_layers,
                    num_steps=args.num_steps,
                    batch_size=args.batch_size,
                    lr=args.lr,
                    seed=args.seed,
                )
                auprc = result["test_metrics"]["auprc"]
                auroc = result["test_metrics"]["auroc"]
                print(f"  → params={result['n_params']:,} | auPRC={auprc:.4f} | auROC={auroc:.4f} | {result['train_time_s']}s")
                all_results.append(result)

    # Save raw results
    with open(output_dir / "grid_results.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Plots
    plot_grid(all_results, output_dir)
    plot_scaling(all_results, output_dir)

    # Summary table
    print(f"\n{'=' * 90}")
    print("GRID SEARCH SUMMARY")
    print(f"{'=' * 90}")
    print(f"{'Model':<6} {'Dim':>5} {'Layers':>6} {'Hidden':>7} {'Params':>10} {'Loss':>7} {'auPRC':>7} {'auROC':>7}")
    print("-" * 65)
    for r in sorted(all_results, key=lambda x: (-x["test_metrics"]["auprc"])):
        t = r["test_metrics"]
        print(f"{r['model_name'].upper():<6} {r['embedding_dim']:>5} {r['num_layers']:>6} {r['hidden_dim']:>7} {r['n_params']:>10,} {t['loss']:>7.4f} {t['auprc']:>7.4f} {t['auroc']:>7.4f}")

    # Best per model
    print(f"\n{'=' * 50}")
    print("BEST PER EMBEDDING MODEL")
    print(f"{'=' * 50}")
    for mn in model_names:
        best = max((r for r in all_results if r["model_name"] == mn),
                   key=lambda x: x["test_metrics"]["auprc"])
        t = best["test_metrics"]
        print(f"  {mn.upper()}: layers={best['num_layers']} hidden={best['hidden_dim']} "
              f"params={best['n_params']:,} → auPRC={t['auprc']:.4f} auROC={t['auroc']:.4f}")

    with open(output_dir / "summary.json", "w") as f:
        json.dump({
            "grid": {"models": model_names, "hidden_dims": hidden_dims, "num_layers": num_layers_list},
            "results": [{k: v for k, v in r.items()} for r in all_results],
        }, f, indent=2)

    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
