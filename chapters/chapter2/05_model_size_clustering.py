#!/usr/bin/env python3
"""Multi-model clustering comparison with uniform protein set.

Loads ESM2 150M, 650M, and 3B embeddings, filters to common proteins,
runs KMeans + UMAP with multiple seeds, and computes silhouette scores
in both the original embedding space and UMAP 2D space.

Usage:
    python 05_model_size_clustering.py --output-dir chapters/chapter2/lab/011_uniform_clustering/results
    python 05_model_size_clustering.py --seeds 42,123,456,789
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from umap import UMAP

PROJECT_ROOT = Path(__file__).parent.parent.parent
DATASETS = PROJECT_ROOT / "assets/proteins/datasets"
DEFAULT_DATA = DATASETS / "cafa3_merged/cafa3_annotations.feather"

MODELS = {
    "150M": {
        "embeddings": DATASETS / "esm2_150m_embeddings.feather",
        "dim": 640,
    },
    "650M": {
        "embeddings": DATASETS / "all_species_embeddings.feather",
        "dim": 1280,
    },
    "3B": {
        "embeddings": DATASETS / "esm2_3b_embeddings.feather",
        "dim": 2560,
    },
}


def load_common_proteins(data_path: Path) -> tuple[set, pd.DataFrame]:
    """Find proteins common to all 3 embedding files and load annotations."""
    print("Finding common protein set across all models...")
    entry_sets = {}
    for name, cfg in MODELS.items():
        ids = pd.read_feather(cfg["embeddings"], columns=["EntryID"])["EntryID"]
        entry_sets[name] = set(ids)
        print(f"  {name}: {len(ids):,} proteins")

    common = entry_sets["150M"] & entry_sets["650M"] & entry_sets["3B"]
    print(f"  Common: {len(common):,} proteins")

    # Load annotations for sampling
    ann = pd.read_feather(data_path)
    protein_cols = ["EntryID", "Length", "taxonomyID", "scientific_name"]
    ann = ann[protein_cols].drop_duplicates(subset="EntryID")
    ann = ann[ann["EntryID"].isin(common)]
    print(f"  With annotations: {len(ann):,}")
    return common, ann


def get_taxa_and_sample(ann: pd.DataFrame, sample_size: int, seed: int,
                        min_taxa_count: int = 500) -> pd.DataFrame:
    """Balanced sample from abundant taxa."""
    taxa_counts = ann.groupby("taxonomyID").size()
    abundant = taxa_counts[taxa_counts >= min_taxa_count].index

    taxid_to_name = ann.groupby("taxonomyID")["scientific_name"].first().to_dict()

    rng = np.random.RandomState(seed)
    parts = []
    for taxid in abundant:
        sub = ann[ann["taxonomyID"] == taxid]
        n = min(sample_size, len(sub))
        parts.append(sub.sample(n=n, random_state=rng.randint(0, 2**31)))

    sample = pd.concat(parts, ignore_index=True)

    # Short species names
    sample["taxon_name"] = sample["taxonomyID"].map(
        lambda t: (lambda n: f"{n.split()[0][0]}. {n.split()[1]}" if n and len(n.split()) >= 2 else f"Tax_{t}")(taxid_to_name.get(t, ""))
    )
    print(f"  Sampled {len(sample):,} proteins from {len(abundant)} taxa (seed={seed})")
    return sample


def load_embeddings_for_sample(model_name: str, entry_ids: set) -> pd.DataFrame:
    """Load embeddings for a specific set of proteins."""
    cfg = MODELS[model_name]
    df = pd.read_feather(cfg["embeddings"])
    df = df[df["EntryID"].isin(entry_ids)]
    ecols = [c for c in df.columns if c.startswith("ME:")]
    return df[["EntryID"] + ecols]


def run_analysis(X: np.ndarray, n_clusters: int, seed: int) -> dict:
    """KMeans + UMAP, return metrics and coordinates."""
    t0 = time.time()

    # KMeans in original space
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    labels = kmeans.fit_predict(X)
    sil_orig = silhouette_score(X, labels, metric="cosine")
    t_kmeans = time.time() - t0

    # UMAP
    t1 = time.time()
    umap = UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                metric="cosine", random_state=seed, verbose=False)
    coords = umap.fit_transform(X)
    t_umap = time.time() - t1

    # Silhouette in UMAP space (same KMeans labels)
    sil_umap = silhouette_score(coords, labels, metric="euclidean")

    return {
        "kmeans_sil_original": float(sil_orig),
        "kmeans_sil_umap": float(sil_umap),
        "n_clusters": n_clusters,
        "n_proteins": int(X.shape[0]),
        "embedding_dim": int(X.shape[1]),
        "seed": seed,
        "time_kmeans_s": round(t_kmeans, 1),
        "time_umap_s": round(t_umap, 1),
        "umap_coords": coords,
        "kmeans_labels": labels,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--seeds", type=str, default="42,123,456,789")
    parser.add_argument("--sample-size", type=int, default=500)
    parser.add_argument("--n-clusters", type=int, default=20)
    parser.add_argument("--data", type=str, default=None)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seeds = [int(s) for s in args.seeds.split(",")]
    data_path = Path(args.data) if args.data else DEFAULT_DATA

    common, ann = load_common_proteins(data_path)

    all_results = []

    for seed in seeds:
        print(f"\n{'='*60}")
        print(f"SEED {seed}")
        print(f"{'='*60}")

        sample = get_taxa_and_sample(ann, args.sample_size, seed)
        sample_ids = set(sample["EntryID"])

        for model_name in ["150M", "650M", "3B"]:
            print(f"\n  --- {model_name} (seed={seed}) ---")
            emb_df = load_embeddings_for_sample(model_name, sample_ids)

            # Align to sample order
            emb_df = emb_df.set_index("EntryID").loc[sample["EntryID"].values].reset_index()
            ecols = [c for c in emb_df.columns if c.startswith("ME:")]
            X = emb_df[ecols].values
            print(f"  Embedding matrix: {X.shape}")

            result = run_analysis(X, args.n_clusters, seed)
            result["model"] = model_name

            print(f"  KMeans sil (original {X.shape[1]}D): {result['kmeans_sil_original']:.4f}")
            print(f"  KMeans sil (UMAP 2D):              {result['kmeans_sil_umap']:.4f}")

            # Save UMAP coordinates
            seed_model_dir = output_dir / f"seed_{seed}" / model_name
            seed_model_dir.mkdir(parents=True, exist_ok=True)
            coords_df = sample[["EntryID", "taxonomyID", "taxon_name"]].copy()
            coords_df["umap_x"] = result["umap_coords"][:, 0]
            coords_df["umap_y"] = result["umap_coords"][:, 1]
            coords_df["kmeans_cluster"] = result["kmeans_labels"]
            coords_df.to_csv(seed_model_dir / "umap_coordinates.csv", index=False)

            # Strip non-serializable fields before collecting
            row = {k: v for k, v in result.items()
                   if k not in ("umap_coords", "kmeans_labels")}
            all_results.append(row)

    # Summary table
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Model':<8} {'Seed':<6} {'Dim':<6} {'KMeans Sil (orig)':<20} {'KMeans Sil (UMAP)':<20}")
    print("-" * 60)
    for r in all_results:
        print(f"{r['model']:<8} {r['seed']:<6} {r['embedding_dim']:<6} "
              f"{r['kmeans_sil_original']:<20.4f} {r['kmeans_sil_umap']:<20.4f}")

    # Averages per model
    print(f"\n{'Model':<8} {'Avg Sil (orig)':<20} {'Avg Sil (UMAP)':<20} {'Std (orig)':<15} {'Std (UMAP)':<15}")
    print("-" * 78)
    for model_name in ["150M", "650M", "3B"]:
        rows = [r for r in all_results if r["model"] == model_name]
        orig = [r["kmeans_sil_original"] for r in rows]
        umap_s = [r["kmeans_sil_umap"] for r in rows]
        print(f"{model_name:<8} {np.mean(orig):<20.4f} {np.mean(umap_s):<20.4f} "
              f"{np.std(orig):<15.4f} {np.std(umap_s):<15.4f}")

    # Save results JSON
    with open(output_dir / "results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {output_dir / 'results.json'}")


if __name__ == "__main__":
    main()
