#!/usr/bin/env python3
"""Investigation of PCA and Whitening on ESM2 embedding clustering.

Loads ESM2 embeddings, applies PCA (with and without whitening),
and computes silhouette scores for both KMeans and root HOG labels.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
DEFAULT_DATA_PATH = PROJECT_ROOT / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "chapters/chapter2/results/whitening"

def load_data(data_path, sample_size=5000, seed=42):
    print(f"Loading merged data (memory-efficient columns)...")
    import pyarrow.feather as feather
    import numpy as np
    import pandas as pd
    
    # Just read the schema to get column names
    schema = feather.read_table(data_path, columns=[]).schema
    emb_cols = [n for n in schema.names if n.startswith("ME:")]
    keep_cols = ["EntryID", "roothog_id"] + emb_cols
    
    # Read the data - we only read the columns we need.
    # On this machine, reading 1280 columns for all rows might be too much if it's 600k rows.
    # Let's try to read it and sample immediately.
    print(f"  Reading {len(keep_cols)} columns...")
    table = feather.read_table(data_path, columns=keep_cols)
    df = table.to_pandas()
    
    # Sample after loading if it fits, or we need a different format (like Parquet) for row sampling.
    if len(df) > sample_size:
        df = df.sample(n=sample_size, random_state=seed)
        
    return df, emb_cols

def run_analysis():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-size', type=int, default=5000)
    parser.add_argument('--n-clusters', type=int, default=20)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    output_dir = DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    df, emb_cols = load_data(DEFAULT_DATA_PATH, sample_size=args.sample_size, seed=args.seed)
    
    # Drop NaNs in embeddings
    df = df.dropna(subset=emb_cols)
    
    sample_df = df
    X = sample_df[emb_cols].values
    hog_labels = sample_df['roothog_id'].values
    
    # Mask for valid HOGs (non-zero)
    valid_hog_mask = (hog_labels != 0) & pd.notna(hog_labels)
    X_hog = X[valid_hog_mask]
    y_hog = hog_labels[valid_hog_mask].astype(int)

    results = []

    configs = [
        {"name": "Original", "pca": False, "whiten": False},
        {"name": "Standardized", "pca": False, "whiten": False, "scale": True},
        {"name": "PCA-50", "pca": True, "n_comp": 50, "whiten": False},
        {"name": "PCA-50-Whitened", "pca": True, "n_comp": 50, "whiten": True},
        {"name": "PCA-100", "pca": True, "n_comp": 100, "whiten": False},
        {"name": "PCA-100-Whitened", "pca": True, "n_comp": 100, "whiten": True},
        {"name": "PCA-256", "pca": True, "n_comp": 256, "whiten": False},
        {"name": "PCA-256-Whitened", "pca": True, "n_comp": 256, "whiten": True},
    ]

    for cfg in configs:
        print(f"\nProcessing: {cfg['name']}")
        X_proc = X.copy()
        X_hog_proc = X_hog.copy()

        if cfg.get("scale"):
            scaler = StandardScaler()
            X_proc = scaler.fit_transform(X_proc)
            X_hog_proc = scaler.transform(X_hog_proc)

        if cfg["pca"]:
            pca = PCA(n_components=cfg["n_comp"], whiten=cfg["whiten"], random_state=args.seed)
            X_proc = pca.fit_transform(X_proc)
            X_hog_proc = pca.transform(X_hog_proc)

        # KMeans Silhouette
        kmeans = KMeans(n_clusters=args.n_clusters, random_state=args.seed, n_init=10)
        km_labels = kmeans.fit_predict(X_proc)
        km_sil = silhouette_score(X_proc, km_labels, metric='cosine' if not cfg['whiten'] else 'euclidean')

        # HOG Silhouette
        hog_sil = silhouette_score(X_hog_proc, y_hog, metric='cosine' if not cfg['whiten'] else 'euclidean')

        res = {
            "config": cfg["name"],
            "kmeans_silhouette": float(km_sil),
            "hog_silhouette": float(hog_sil)
        }
        results.append(res)
        print(f"  KMeans Sil: {km_sil:.4f}")
        print(f"  HOG Sil:    {hog_sil:.4f}")

    # Save results
    with open(output_dir / "whitening_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    # Output summary
    res_df = pd.DataFrame(results)
    print("\nResults Summary:")
    print(res_df.to_string(index=False))

if __name__ == "__main__":
    run_analysis()
