#!/usr/bin/env python3
"""Clustering analysis of ESM2 protein embeddings.

Loads the merged CAFA3 dataset, performs KMeans clustering in the original
640-dimensional embedding space, computes silhouette scores, runs UMAP for
visualization, and generates plots colored by taxa and HOGs.

Usage:
    python 03_clustering_analysis.py                      # Default settings
    python 03_clustering_analysis.py --sample-size 1000   # More samples per taxon
    python 03_clustering_analysis.py --n-clusters 30      # More KMeans clusters
    python 03_clustering_analysis.py --output-dir ./out   # Custom output
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from umap import UMAP

# Paths relative to this script's location (chapters/chapter2/)
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_PATH = PROJECT_ROOT / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "chapters/chapter2/results"


# =============================================================================
# Data Loading & Preparation
# =============================================================================

def load_dataset(data_path: Path) -> pd.DataFrame:
    """Load the CAFA3 merged dataset."""
    print(f"Loading dataset from {data_path.name}...")
    df = pd.read_feather(data_path)
    print(f"  Loaded {len(df):,} rows, {df['EntryID'].nunique():,} unique proteins")
    return df


def get_unique_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate to unique proteins with embeddings."""
    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    protein_cols = ['EntryID', 'Length', 'taxonomyID', 'scientific_name',
                    'oma_id', 'hog_id', 'roothog_id'] + embedding_cols

    protein_df = df[protein_cols].drop_duplicates(subset='EntryID')
    protein_df = protein_df[protein_df[embedding_cols[0]].notna()]

    print(f"  Unique proteins with embeddings: {len(protein_df):,}")
    return protein_df


def get_abundant_taxa(protein_df: pd.DataFrame, min_count: int = 500) -> dict:
    """Get taxa with at least min_count proteins, mapped to short names."""
    taxa_counts = protein_df.groupby('taxonomyID').size().sort_values(ascending=False)
    abundant = taxa_counts[taxa_counts >= min_count]

    taxid_to_name = protein_df.groupby('taxonomyID')['scientific_name'].first().to_dict()

    taxa_map = {}
    for taxid, count in abundant.items():
        full_name = taxid_to_name.get(taxid)
        if full_name and isinstance(full_name, str):
            parts = full_name.split()
            short_name = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else full_name[:15]
        else:
            short_name = f"Tax_{taxid}"
        taxa_map[taxid] = short_name

    print(f"\nAbundant taxa (>= {min_count} proteins): {len(taxa_map)}")
    for taxid, name in list(taxa_map.items())[:10]:
        print(f"  {name}: {abundant.loc[taxid]:,}")
    if len(taxa_map) > 10:
        print(f"  ... and {len(taxa_map) - 10} more")

    return taxa_map


def get_abundant_hogs(protein_df: pd.DataFrame, min_count: int = 100) -> pd.Series:
    """Get root HOGs with at least min_count proteins."""
    hog_counts = protein_df['roothog_id'].value_counts()
    abundant = hog_counts[hog_counts >= min_count]
    abundant = abundant[abundant.index != 0]

    print(f"\nAbundant root HOGs (>= {min_count} proteins): {len(abundant)}")
    return abundant


def sample_proteins(
    protein_df: pd.DataFrame,
    taxa_map: dict,
    sample_size: int = 500,
    seed: int = 42
) -> pd.DataFrame:
    """Balanced sample of proteins from each taxon."""
    np.random.seed(seed)

    sampled_dfs = []
    for taxid in taxa_map:
        taxon_df = protein_df[protein_df['taxonomyID'] == taxid]
        n_sample = min(sample_size, len(taxon_df))
        sampled_dfs.append(taxon_df.sample(n=n_sample, random_state=seed))

    sample_df = pd.concat(sampled_dfs, ignore_index=True)
    sample_df['taxon_name'] = sample_df['taxonomyID'].map(taxa_map)

    print(f"\nSampled {len(sample_df):,} proteins from {len(taxa_map)} taxa")
    return sample_df


# =============================================================================
# Clustering & Metrics
# =============================================================================

def run_umap(X: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.1,
             seed: int = 42) -> np.ndarray:
    """Run UMAP dimensionality reduction."""
    print("\nRunning UMAP...")
    umap_model = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        metric='cosine',
        random_state=seed,
        verbose=True
    )
    return umap_model.fit_transform(X)


def compute_clustering_metrics(
    X: np.ndarray,
    hog_labels: np.ndarray,
    n_clusters: int = 20,
    seed: int = 42,
) -> dict:
    """Compute KMeans clustering and silhouette scores in original embedding space.

    Returns metrics dict with KMeans labels included.
    """
    dim = X.shape[1]
    print(f"\nRunning KMeans (k={n_clusters}) in {X.shape[1]}D space...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    kmeans_labels = kmeans.fit_predict(X)

    print("Computing silhouette scores (cosine metric, original space)...")

    # Silhouette for KMeans clusters
    kmeans_silhouette = silhouette_score(X, kmeans_labels, metric='cosine')

    # Silhouette for HOG labels (only if enough valid labels)
    valid_hog_mask = pd.notna(hog_labels) & (hog_labels != 0)
    hog_silhouette = None
    if valid_hog_mask.sum() > 100:
        X_valid = X[valid_hog_mask]
        hog_valid = hog_labels[valid_hog_mask].astype(int)
        n_unique = len(np.unique(hog_valid))
        if 2 <= n_unique < len(X_valid):
            hog_silhouette = silhouette_score(X_valid, hog_valid, metric='cosine')

    metrics = {
        "n_proteins": int(X.shape[0]),
        "embedding_dim": int(X.shape[1]),
        "n_kmeans_clusters": n_clusters,
        "kmeans_silhouette": float(kmeans_silhouette),
        "hog_silhouette": float(hog_silhouette) if hog_silhouette is not None else None,
        "n_proteins_with_hog": int(valid_hog_mask.sum()),
        "n_unique_hogs": int(len(np.unique(hog_labels[valid_hog_mask].astype(int)))) if valid_hog_mask.sum() > 0 else 0,
    }

    print(f"  KMeans silhouette ({dim}D, cosine): {kmeans_silhouette:.4f}")
    if hog_silhouette is not None:
        print(f"  HOG silhouette ({dim}D, cosine):    {hog_silhouette:.4f}")

    return metrics, kmeans_labels


# =============================================================================
# Visualization
# =============================================================================

def plot_by_taxa(sample_df: pd.DataFrame, output_path: Path):
    """UMAP plot colored by taxa."""
    fig, ax = plt.subplots(figsize=(14, 10))

    taxa_order = sample_df['taxon_name'].value_counts().index.tolist()
    colors = plt.cm.tab20.colors + plt.cm.tab20b.colors

    for i, taxon in enumerate(taxa_order):
        group = sample_df[sample_df['taxon_name'] == taxon]
        ax.scatter(
            group['umap_x'], group['umap_y'],
            label=f"{taxon} (n={len(group)})",
            alpha=0.6, s=15,
            c=[colors[i % len(colors)]]
        )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('ESM2 Protein Embeddings - UMAP by Taxa', fontsize=14)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_by_hog(sample_df: pd.DataFrame, abundant_hogs: pd.Series, output_path: Path):
    """UMAP plot colored by top root HOGs."""
    fig, ax = plt.subplots(figsize=(14, 10))

    top_hog_ids = abundant_hogs.index.tolist()[:15]

    other_mask = ~sample_df['roothog_id'].isin(top_hog_ids)
    other_df = sample_df[other_mask]
    ax.scatter(
        other_df['umap_x'], other_df['umap_y'],
        label=f'Other/No HOG (n={len(other_df)})',
        alpha=0.2, s=10, c='lightgray'
    )

    colors = plt.cm.tab20.colors
    for i, hog_id in enumerate(top_hog_ids):
        hog_df = sample_df[sample_df['roothog_id'] == hog_id]
        if len(hog_df) > 0:
            ax.scatter(
                hog_df['umap_x'], hog_df['umap_y'],
                label=f'HOG {int(hog_id)} (n={len(hog_df)})',
                alpha=0.7, s=25,
                c=[colors[i % len(colors)]]
            )

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('ESM2 Protein Embeddings - UMAP by Root HOG', fontsize=14)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_taxa_hog_grid(sample_df: pd.DataFrame, abundant_hogs: pd.Series, output_path: Path):
    """Heatmap showing HOG distribution across taxa."""
    top_taxa = sample_df['taxon_name'].value_counts().head(12).index.tolist()
    top_hogs = abundant_hogs.head(12).index.tolist()

    subset = sample_df[
        sample_df['taxon_name'].isin(top_taxa) &
        sample_df['roothog_id'].isin(top_hogs)
    ]

    crosstab = pd.crosstab(subset['taxon_name'], subset['roothog_id'])
    crosstab = crosstab.reindex(index=top_taxa, columns=top_hogs, fill_value=0)
    crosstab.columns = [f'HOG {int(h)}' for h in crosstab.columns]

    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(crosstab, annot=True, fmt='d', cmap='YlOrRd', ax=ax)
    ax.set_title('HOG Distribution Across Taxa (sampled proteins)', fontsize=14)
    ax.set_xlabel('Root HOG ID', fontsize=12)
    ax.set_ylabel('Species', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Clustering analysis of ESM2 protein embeddings")
    parser.add_argument('--sample-size', type=int, default=500,
                        help='Proteins to sample per taxon (default: 500)')
    parser.add_argument('--min-taxa', type=int, default=500,
                        help='Minimum proteins for taxon inclusion (default: 500)')
    parser.add_argument('--min-hog', type=int, default=100,
                        help='Minimum proteins for HOG inclusion (default: 100)')
    parser.add_argument('--n-clusters', type=int, default=20,
                        help='Number of KMeans clusters (default: 20)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for plots and metrics')
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load and prepare data
    df = load_dataset(DATA_PATH)
    protein_df = get_unique_proteins(df)

    taxa_map = get_abundant_taxa(protein_df, min_count=args.min_taxa)
    abundant_hogs = get_abundant_hogs(protein_df, min_count=args.min_hog)

    # Sample proteins
    sample_df = sample_proteins(protein_df, taxa_map, sample_size=args.sample_size)

    # Extract embeddings
    embedding_cols = [c for c in sample_df.columns if c.startswith('ME:')]
    X = sample_df[embedding_cols].values
    print(f"\nEmbedding matrix: {X.shape}")

    # Clustering metrics in original embedding space
    hog_labels = sample_df['roothog_id'].values
    metrics, kmeans_labels = compute_clustering_metrics(
        X, hog_labels, n_clusters=args.n_clusters
    )
    sample_df['kmeans_cluster'] = kmeans_labels

    # UMAP for visualization
    umap_coords = run_umap(X)
    sample_df['umap_x'] = umap_coords[:, 0]
    sample_df['umap_y'] = umap_coords[:, 1]

    # Generate plots
    print("\nGenerating plots...")
    plot_by_taxa(sample_df, output_dir / 'umap_by_taxa.png')
    plot_by_hog(sample_df, abundant_hogs, output_dir / 'umap_by_hog.png')
    plot_taxa_hog_grid(sample_df, abundant_hogs, output_dir / 'taxa_hog_heatmap.png')

    # Save UMAP coordinates
    coords_df = sample_df[['EntryID', 'taxonomyID', 'taxon_name', 'roothog_id',
                            'kmeans_cluster', 'umap_x', 'umap_y']]
    coords_df.to_csv(output_dir / 'umap_coordinates.csv', index=False)
    print(f"Saved UMAP coordinates: {output_dir / 'umap_coordinates.csv'}")

    # Save metrics
    with open(output_dir / 'clustering_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved metrics: {output_dir / 'clustering_metrics.json'}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Taxa analyzed: {len(taxa_map)}")
    print(f"Proteins sampled: {len(sample_df):,}")
    print(f"KMeans clusters: {args.n_clusters}")
    dim = metrics['embedding_dim']
    print(f"KMeans silhouette ({dim}D): {metrics['kmeans_silhouette']:.4f}")
    if metrics['hog_silhouette'] is not None:
        print(f"HOG silhouette ({dim}D):    {metrics['hog_silhouette']:.4f}")
    print(f"Output directory: {output_dir}")
    print("\nDone!")


if __name__ == "__main__":
    main()
