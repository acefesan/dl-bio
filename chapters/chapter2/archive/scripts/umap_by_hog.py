#!/usr/bin/env python3
"""UMAP analysis of protein embeddings grouped by Hierarchical Orthologous Groups (HOGs).

This script explores intra-HOG variance in ESM2 embeddings to understand how
well protein language model embeddings capture evolutionary relationships.

HOGs (Hierarchical Orthologous Groups) from OMA browser represent proteins
descended from a common ancestor. Root HOGs trace back to LUCA (Last Universal
Common Ancestor) level.

If embeddings capture evolutionary/functional relationships well, proteins
within the same HOG should cluster together despite taxonomic distance.

Usage:
    python umap_by_hog.py                    # Default analysis
    python umap_by_hog.py --min-hog-size 50  # Include smaller HOGs
    python umap_by_hog.py --top-n 30         # Analyze top 30 HOGs
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from umap import UMAP

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
OUTPUT_DIR = BASE_PATH / "assets/proteins/analysis/hog_analysis"


def load_data() -> pd.DataFrame:
    """Load the merged CAFA3 dataset with embeddings and HOG annotations."""
    print("Loading dataset...")
    df = pd.read_feather(DATA_PATH)

    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    keep_cols = ['EntryID', 'taxonomyID', 'scientific_name', 'Length',
                 'roothog_id', 'hog_id', 'Sequence'] + embedding_cols

    # Get unique proteins
    protein_df = df[keep_cols].drop_duplicates(subset='EntryID')

    # Filter to proteins with embeddings and valid root HOGs
    protein_df = protein_df[protein_df[embedding_cols[0]].notna()]
    protein_df = protein_df[protein_df['roothog_id'].notna()]
    protein_df = protein_df[protein_df['roothog_id'] != 0]  # Exclude unassigned

    protein_df['roothog_id'] = protein_df['roothog_id'].astype(int)

    print(f"  Total proteins with embeddings: {len(protein_df):,}")
    print(f"  Unique root HOGs: {protein_df['roothog_id'].nunique():,}")

    return protein_df, embedding_cols


def get_abundant_hogs(protein_df: pd.DataFrame, min_size: int = 100) -> pd.DataFrame:
    """Get root HOGs with at least min_size proteins."""
    hog_counts = protein_df.groupby('roothog_id').agg(
        n_proteins=('EntryID', 'count'),
        n_taxa=('taxonomyID', 'nunique'),
        mean_length=('Length', 'mean')
    ).sort_values('n_proteins', ascending=False)

    abundant = hog_counts[hog_counts['n_proteins'] >= min_size]

    print(f"\nRoot HOGs with >= {min_size} proteins: {len(abundant)}")
    print(abundant.head(15))

    return abundant


def compute_hog_metrics(protein_df: pd.DataFrame, embedding_cols: list,
                        hog_ids: list) -> pd.DataFrame:
    """Compute embedding metrics for each HOG.

    Metrics:
    - centroid: Mean embedding of HOG members
    - intra_variance: Mean L2 distance from centroid (how spread out)
    - diameter: Max pairwise distance within HOG
    - cross_taxa_coherence: How well proteins cluster despite taxonomic distance
    """
    print("\nComputing HOG embedding metrics...")

    X_all = protein_df[embedding_cols].values
    global_centroid = X_all.mean(axis=0)

    metrics = []
    for hog_id in hog_ids:
        mask = protein_df['roothog_id'] == hog_id
        X_hog = protein_df.loc[mask, embedding_cols].values
        taxa = protein_df.loc[mask, 'taxonomyID'].values

        if len(X_hog) < 2:
            continue

        # Centroid
        centroid = X_hog.mean(axis=0)

        # Distance from global centroid
        dist_from_global = np.linalg.norm(centroid - global_centroid)

        # Intra-HOG variance (mean distance from centroid)
        dists_from_centroid = np.linalg.norm(X_hog - centroid, axis=1)
        intra_variance = dists_from_centroid.mean()
        intra_std = dists_from_centroid.std()

        # Diameter (max pairwise distance) - sample if too large
        if len(X_hog) > 200:
            idx = np.random.choice(len(X_hog), 200, replace=False)
            X_sample = X_hog[idx]
        else:
            X_sample = X_hog

        from scipy.spatial.distance import pdist
        pairwise_dists = pdist(X_sample)
        diameter = pairwise_dists.max()
        median_pairwise = np.median(pairwise_dists)

        # Cross-taxa coherence: ratio of intra-HOG to expected random variance
        # Lower is better (HOG is tighter than random)
        n_taxa = len(np.unique(taxa))

        metrics.append({
            'roothog_id': hog_id,
            'n_proteins': len(X_hog),
            'n_taxa': n_taxa,
            'dist_from_global': dist_from_global,
            'intra_variance': intra_variance,
            'intra_std': intra_std,
            'diameter': diameter,
            'median_pairwise': median_pairwise,
            'coherence': intra_variance / dist_from_global if dist_from_global > 0 else np.inf,
        })

    return pd.DataFrame(metrics)


def sample_hog_proteins(protein_df: pd.DataFrame, hog_ids: list,
                        max_per_hog: int = 100, seed: int = 42) -> pd.DataFrame:
    """Sample proteins from each HOG for visualization."""
    np.random.seed(seed)

    sampled = []
    for hog_id in hog_ids:
        hog_df = protein_df[protein_df['roothog_id'] == hog_id]
        n_sample = min(max_per_hog, len(hog_df))
        sampled.append(hog_df.sample(n=n_sample, random_state=seed))

    return pd.concat(sampled, ignore_index=True)


def run_umap(X: np.ndarray, **kwargs) -> np.ndarray:
    """Run UMAP with default parameters."""
    defaults = {
        'n_neighbors': 15,
        'min_dist': 0.1,
        'metric': 'cosine',
        'random_state': 42,
        'verbose': False
    }
    defaults.update(kwargs)

    print(f"Running UMAP on {X.shape[0]} samples...")
    umap_model = UMAP(**defaults)
    return umap_model.fit_transform(X)


def plot_hog_umap(sample_df: pd.DataFrame, Y: np.ndarray,
                  top_n: int = 10, output_path: Path = None):
    """Create UMAP visualization colored by HOG."""
    fig, ax = plt.subplots(figsize=(14, 10))

    # Get top HOGs by sample count
    hog_counts = sample_df['roothog_id'].value_counts()
    top_hogs = hog_counts.head(top_n).index.tolist()

    # Plot "other" HOGs as background
    other_mask = ~sample_df['roothog_id'].isin(top_hogs)
    ax.scatter(Y[other_mask, 0], Y[other_mask, 1],
               c='lightgray', s=10, alpha=0.3, label='Other HOGs')

    # Plot top HOGs with distinct colors
    colors = plt.cm.tab20.colors
    for i, hog_id in enumerate(top_hogs):
        mask = sample_df['roothog_id'] == hog_id
        n = mask.sum()
        ax.scatter(Y[mask, 0], Y[mask, 1],
                   c=[colors[i % len(colors)]], s=25, alpha=0.7,
                   label=f'HOG {hog_id} (n={n})')

    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('ESM2 Embeddings Colored by Root HOG\n(Orthologous proteins should cluster together)',
                 fontsize=14)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    plt.close()


def plot_hog_taxa_distribution(sample_df: pd.DataFrame, Y: np.ndarray,
                               hog_id: int, output_path: Path = None):
    """Show how one HOG distributes across taxa in UMAP space."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    hog_mask = sample_df['roothog_id'] == hog_id
    hog_df = sample_df[hog_mask]

    # Left: HOG highlighted in full UMAP
    ax = axes[0]
    ax.scatter(Y[~hog_mask, 0], Y[~hog_mask, 1], c='lightgray', s=5, alpha=0.2)
    ax.scatter(Y[hog_mask, 0], Y[hog_mask, 1], c='red', s=30, alpha=0.7)
    ax.set_title(f'HOG {hog_id} in UMAP Space (n={hog_mask.sum()})', fontsize=12)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    # Right: Same HOG colored by taxonomy
    ax = axes[1]
    taxa = hog_df['taxonomyID'].values
    unique_taxa = np.unique(taxa)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_taxa)))
    taxa_to_color = dict(zip(unique_taxa, colors))

    for taxid in unique_taxa:
        tax_mask = hog_df['taxonomyID'] == taxid
        idx = hog_df[tax_mask].index
        # Get positions in Y array
        pos = [sample_df.index.get_loc(i) for i in idx]
        name = hog_df[tax_mask]['scientific_name'].iloc[0].split()[0][:8]
        ax.scatter(Y[pos, 0], Y[pos, 1], c=[taxa_to_color[taxid]],
                   s=40, alpha=0.8, label=f'{name} (n={tax_mask.sum()})')

    ax.set_title(f'HOG {hog_id} by Taxonomy ({len(unique_taxa)} taxa)', fontsize=12)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    plt.close()


def plot_intra_hog_variance(metrics_df: pd.DataFrame, output_path: Path = None):
    """Visualize intra-HOG variance metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Variance vs HOG size
    ax = axes[0, 0]
    ax.scatter(metrics_df['n_proteins'], metrics_df['intra_variance'],
               c=metrics_df['n_taxa'], cmap='viridis', alpha=0.6)
    ax.set_xlabel('HOG Size (n proteins)')
    ax.set_ylabel('Intra-HOG Variance')
    ax.set_title('Variance vs HOG Size\n(color = n taxa)')
    ax.set_xscale('log')
    plt.colorbar(ax.collections[0], ax=ax, label='N Taxa')

    # 2. Variance vs n_taxa
    ax = axes[0, 1]
    ax.scatter(metrics_df['n_taxa'], metrics_df['intra_variance'],
               c=metrics_df['n_proteins'], cmap='plasma', alpha=0.6)
    ax.set_xlabel('Number of Taxa in HOG')
    ax.set_ylabel('Intra-HOG Variance')
    ax.set_title('Variance vs Taxonomic Breadth\n(color = HOG size)')
    plt.colorbar(ax.collections[0], ax=ax, label='N Proteins')

    # 3. Distribution of variance
    ax = axes[1, 0]
    ax.hist(metrics_df['intra_variance'], bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(metrics_df['intra_variance'].median(), color='red',
               linestyle='--', label=f"Median: {metrics_df['intra_variance'].median():.2f}")
    ax.set_xlabel('Intra-HOG Variance')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Intra-HOG Variance')
    ax.legend()

    # 4. Coherence (variance / distance from global)
    ax = axes[1, 1]
    coherent = metrics_df[metrics_df['coherence'] < np.inf].copy()
    ax.scatter(coherent['dist_from_global'], coherent['intra_variance'],
               c=coherent['n_taxa'], cmap='viridis', alpha=0.6)

    # Add diagonal reference line (coherence = 1)
    max_val = max(coherent['dist_from_global'].max(), coherent['intra_variance'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Coherence = 1')

    ax.set_xlabel('Distance from Global Centroid')
    ax.set_ylabel('Intra-HOG Variance')
    ax.set_title('Coherence: Variance vs Distance\n(below line = coherent, above = dispersed)')
    ax.legend()

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="UMAP analysis by HOG")
    parser.add_argument('--min-hog-size', type=int, default=100,
                        help='Minimum proteins per HOG (default: 100)')
    parser.add_argument('--top-n', type=int, default=20,
                        help='Number of top HOGs to analyze (default: 20)')
    parser.add_argument('--max-per-hog', type=int, default=100,
                        help='Max proteins to sample per HOG (default: 100)')
    args = parser.parse_args()

    # Setup output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    protein_df, embedding_cols = load_data()

    # Get abundant HOGs
    abundant_hogs = get_abundant_hogs(protein_df, min_size=args.min_hog_size)
    top_hog_ids = abundant_hogs.head(args.top_n).index.tolist()

    # Compute metrics for abundant HOGs
    metrics_df = compute_hog_metrics(protein_df, embedding_cols, top_hog_ids)
    metrics_df = metrics_df.merge(
        abundant_hogs.reset_index()[['roothog_id', 'n_proteins', 'n_taxa', 'mean_length']],
        on='roothog_id', suffixes=('', '_full')
    )

    print("\nHOG Metrics Summary:")
    print(metrics_df[['roothog_id', 'n_proteins', 'n_taxa', 'intra_variance',
                      'dist_from_global', 'coherence']].to_string(index=False))

    # Save metrics
    metrics_df.to_csv(OUTPUT_DIR / 'hog_metrics.csv', index=False)
    print(f"\nSaved metrics to: {OUTPUT_DIR / 'hog_metrics.csv'}")

    # Sample proteins for UMAP
    sample_df = sample_hog_proteins(protein_df, top_hog_ids,
                                    max_per_hog=args.max_per_hog)
    print(f"\nSampled {len(sample_df)} proteins from {len(top_hog_ids)} HOGs")

    # Run UMAP
    X = sample_df[embedding_cols].values
    Y = run_umap(X)

    sample_df['umap_x'] = Y[:, 0]
    sample_df['umap_y'] = Y[:, 1]

    # Save coordinates
    sample_df[['EntryID', 'roothog_id', 'taxonomyID', 'umap_x', 'umap_y']].to_csv(
        OUTPUT_DIR / 'hog_umap_coordinates.csv', index=False)

    # Generate plots
    print("\nGenerating visualizations...")

    plot_hog_umap(sample_df, Y, top_n=15,
                  output_path=OUTPUT_DIR / 'umap_by_hog.png')

    plot_intra_hog_variance(metrics_df,
                            output_path=OUTPUT_DIR / 'intra_hog_variance.png')

    # Detailed view of top 3 most coherent HOGs
    coherent_hogs = metrics_df.nsmallest(3, 'coherence')['roothog_id'].tolist()
    for hog_id in coherent_hogs:
        plot_hog_taxa_distribution(sample_df, Y, hog_id,
                                   output_path=OUTPUT_DIR / f'hog_{hog_id}_taxa.png')

    # Summary statistics
    print("\n" + "="*70)
    print("INTRA-HOG VARIANCE ANALYSIS SUMMARY")
    print("="*70)

    print(f"""
Dataset:
  Proteins with HOG assignments: {len(protein_df):,}
  Unique root HOGs: {protein_df['roothog_id'].nunique():,}
  HOGs analyzed (>= {args.min_hog_size} proteins): {len(metrics_df)}

Intra-HOG Variance:
  Mean: {metrics_df['intra_variance'].mean():.3f}
  Median: {metrics_df['intra_variance'].median():.3f}
  Range: [{metrics_df['intra_variance'].min():.3f}, {metrics_df['intra_variance'].max():.3f}]

Coherence (variance / distance from global):
  Mean: {metrics_df['coherence'].mean():.3f}
  < 1 means HOG is tighter than expected by distance
  > 1 means HOG is more dispersed than expected

Correlation with taxonomic breadth:
  Variance vs N_taxa: r = {spearmanr(metrics_df['n_taxa'], metrics_df['intra_variance'])[0]:.3f}
  Variance vs HOG_size: r = {spearmanr(metrics_df['n_proteins'], metrics_df['intra_variance'])[0]:.3f}

Output directory: {OUTPUT_DIR}
""")

    # Save summary (convert numpy types to Python native for JSON serialization)
    summary = {
        'n_proteins_with_hog': int(len(protein_df)),
        'n_unique_hogs': int(protein_df['roothog_id'].nunique()),
        'n_hogs_analyzed': int(len(metrics_df)),
        'min_hog_size': int(args.min_hog_size),
        'variance_mean': float(metrics_df['intra_variance'].mean()),
        'variance_median': float(metrics_df['intra_variance'].median()),
        'coherence_mean': float(metrics_df['coherence'].mean()),
        'corr_variance_ntaxa': float(spearmanr(metrics_df['n_taxa'], metrics_df['intra_variance'])[0]),
    }

    with open(OUTPUT_DIR / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print("Done!")


if __name__ == "__main__":
    main()
