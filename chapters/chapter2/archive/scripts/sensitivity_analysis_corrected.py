#!/usr/bin/env python3
"""CORRECTED Sensitivity Analysis with Fixed Metrics.

This script fixes the methodological issues identified in CRITICAL_REVIEW.md:

1. **ARI Bug Fixed**: Compares KMeans cluster assignments (not HOG labels)
2. **Original Space Metrics**: Silhouette computed in 640D embedding space
3. **Proper Methodology**: Tests actual clustering stability, not label identity

Key Fixes:
- ARI now compares cluster assignments between different samples
- Silhouette scores computed in original 640D space using cosine distance
- Added comparison of both UMAP-space and original-space metrics

Usage:
    python sensitivity_analysis_corrected.py
"""

import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score, calinski_harabasz_score
from umap import UMAP

warnings.filterwarnings('ignore')

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
OUTPUT_DIR = BASE_PATH / "assets/proteins/analysis/sensitivity_analysis_corrected"


def load_data() -> tuple[pd.DataFrame, list]:
    """Load the merged CAFA3 dataset with embeddings and HOG annotations."""
    print("Loading dataset (this may take a while for 1.8GB file)...")
    df = pd.read_feather(DATA_PATH)
    
    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    keep_cols = ['EntryID', 'taxonomyID', 'scientific_name', 'Length',
                 'roothog_id', 'hog_id', 'Sequence'] + embedding_cols
    
    protein_df = df[keep_cols].drop_duplicates(subset='EntryID')
    protein_df = protein_df[protein_df[embedding_cols[0]].notna()]
    protein_df = protein_df[protein_df['roothog_id'].notna()]
    protein_df = protein_df[protein_df['roothog_id'] != 0]
    
    protein_df['roothog_id'] = protein_df['roothog_id'].astype(int)
    protein_df['taxonomyID'] = protein_df['taxonomyID'].astype(int)
    
    print(f"  Total proteins with embeddings: {len(protein_df):,}")
    print(f"  Unique root HOGs: {protein_df['roothog_id'].nunique():,}")
    print(f"  Unique taxa: {protein_df['taxonomyID'].nunique():,}")
    
    return protein_df, embedding_cols


def run_umap(X: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.1,
             random_state: int = 42) -> np.ndarray:
    """Run UMAP dimensionality reduction."""
    umap_model = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric='cosine',
        random_state=random_state,
        verbose=False
    )
    return umap_model.fit_transform(X)


def sample_hog_proteins(protein_df: pd.DataFrame, hog_ids: list,
                        max_per_hog: int = 100, seed: int = 42) -> pd.DataFrame:
    """Sample proteins from each HOG."""
    np.random.seed(seed)
    sampled = []
    for hog_id in hog_ids:
        hog_df = protein_df[protein_df['roothog_id'] == hog_id]
        n_sample = min(max_per_hog, len(hog_df))
        sampled.append(hog_df.sample(n=n_sample, random_state=seed))
    return pd.concat(sampled, ignore_index=True)


def compute_cluster_metrics_both_spaces(X_orig: np.ndarray, X_umap: np.ndarray, 
                                        labels: np.ndarray) -> dict:
    """
    Compute clustering metrics in BOTH original and UMAP space.
    
    This demonstrates the difference between the two approaches.
    """
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    if n_clusters < 2 or n_clusters >= len(X_orig):
        return {
            'n_clusters': n_clusters,
            'silhouette_original': np.nan,
            'silhouette_umap': np.nan,
            'calinski_harabasz_original': np.nan,
            'calinski_harabasz_umap': np.nan
        }
    
    return {
        'n_clusters': n_clusters,
        'silhouette_original': silhouette_score(X_orig, labels, metric='cosine'),
        'silhouette_umap': silhouette_score(X_umap, labels, metric='euclidean'),
        'calinski_harabasz_original': calinski_harabasz_score(X_orig, labels),
        'calinski_harabasz_umap': calinski_harabasz_score(X_umap, labels)
    }


def sensitivity_sampling_stability_corrected(
    protein_df: pd.DataFrame, embedding_cols: list,
    n_iterations: int = 10, n_hogs_per_sample: int = 500,
    min_hog_size: int = 5, max_per_hog: int = 50
) -> dict:
    """
    CORRECTED: Test sampling stability with proper ARI computation.
    
    Key fixes:
    1. ARI compares CLUSTER ASSIGNMENTS (not HOG labels)
    2. Silhouette computed in ORIGINAL 640D space
    3. Both original-space and UMAP-space metrics reported
    
    Question: Do different random samples produce consistent cluster structures?
    """
    print("\n" + "="*70)
    print("CORRECTED SENSITIVITY ANALYSIS: Sampling Stability")
    print("="*70)
    print("\nKey methodological improvements:")
    print("  ✓ ARI compares cluster assignments (not HOG labels)")
    print("  ✓ Silhouette computed in original 640D embedding space")
    print("  ✓ UMAP-space metrics also reported for comparison")
    
    # Get eligible HOGs
    hog_counts = protein_df['roothog_id'].value_counts()
    eligible_hogs = hog_counts[hog_counts >= min_hog_size].index.tolist()
    print(f"\nEligible HOGs (>= {min_hog_size} proteins): {len(eligible_hogs)}")
    
    n_hogs_per_sample = min(n_hogs_per_sample, len(eligible_hogs))
    print(f"Sampling {n_hogs_per_sample} HOGs per iteration")
    
    results = {
        'n_iterations': n_iterations,
        'n_hogs_per_sample': n_hogs_per_sample,
        'min_hog_size': min_hog_size,
        'methodology': 'Corrected: ARI on cluster assignments, silhouette in original space',
        'iterations': [],
        'pairwise_ari': []
    }
    
    all_samples = []
    all_embeddings = []  # Store original embeddings
    all_umaps = []
    all_kmeans_labels = []
    
    for i in range(n_iterations):
        print(f"\nIteration {i+1}/{n_iterations}...")
        seed = 42 + i * 100
        np.random.seed(seed)
        
        # Sample random HOGs
        sampled_hogs = np.random.choice(eligible_hogs, n_hogs_per_sample, replace=False)
        sample_df = sample_hog_proteins(protein_df, sampled_hogs, max_per_hog, seed)
        print(f"  Sampled {len(sample_df)} proteins from {n_hogs_per_sample} HOGs")
        
        # Get embeddings (ORIGINAL space)
        X = sample_df[embedding_cols].values
        
        # Run UMAP
        Y = run_umap(X, random_state=seed)
        
        # Cluster in UMAP space
        n_clusters = min(20, n_hogs_per_sample // 10)
        kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
        cluster_labels = kmeans.fit_predict(Y)
        
        # Compute metrics in BOTH spaces
        hog_labels = sample_df['roothog_id'].values
        
        # HOG-based metrics (comparing to ground truth HOG labels)
        hog_metrics_orig = silhouette_score(X, hog_labels, metric='cosine')
        hog_metrics_umap = silhouette_score(Y, hog_labels, metric='euclidean')
        
        # KMeans cluster metrics
        kmeans_metrics_orig = silhouette_score(X, cluster_labels, metric='cosine')
        kmeans_metrics_umap = silhouette_score(Y, cluster_labels, metric='euclidean')
        
        iteration_result = {
            'seed': seed,
            'n_proteins': len(sample_df),
            'n_unique_hogs': sample_df['roothog_id'].nunique(),
            'n_unique_taxa': sample_df['taxonomyID'].nunique(),
            'n_clusters': n_clusters,
            # HOG-based silhouette (in both spaces)
            'hog_silhouette_original': float(hog_metrics_orig),
            'hog_silhouette_umap': float(hog_metrics_umap),
            # KMeans cluster silhouette (in both spaces)
            'kmeans_silhouette_original': float(kmeans_metrics_orig),
            'kmeans_silhouette_umap': float(kmeans_metrics_umap),
        }
        results['iterations'].append(iteration_result)
        
        all_samples.append(sample_df)
        all_embeddings.append(X)
        all_umaps.append(Y)
        all_kmeans_labels.append(cluster_labels)
        
        print(f"  HOG Silhouette (640D): {hog_metrics_orig:.4f}")
        print(f"  HOG Silhouette (UMAP): {hog_metrics_umap:.4f}")
        print(f"  KMeans Silhouette (640D): {kmeans_metrics_orig:.4f}")
        print(f"  KMeans Silhouette (UMAP): {kmeans_metrics_umap:.4f}")
    
    # CORRECTED ARI COMPUTATION: Compare CLUSTER ASSIGNMENTS
    print("\n" + "="*70)
    print("COMPUTING CORRECTED ARI (Cluster Assignments, Not HOG Labels)")
    print("="*70)
    
    ari_matrix = np.zeros((n_iterations, n_iterations))
    
    for i in range(n_iterations):
        for j in range(i+1, n_iterations):
            # Find overlapping proteins
            entryids_i = set(all_samples[i]['EntryID'])
            entryids_j = set(all_samples[j]['EntryID'])
            overlap = entryids_i & entryids_j
            
            if len(overlap) > 50:  # Only compute if sufficient overlap
                # Add cluster labels to full samples first
                df_i = all_samples[i].copy()
                df_j = all_samples[j].copy()
                df_i['cluster'] = all_kmeans_labels[i]
                df_j['cluster'] = all_kmeans_labels[j]
                
                # Filter to overlapping proteins
                df_i_overlap = df_i[df_i['EntryID'].isin(overlap)][['EntryID', 'cluster']].set_index('EntryID')
                df_j_overlap = df_j[df_j['EntryID'].isin(overlap)][['EntryID', 'cluster']].set_index('EntryID')
                
                # Align by EntryID
                merged = df_i_overlap.join(df_j_overlap, lsuffix='_i', rsuffix='_j')
                
                # CORRECT: Compare cluster assignments
                ari = adjusted_rand_score(merged['cluster_i'], merged['cluster_j'])
                ari_matrix[i, j] = ari
                ari_matrix[j, i] = ari
                results['pairwise_ari'].append({
                    'i': i, 
                    'j': j, 
                    'ari': float(ari), 
                    'n_overlap': len(overlap)
                })
                
                print(f"  Samples {i} vs {j}: ARI = {ari:.4f} ({len(overlap)} overlapping proteins)")
    
    np.fill_diagonal(ari_matrix, 1.0)
    
    # Summary statistics
    metrics_df = pd.DataFrame(results['iterations'])
    
    ari_values = [x['ari'] for x in results['pairwise_ari']]
    
    results['summary'] = {
        'hog_silhouette_original_mean': float(metrics_df['hog_silhouette_original'].mean()),
        'hog_silhouette_original_std': float(metrics_df['hog_silhouette_original'].std()),
        'hog_silhouette_umap_mean': float(metrics_df['hog_silhouette_umap'].mean()),
        'hog_silhouette_umap_std': float(metrics_df['hog_silhouette_umap'].std()),
        'kmeans_silhouette_original_mean': float(metrics_df['kmeans_silhouette_original'].mean()),
        'kmeans_silhouette_original_std': float(metrics_df['kmeans_silhouette_original'].std()),
        'kmeans_silhouette_umap_mean': float(metrics_df['kmeans_silhouette_umap'].mean()),
        'kmeans_silhouette_umap_std': float(metrics_df['kmeans_silhouette_umap'].std()),
        'ari_mean': float(np.mean(ari_values)) if ari_values else np.nan,
        'ari_std': float(np.std(ari_values)) if ari_values else np.nan,
        'ari_min': float(np.min(ari_values)) if ari_values else np.nan,
        'ari_max': float(np.max(ari_values)) if ari_values else np.nan,
    }
    
    print("\n" + "="*70)
    print("CORRECTED RESULTS SUMMARY")
    print("="*70)
    print("\nHOG-based Silhouette Scores:")
    print(f"  Original 640D space: {results['summary']['hog_silhouette_original_mean']:.4f} ± {results['summary']['hog_silhouette_original_std']:.4f}")
    print(f"  UMAP 2D space:       {results['summary']['hog_silhouette_umap_mean']:.4f} ± {results['summary']['hog_silhouette_umap_std']:.4f}")
    
    print("\nKMeans Cluster Silhouette Scores:")
    print(f"  Original 640D space: {results['summary']['kmeans_silhouette_original_mean']:.4f} ± {results['summary']['kmeans_silhouette_original_std']:.4f}")
    print(f"  UMAP 2D space:       {results['summary']['kmeans_silhouette_umap_mean']:.4f} ± {results['summary']['kmeans_silhouette_umap_std']:.4f}")
    
    print(f"\nAdjusted Rand Index (Cluster Assignment Stability):")
    print(f"  Mean: {results['summary']['ari_mean']:.4f} ± {results['summary']['ari_std']:.4f}")
    print(f"  Range: [{results['summary']['ari_min']:.4f}, {results['summary']['ari_max']:.4f}]")
    
    # Interpretation
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    if results['summary']['ari_mean'] > 0.5:
        print("✅ GOOD: Cluster assignments are stable across samples (ARI > 0.5)")
    elif results['summary']['ari_mean'] > 0.3:
        print("🟡 MODERATE: Some stability in cluster assignments (ARI 0.3-0.5)")
    else:
        print("❌ POOR: Cluster assignments are unstable (ARI < 0.3)")
    
    orig_hog_sil = results['summary']['hog_silhouette_original_mean']
    if orig_hog_sil > 0.5:
        print("✅ GOOD: HOGs form well-separated clusters in original space")
    elif orig_hog_sil > 0.25:
        print("🟡 MODERATE: HOGs form weak clusters in original space")
    else:
        print("❌ POOR: HOGs do not cluster in original space")
    
    # Store for visualization
    results['ari_matrix'] = ari_matrix.tolist()
    results['last_umap'] = all_umaps[-1]
    results['last_sample_df'] = all_samples[-1]
    results['last_kmeans_labels'] = all_kmeans_labels[-1]
    
    return results


def plot_corrected_results(results: dict, output_dir: Path):
    """Generate visualizations comparing corrected vs original metrics."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: ARI matrix heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    ari_matrix = np.array(results['ari_matrix'])
    
    sns.heatmap(ari_matrix, annot=True, fmt='.3f', cmap='RdYlGn',
                vmin=0, vmax=1, ax=ax, square=True)
    ax.set_title('Adjusted Rand Index Matrix\n(Cluster Assignment Similarity Between Samples)',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Sample Index')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'ari_matrix_corrected.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'ari_matrix_corrected.png'}")
    plt.close()
    
    # Plot 2: Silhouette comparison (Original vs UMAP space)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    iterations_df = pd.DataFrame(results['iterations'])
    
    # HOG silhouette
    ax = axes[0]
    x = np.arange(len(iterations_df))
    ax.plot(x, iterations_df['hog_silhouette_original'], 'o-', label='Original 640D', linewidth=2)
    ax.plot(x, iterations_df['hog_silhouette_umap'], 's--', label='UMAP 2D', linewidth=2)
    ax.axhline(results['summary']['hog_silhouette_original_mean'], 
               color='blue', linestyle=':', alpha=0.5, label='Mean (640D)')
    ax.axhline(results['summary']['hog_silhouette_umap_mean'], 
               color='orange', linestyle=':', alpha=0.5, label='Mean (UMAP)')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Silhouette Score')
    ax.set_title('HOG-Based Silhouette\n(Original vs UMAP Space)', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # KMeans silhouette
    ax = axes[1]
    ax.plot(x, iterations_df['kmeans_silhouette_original'], 'o-', label='Original 640D', linewidth=2)
    ax.plot(x, iterations_df['kmeans_silhouette_umap'], 's--', label='UMAP 2D', linewidth=2)
    ax.axhline(results['summary']['kmeans_silhouette_original_mean'], 
               color='blue', linestyle=':', alpha=0.5, label='Mean (640D)')
    ax.axhline(results['summary']['kmeans_silhouette_umap_mean'], 
               color='orange', linestyle=':', alpha=0.5, label='Mean (UMAP)')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Silhouette Score')
    ax.set_title('KMeans Cluster Silhouette\n(Original vs UMAP Space)', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'silhouette_comparison.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'silhouette_comparison.png'}")
    plt.close()
    
    # Plot 3: Summary comparison table
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('off')
    
    table_data = [
        ['Metric', 'Original 640D', 'UMAP 2D', 'Difference'],
        ['', '', '', ''],
        ['HOG Silhouette (mean)', 
         f"{results['summary']['hog_silhouette_original_mean']:.4f}",
         f"{results['summary']['hog_silhouette_umap_mean']:.4f}",
         f"{results['summary']['hog_silhouette_umap_mean'] - results['summary']['hog_silhouette_original_mean']:+.4f}"],
        ['HOG Silhouette (std)', 
         f"{results['summary']['hog_silhouette_original_std']:.4f}",
         f"{results['summary']['hog_silhouette_umap_std']:.4f}",
         f"{results['summary']['hog_silhouette_umap_std'] - results['summary']['hog_silhouette_original_std']:+.4f}"],
        ['', '', '', ''],
        ['KMeans Silhouette (mean)', 
         f"{results['summary']['kmeans_silhouette_original_mean']:.4f}",
         f"{results['summary']['kmeans_silhouette_umap_mean']:.4f}",
         f"{results['summary']['kmeans_silhouette_umap_mean'] - results['summary']['kmeans_silhouette_original_mean']:+.4f}"],
        ['KMeans Silhouette (std)', 
         f"{results['summary']['kmeans_silhouette_original_std']:.4f}",
         f"{results['summary']['kmeans_silhouette_umap_std']:.4f}",
         f"{results['summary']['kmeans_silhouette_umap_std'] - results['summary']['kmeans_silhouette_original_std']:+.4f}"],
        ['', '', '', ''],
        ['Cluster ARI (mean)', 
         f"{results['summary']['ari_mean']:.4f}", 
         'N/A', 
         ''],
        ['Cluster ARI (std)', 
         f"{results['summary']['ari_std']:.4f}", 
         'N/A', 
         ''],
    ]
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.4, 0.2, 0.2, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header
    for i in range(4):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax.set_title('Corrected Sensitivity Analysis: Metric Comparison',
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'metrics_summary_table.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'metrics_summary_table.png'}")
    plt.close()


def main():
    """Run corrected sensitivity analysis."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("CORRECTED SENSITIVITY ANALYSIS")
    print("="*70)
    print("\nThis analysis fixes the methodological issues in the original:")
    print("  1. ARI compares cluster assignments (not HOG labels)")
    print("  2. Silhouette computed in original 640D space")
    print("  3. Both original and UMAP metrics reported for comparison")
    
    # Load data
    protein_df, embedding_cols = load_data()
    
    # Run corrected sensitivity analysis
    results = sensitivity_sampling_stability_corrected(
        protein_df, embedding_cols,
        n_iterations=10,
        n_hogs_per_sample=500,
        min_hog_size=5,
        max_per_hog=50
    )
    
    # Generate visualizations
    print("\n" + "="*70)
    print("GENERATING VISUALIZATIONS")
    print("="*70)
    plot_corrected_results(results, OUTPUT_DIR)
    
    # Save results
    results_to_save = {k: v for k, v in results.items() 
                       if k not in ['last_umap', 'last_sample_df', 'last_kmeans_labels']}
    
    with open(OUTPUT_DIR / 'corrected_results.json', 'w') as f:
        json.dump(results_to_save, f, indent=2)
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nDone!")


if __name__ == "__main__":
    main()
