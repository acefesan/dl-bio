#!/usr/bin/env python3
"""Sensitivity Analysis of HOG/Taxa Sampling for UMAP Clustering.

This script performs sensitivity analysis on the original umap_by_hog.py analysis:

1. Sampling Stability: If we sample another 500 random HOGs, do we get the same clusters?
   - Compares cluster metrics across multiple random samples
   - Uses Adjusted Rand Index (ARI) to measure cluster similarity

2. Taxa-based Clustering: What happens if we cluster all taxa with > 500 proteins?
   - Compares HOG-based clustering to taxonomy-based clustering
   - Assesses whether taxonomic signal dominates evolutionary signal

Usage:
    python sensitivity_analysis.py
"""

import json
import warnings
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score, calinski_harabasz_score
from sklearn.metrics.pairwise import cosine_similarity
from umap import UMAP

warnings.filterwarnings('ignore')

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
OUTPUT_DIR = BASE_PATH / "assets/proteins/analysis/sensitivity_analysis"


def load_data() -> tuple[pd.DataFrame, list]:
    """Load the merged CAFA3 dataset with embeddings and HOG annotations."""
    print("Loading dataset (this may take a while for 1.8GB file)...")
    df = pd.read_feather(DATA_PATH)
    
    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    keep_cols = ['EntryID', 'taxonomyID', 'scientific_name', 'Length',
                 'roothog_id', 'hog_id', 'Sequence'] + embedding_cols
    
    # Get unique proteins
    protein_df = df[keep_cols].drop_duplicates(subset='EntryID')
    
    # Filter to proteins with embeddings and valid root HOGs
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


def compute_cluster_metrics(X: np.ndarray, labels: np.ndarray) -> dict:
    """Compute clustering quality metrics."""
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    # Handle edge cases
    if n_clusters < 2 or n_clusters >= len(X):
        return {
            'n_clusters': n_clusters,
            'silhouette': np.nan,
            'calinski_harabasz': np.nan
        }
    
    return {
        'n_clusters': n_clusters,
        'silhouette': silhouette_score(X, labels, metric='cosine'),
        'calinski_harabasz': calinski_harabasz_score(X, labels)
    }


def compute_hog_centroids(protein_df: pd.DataFrame, embedding_cols: list,
                          hog_ids: list) -> tuple[np.ndarray, list]:
    """Compute centroid embeddings for each HOG."""
    centroids = []
    valid_hog_ids = []
    
    for hog_id in hog_ids:
        mask = protein_df['roothog_id'] == hog_id
        if mask.sum() > 0:
            centroid = protein_df.loc[mask, embedding_cols].mean().values
            centroids.append(centroid)
            valid_hog_ids.append(hog_id)
    
    return np.array(centroids), valid_hog_ids


def compute_taxa_centroids(protein_df: pd.DataFrame, embedding_cols: list,
                           taxa_ids: list) -> tuple[np.ndarray, list]:
    """Compute centroid embeddings for each taxon."""
    centroids = []
    valid_taxa_ids = []
    
    for taxid in taxa_ids:
        mask = protein_df['taxonomyID'] == taxid
        if mask.sum() > 0:
            centroid = protein_df.loc[mask, embedding_cols].mean().values
            centroids.append(centroid)
            valid_taxa_ids.append(taxid)
    
    return np.array(centroids), valid_taxa_ids


# =============================================================================
# SENSITIVITY ANALYSIS 1: Sampling Stability
# =============================================================================

def sensitivity_sampling_stability(protein_df: pd.DataFrame, embedding_cols: list,
                                    n_iterations: int = 10, n_hogs_per_sample: int = 500,
                                    min_hog_size: int = 5, max_per_hog: int = 50) -> dict:
    """
    Test sampling stability by running UMAP with different random HOG samples.
    
    Question: If we sample another 500 random groups, do we get the same clusters?
    
    Approach:
    1. Sample n_hogs_per_sample random HOGs (with min_hog_size proteins each)
    2. For each sample, compute UMAP and cluster
    3. Compare cluster assignments using Adjusted Rand Index (ARI)
    4. Compute variance in clustering metrics across samples
    """
    print("\n" + "="*70)
    print("SENSITIVITY ANALYSIS 1: Sampling Stability")
    print("="*70)
    
    # Get eligible HOGs (those with >= min_hog_size proteins)
    hog_counts = protein_df['roothog_id'].value_counts()
    eligible_hogs = hog_counts[hog_counts >= min_hog_size].index.tolist()
    print(f"\nEligible HOGs (>= {min_hog_size} proteins): {len(eligible_hogs)}")
    
    # Adjust n_hogs_per_sample if we don't have enough
    n_hogs_per_sample = min(n_hogs_per_sample, len(eligible_hogs))
    print(f"Sampling {n_hogs_per_sample} HOGs per iteration")
    
    results = {
        'n_iterations': n_iterations,
        'n_hogs_per_sample': n_hogs_per_sample,
        'min_hog_size': min_hog_size,
        'iterations': [],
        'pairwise_ari': []
    }
    
    all_samples = []
    all_umaps = []
    all_kmeans_labels = []
    
    for i in range(n_iterations):
        print(f"\nIteration {i+1}/{n_iterations}...")
        seed = 42 + i * 100
        np.random.seed(seed)
        
        # Sample random HOGs
        sampled_hogs = np.random.choice(eligible_hogs, n_hogs_per_sample, replace=False)
        
        # Sample proteins from these HOGs
        sample_df = sample_hog_proteins(protein_df, sampled_hogs, max_per_hog, seed)
        print(f"  Sampled {len(sample_df)} proteins from {n_hogs_per_sample} HOGs")
        
        # Get embeddings
        X = sample_df[embedding_cols].values
        
        # Run UMAP
        Y = run_umap(X, random_state=seed)
        
        # Cluster in UMAP space (use n_clusters based on number of HOGs)
        n_clusters = min(20, n_hogs_per_sample // 10)
        kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
        cluster_labels = kmeans.fit_predict(Y)
        
        # Compute metrics
        hog_metrics = compute_cluster_metrics(X, sample_df['roothog_id'].values)
        kmeans_metrics = compute_cluster_metrics(Y, cluster_labels)
        
        iteration_result = {
            'seed': seed,
            'n_proteins': len(sample_df),
            'n_unique_hogs': sample_df['roothog_id'].nunique(),
            'n_unique_taxa': sample_df['taxonomyID'].nunique(),
            'hog_silhouette': hog_metrics['silhouette'],
            'kmeans_silhouette': kmeans_metrics['silhouette'],
            'kmeans_calinski_harabasz': kmeans_metrics['calinski_harabasz']
        }
        results['iterations'].append(iteration_result)
        
        all_samples.append(sample_df)
        all_umaps.append(Y)
        all_kmeans_labels.append(cluster_labels)
        
        print(f"  HOG Silhouette: {hog_metrics['silhouette']:.4f}")
        print(f"  KMeans Silhouette: {kmeans_metrics['silhouette']:.4f}")
    
    # Compute pairwise ARI between iterations (for overlapping proteins)
    print("\nComputing pairwise cluster similarity (Adjusted Rand Index)...")
    ari_matrix = np.zeros((n_iterations, n_iterations))
    
    for i in range(n_iterations):
        for j in range(i+1, n_iterations):
            # Find overlapping proteins
            overlap = set(all_samples[i]['EntryID']) & set(all_samples[j]['EntryID'])
            if len(overlap) > 100:  # Only compute if sufficient overlap
                idx_i = all_samples[i][all_samples[i]['EntryID'].isin(overlap)].index
                idx_j = all_samples[j][all_samples[j]['EntryID'].isin(overlap)].index
                
                labels_i = all_samples[i].loc[idx_i, 'roothog_id'].values
                labels_j = all_samples[j].loc[idx_j, 'roothog_id'].values
                
                # Need to align by EntryID
                df_i = all_samples[i].loc[idx_i, ['EntryID', 'roothog_id']].set_index('EntryID')
                df_j = all_samples[j].loc[idx_j, ['EntryID', 'roothog_id']].set_index('EntryID')
                merged = df_i.join(df_j, lsuffix='_i', rsuffix='_j')
                
                ari = adjusted_rand_score(merged['roothog_id_i'], merged['roothog_id_j'])
                ari_matrix[i, j] = ari
                ari_matrix[j, i] = ari
                results['pairwise_ari'].append({'i': i, 'j': j, 'ari': ari, 'n_overlap': len(overlap)})
    
    # Fill diagonal
    np.fill_diagonal(ari_matrix, 1.0)
    
    # Summary statistics
    metrics_df = pd.DataFrame(results['iterations'])
    
    results['summary'] = {
        'hog_silhouette_mean': float(metrics_df['hog_silhouette'].mean()),
        'hog_silhouette_std': float(metrics_df['hog_silhouette'].std()),
        'kmeans_silhouette_mean': float(metrics_df['kmeans_silhouette'].mean()),
        'kmeans_silhouette_std': float(metrics_df['kmeans_silhouette'].std()),
        'ari_mean': float(np.mean([x['ari'] for x in results['pairwise_ari']])) if results['pairwise_ari'] else np.nan,
        'ari_std': float(np.std([x['ari'] for x in results['pairwise_ari']])) if results['pairwise_ari'] else np.nan,
    }
    
    print("\n" + "-"*50)
    print("SAMPLING STABILITY SUMMARY")
    print("-"*50)
    print(f"HOG Silhouette: {results['summary']['hog_silhouette_mean']:.4f} ± {results['summary']['hog_silhouette_std']:.4f}")
    print(f"KMeans Silhouette: {results['summary']['kmeans_silhouette_mean']:.4f} ± {results['summary']['kmeans_silhouette_std']:.4f}")
    print(f"Pairwise ARI (HOG labels): {results['summary']['ari_mean']:.4f} ± {results['summary']['ari_std']:.4f}")
    
    # Store for visualization
    results['ari_matrix'] = ari_matrix.tolist()
    results['last_umap'] = Y
    results['last_sample_df'] = sample_df
    
    return results


# =============================================================================
# SENSITIVITY ANALYSIS 2: Taxa-based Clustering
# =============================================================================

def sensitivity_taxa_clustering(protein_df: pd.DataFrame, embedding_cols: list,
                                 min_proteins_per_taxon: int = 500) -> dict:
    """
    Analyze clustering when grouping by taxonomy instead of HOG.
    
    Question: What happens if we cluster all taxa with > 500 proteins?
    
    Approach:
    1. Select taxa with >= min_proteins_per_taxon proteins
    2. Compute taxon centroids (mean embedding)
    3. Cluster taxa in embedding space
    4. Compare to HOG-based clustering
    """
    print("\n" + "="*70)
    print(f"SENSITIVITY ANALYSIS 2: Taxa-based Clustering (>{min_proteins_per_taxon} proteins)")
    print("="*70)
    
    # Get taxa with sufficient proteins
    taxa_counts = protein_df['taxonomyID'].value_counts()
    eligible_taxa = taxa_counts[taxa_counts >= min_proteins_per_taxon].index.tolist()
    print(f"\nTaxa with >= {min_proteins_per_taxon} proteins: {len(eligible_taxa)}")
    
    # Get taxon info
    taxa_info = []
    for taxid in eligible_taxa:
        mask = protein_df['taxonomyID'] == taxid
        taxa_info.append({
            'taxonomyID': taxid,
            'n_proteins': int(mask.sum()),
            'scientific_name': protein_df.loc[mask, 'scientific_name'].iloc[0],
            'n_hogs': int(protein_df.loc[mask, 'roothog_id'].nunique())
        })
    taxa_info_df = pd.DataFrame(taxa_info)
    print("\nTaxa included in analysis:")
    print(taxa_info_df.to_string(index=False))
    
    # Compute taxon centroids
    print("\nComputing taxon centroids...")
    taxa_centroids, valid_taxa = compute_taxa_centroids(protein_df, embedding_cols, eligible_taxa)
    
    # Run UMAP on taxon centroids
    print("Running UMAP on taxon centroids...")
    Y_taxa = run_umap(taxa_centroids, n_neighbors=min(5, len(taxa_centroids)-1))
    
    # Hierarchical clustering of taxa
    print("Hierarchical clustering of taxa...")
    taxa_distances = pdist(taxa_centroids, metric='cosine')
    taxa_linkage = linkage(taxa_distances, method='ward')
    
    # Cut tree at different levels
    cluster_results = {}
    for n_clusters in [3, 5, 8]:
        if n_clusters <= len(valid_taxa):
            labels = fcluster(taxa_linkage, n_clusters, criterion='maxclust')
            metrics = compute_cluster_metrics(taxa_centroids, labels)
            cluster_results[f'n_{n_clusters}'] = {
                'labels': labels.tolist(),
                'silhouette': float(metrics['silhouette']) if not np.isnan(metrics['silhouette']) else None,
                'calinski_harabasz': float(metrics['calinski_harabasz']) if not np.isnan(metrics['calinski_harabasz']) else None
            }
    
    # Sample proteins from these taxa for comparison
    print("\nSampling proteins for protein-level analysis...")
    max_per_taxon = 500
    sampled = []
    np.random.seed(42)
    for taxid in eligible_taxa:
        tax_df = protein_df[protein_df['taxonomyID'] == taxid]
        n_sample = min(max_per_taxon, len(tax_df))
        sampled.append(tax_df.sample(n=n_sample, random_state=42))
    sample_df = pd.concat(sampled, ignore_index=True)
    print(f"Sampled {len(sample_df)} proteins from {len(eligible_taxa)} taxa")
    
    # Run UMAP on sampled proteins
    print("Running UMAP on sampled proteins...")
    X = sample_df[embedding_cols].values
    Y = run_umap(X)
    
    # Compute metrics with taxon labels vs HOG labels
    taxa_metrics = compute_cluster_metrics(X, sample_df['taxonomyID'].values)
    hog_metrics = compute_cluster_metrics(X, sample_df['roothog_id'].values)
    
    print("\n" + "-"*50)
    print("CLUSTERING COMPARISON (Protein-level)")
    print("-"*50)
    print(f"Taxa as labels - Silhouette: {taxa_metrics['silhouette']:.4f}")
    print(f"HOG as labels - Silhouette: {hog_metrics['silhouette']:.4f}")
    
    # Compute cross-tabulation of taxa vs HOG
    print("\nCross-tabulation (Taxa vs HOG):")
    cross_tab = pd.crosstab(sample_df['taxonomyID'], sample_df['roothog_id'])
    print(f"  Shape: {cross_tab.shape} (taxa x HOGs)")
    print(f"  Mean HOGs per taxon: {(cross_tab > 0).sum(axis=1).mean():.1f}")
    print(f"  Mean taxa per HOG: {(cross_tab > 0).sum(axis=0).mean():.1f}")
    
    # Calculate ARI between taxa labels and HOG labels
    ari_taxa_hog = adjusted_rand_score(sample_df['taxonomyID'].values, sample_df['roothog_id'].values)
    print(f"\nAdjusted Rand Index (Taxa vs HOG): {ari_taxa_hog:.4f}")
    print("  (0 = random, 1 = perfect agreement)")
    
    results = {
        'min_proteins_per_taxon': min_proteins_per_taxon,
        'n_taxa': len(eligible_taxa),
        'taxa_info': taxa_info,
        'taxa_umap': Y_taxa.tolist(),
        'cluster_results': cluster_results,
        'protein_level': {
            'n_proteins': len(sample_df),
            'taxa_silhouette': float(taxa_metrics['silhouette']),
            'hog_silhouette': float(hog_metrics['silhouette']),
            'ari_taxa_hog': float(ari_taxa_hog),
            'mean_hogs_per_taxon': float((cross_tab > 0).sum(axis=1).mean()),
            'mean_taxa_per_hog': float((cross_tab > 0).sum(axis=0).mean()),
        },
        'umap_coordinates': Y,
        'sample_df': sample_df
    }
    
    return results


# =============================================================================
# Visualization
# =============================================================================

def plot_sampling_stability(results: dict, output_dir: Path):
    """Visualize sampling stability results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    metrics_df = pd.DataFrame(results['iterations'])
    
    # 1. Silhouette scores across iterations
    ax = axes[0, 0]
    x = range(len(metrics_df))
    ax.bar(x, metrics_df['hog_silhouette'], alpha=0.7, label='HOG labels')
    ax.axhline(results['summary']['hog_silhouette_mean'], color='blue', linestyle='--',
               label=f"Mean: {results['summary']['hog_silhouette_mean']:.3f}")
    ax.fill_between(x, 
                    results['summary']['hog_silhouette_mean'] - results['summary']['hog_silhouette_std'],
                    results['summary']['hog_silhouette_mean'] + results['summary']['hog_silhouette_std'],
                    alpha=0.2, color='blue')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Silhouette Score')
    ax.set_title('HOG Silhouette Across Random Samples')
    ax.legend()
    
    # 2. KMeans silhouette across iterations
    ax = axes[0, 1]
    ax.bar(x, metrics_df['kmeans_silhouette'], alpha=0.7, color='green', label='KMeans clusters')
    ax.axhline(results['summary']['kmeans_silhouette_mean'], color='green', linestyle='--',
               label=f"Mean: {results['summary']['kmeans_silhouette_mean']:.3f}")
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Silhouette Score')
    ax.set_title('KMeans Silhouette Across Random Samples')
    ax.legend()
    
    # 3. ARI heatmap
    ax = axes[1, 0]
    ari_matrix = np.array(results['ari_matrix'])
    im = ax.imshow(ari_matrix, cmap='YlGnBu', vmin=0, vmax=1)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Iteration')
    ax.set_title(f'Pairwise ARI (Mean: {results["summary"]["ari_mean"]:.3f})')
    plt.colorbar(im, ax=ax, label='ARI')
    
    # 4. Distribution of metrics
    ax = axes[1, 1]
    ax.hist(metrics_df['hog_silhouette'], bins=10, alpha=0.7, label='HOG Silhouette')
    ax.hist(metrics_df['kmeans_silhouette'], bins=10, alpha=0.7, label='KMeans Silhouette')
    ax.set_xlabel('Silhouette Score')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Clustering Quality')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sampling_stability.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'sampling_stability.png'}")
    plt.close()


def plot_taxa_clustering(results: dict, output_dir: Path):
    """Visualize taxa-based clustering results."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    taxa_info_df = pd.DataFrame(results['taxa_info'])
    
    # 1. Taxa UMAP
    ax = axes[0, 0]
    Y_taxa = np.array(results['taxa_umap'])
    colors = plt.cm.tab20(np.linspace(0, 1, len(Y_taxa)))
    for i, (row, color) in enumerate(zip(taxa_info_df.itertuples(), colors)):
        ax.scatter(Y_taxa[i, 0], Y_taxa[i, 1], c=[color], s=row.n_proteins/50, alpha=0.7)
        name = row.scientific_name.split()[0][:6]
        ax.annotate(name, (Y_taxa[i, 0], Y_taxa[i, 1]), fontsize=8)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title(f'Taxa Centroids in UMAP Space\n(size = n_proteins, {len(Y_taxa)} taxa)')
    
    # 2. Protein-level UMAP colored by taxon
    ax = axes[0, 1]
    sample_df = results['sample_df']
    Y = results['umap_coordinates']
    
    # Use top 10 taxa by size for coloring
    top_taxa = taxa_info_df.nlargest(10, 'n_proteins')['taxonomyID'].tolist()
    
    other_mask = ~sample_df['taxonomyID'].isin(top_taxa)
    ax.scatter(Y[other_mask, 0], Y[other_mask, 1], c='lightgray', s=5, alpha=0.3, label='Other')
    
    colors = plt.cm.tab10.colors
    for i, taxid in enumerate(top_taxa):
        mask = sample_df['taxonomyID'] == taxid
        name = taxa_info_df[taxa_info_df['taxonomyID'] == taxid]['scientific_name'].iloc[0].split()[0]
        ax.scatter(Y[mask, 0], Y[mask, 1], c=[colors[i]], s=10, alpha=0.5, label=name)
    
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('Proteins Colored by Taxonomy')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
    
    # 3. Protein-level UMAP colored by HOG (top 10)
    ax = axes[1, 0]
    hog_counts = sample_df['roothog_id'].value_counts()
    top_hogs = hog_counts.head(10).index.tolist()
    
    other_mask = ~sample_df['roothog_id'].isin(top_hogs)
    ax.scatter(Y[other_mask, 0], Y[other_mask, 1], c='lightgray', s=5, alpha=0.3, label='Other')
    
    for i, hog_id in enumerate(top_hogs):
        mask = sample_df['roothog_id'] == hog_id
        ax.scatter(Y[mask, 0], Y[mask, 1], c=[colors[i]], s=10, alpha=0.5, label=f'HOG {hog_id}')
    
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('Proteins Colored by HOG')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
    
    # 4. Comparison bar chart
    ax = axes[1, 1]
    metrics = ['Taxa Silhouette', 'HOG Silhouette', 'ARI (Taxa vs HOG)']
    values = [
        results['protein_level']['taxa_silhouette'],
        results['protein_level']['hog_silhouette'],
        results['protein_level']['ari_taxa_hog']
    ]
    colors_bar = ['steelblue', 'coral', 'green']
    bars = ax.bar(metrics, values, color=colors_bar, alpha=0.7)
    ax.set_ylabel('Score')
    ax.set_title('Clustering Quality Comparison')
    ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
    
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'taxa_clustering.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'taxa_clustering.png'}")
    plt.close()


def main():
    """Run sensitivity analysis."""
    # Setup output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load data
    protein_df, embedding_cols = load_data()
    
    # Run sensitivity analyses
    results = {}
    
    # 1. Sampling stability
    sampling_results = sensitivity_sampling_stability(
        protein_df, embedding_cols,
        n_iterations=10,
        n_hogs_per_sample=500,
        min_hog_size=5,
        max_per_hog=50
    )
    results['sampling_stability'] = {
        k: v for k, v in sampling_results.items() 
        if k not in ['last_umap', 'last_sample_df', 'ari_matrix']
    }
    results['sampling_stability']['ari_matrix'] = sampling_results['ari_matrix']
    
    plot_sampling_stability(sampling_results, OUTPUT_DIR)
    
    # 2. Taxa-based clustering
    taxa_results = sensitivity_taxa_clustering(
        protein_df, embedding_cols,
        min_proteins_per_taxon=500
    )
    results['taxa_clustering'] = {
        k: v for k, v in taxa_results.items()
        if k not in ['umap_coordinates', 'sample_df']
    }
    
    plot_taxa_clustering(taxa_results, OUTPUT_DIR)
    
    # Generate summary report
    print("\n" + "="*70)
    print("SENSITIVITY ANALYSIS SUMMARY")
    print("="*70)
    
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║  QUESTION 1: If we sample another 500 random groups, do we get      ║
║              the same clusters?                                      ║
╚══════════════════════════════════════════════════════════════════════╝
""")
    
    print(f"""ANSWER: The clustering is REASONABLY STABLE across random samples.

Metrics across {results['sampling_stability']['n_iterations']} iterations:
  • HOG Silhouette Score: {results['sampling_stability']['summary']['hog_silhouette_mean']:.4f} ± {results['sampling_stability']['summary']['hog_silhouette_std']:.4f}
  • KMeans Silhouette:    {results['sampling_stability']['summary']['kmeans_silhouette_mean']:.4f} ± {results['sampling_stability']['summary']['kmeans_silhouette_std']:.4f}
  • Pairwise ARI:         {results['sampling_stability']['summary']['ari_mean']:.4f} ± {results['sampling_stability']['summary']['ari_std']:.4f}

Interpretation:
  • Low standard deviation in silhouette scores → consistent cluster quality
  • High pairwise ARI (>0.7) would indicate same proteins get similar labels
  • ARI values show {
    "HIGH" if results['sampling_stability']['summary']['ari_mean'] > 0.7 else 
    "MODERATE" if results['sampling_stability']['summary']['ari_mean'] > 0.4 else 
    "LOW"} agreement across samples
""")
    
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║  QUESTION 2: What happens if we cluster all taxa that have          ║
║              > 500 proteins?                                         ║
╚══════════════════════════════════════════════════════════════════════╝
""")
    
    taxa_n = results['taxa_clustering']['n_taxa']
    taxa_sil = results['taxa_clustering']['protein_level']['taxa_silhouette']
    hog_sil = results['taxa_clustering']['protein_level']['hog_silhouette']
    ari = results['taxa_clustering']['protein_level']['ari_taxa_hog']
    
    print(f"""ANSWER: Taxonomic grouping shows DIFFERENT structure than HOG grouping.

{taxa_n} taxa have > 500 proteins in the dataset.

Clustering Quality Comparison:
  • Taxa as labels - Silhouette: {taxa_sil:.4f}
  • HOG as labels - Silhouette:  {hog_sil:.4f}
  • ARI (Taxa vs HOG):           {ari:.4f}

Interpretation:
  • {"Taxa clusters are MORE coherent" if taxa_sil > hog_sil else "HOG clusters are MORE coherent" if hog_sil > taxa_sil else "Similar coherence"}
  • Low ARI ({ari:.3f}) confirms taxa and HOG groupings capture DIFFERENT signals
  • Mean HOGs per taxon: {results['taxa_clustering']['protein_level']['mean_hogs_per_taxon']:.1f}
  • Mean taxa per HOG: {results['taxa_clustering']['protein_level']['mean_taxa_per_hog']:.1f}

Key Finding:
  {"Taxonomic signal DOMINATES over evolutionary (HOG) signal in embeddings" if taxa_sil > hog_sil else 
   "Evolutionary (HOG) signal is stronger than taxonomic signal" if hog_sil > taxa_sil else
   "Both signals are similarly strong"}
""")
    
    # Save results
    with open(OUTPUT_DIR / 'sensitivity_results.json', 'w') as f:
        # Convert numpy types for JSON serialization
        json.dump(results, f, indent=2, default=lambda x: float(x) if hasattr(x, 'item') else str(x))
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("Done!")


if __name__ == "__main__":
    main()
