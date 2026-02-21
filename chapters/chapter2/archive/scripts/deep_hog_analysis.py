#!/usr/bin/env python3
"""Deep HOG Hierarchy Analysis - Addressing Issues in Original Analysis

This script addresses the methodological issues identified in the critical review:

1. Analyzes sub-HOG hierarchies (not just root HOGs)
2. Computes metrics in ORIGINAL embedding space (not UMAP projections)
3. Tests whether deeper hierarchy levels cluster more tightly
4. Correlates embedding distance with sequence identity (when available)

Key insight: If embeddings capture evolutionary relationships, proteins at 
deeper hierarchy levels (more recent divergence) should cluster MORE tightly.

Usage:
    python deep_hog_analysis.py
"""

import json
import warnings
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import KMeans
from Levenshtein import ratio as levenshtein_ratio

warnings.filterwarnings('ignore')

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
HOG_CACHE_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/hog_cache.csv"
OUTPUT_DIR = BASE_PATH / "assets/proteins/analysis/deep_hog_analysis"


def load_data():
    """Load dataset and HOG cache."""
    print("Loading dataset...")
    df = pd.read_feather(DATA_PATH)
    
    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    keep_cols = ['EntryID', 'taxonomyID', 'scientific_name', 'Length',
                 'roothog_id', 'hog_id', 'Sequence'] + embedding_cols
    
    protein_df = df[keep_cols].drop_duplicates(subset='EntryID')
    protein_df = protein_df[protein_df[embedding_cols[0]].notna()]
    protein_df = protein_df[protein_df['roothog_id'].notna()]
    protein_df = protein_df[protein_df['roothog_id'] != 0]
    
    # Load HOG cache with hierarchy info
    hog_cache = pd.read_csv(HOG_CACHE_PATH)
    
    print(f"  Proteins with embeddings: {len(protein_df):,}")
    print(f"  HOG cache entries: {len(hog_cache):,}")
    
    return protein_df, embedding_cols, hog_cache


def get_hog_depth(hog_id):
    """Get hierarchy depth from HOG ID string."""
    if pd.isna(hog_id) or hog_id == '' or hog_id == '0':
        return 0
    return str(hog_id).count('.')


def compute_silhouette_in_original_space(embeddings, labels, metric='cosine', sample_size=None):
    """
    Compute silhouette score in ORIGINAL embedding space.
    
    This is the correct approach - NOT computing in UMAP space.
    """
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    if n_clusters < 2:
        return np.nan, np.array([])
    
    # Sample if too large (silhouette is O(n²))
    if sample_size and len(embeddings) > sample_size:
        idx = np.random.choice(len(embeddings), sample_size, replace=False)
        embeddings = embeddings[idx]
        labels = labels[idx]
    
    # Compute in original high-dimensional space
    score = silhouette_score(embeddings, labels, metric=metric)
    sample_scores = silhouette_samples(embeddings, labels, metric=metric)
    
    return score, sample_scores


def analyze_hierarchy_depth_clustering(protein_df, embedding_cols, hog_cache, roothog_id):
    """
    Analyze how clustering quality varies with HOG hierarchy depth.
    
    Key hypothesis: Deeper sub-HOGs (more recent divergence) should cluster
    more tightly than shallow sub-HOGs (ancient divergence).
    """
    print(f"\n{'='*70}")
    print(f"ANALYZING ROOT HOG {roothog_id}")
    print(f"{'='*70}")
    
    # Get all proteins in this root HOG
    protein_mask = protein_df['roothog_id'] == roothog_id
    hog_proteins = protein_df[protein_mask].copy()
    
    # Merge with HOG cache to get full hog_id
    hog_info = hog_cache[hog_cache['roothog_id'] == roothog_id][['EntryID', 'hog_id']]
    hog_proteins = hog_proteins.merge(hog_info, on='EntryID', how='left', suffixes=('', '_cache'))
    
    # Compute depth for each protein
    hog_proteins['depth'] = hog_proteins['hog_id_cache'].apply(get_hog_depth)
    
    print(f"  Total proteins: {len(hog_proteins)}")
    print(f"  Depth distribution:")
    depth_counts = hog_proteins['depth'].value_counts().sort_index()
    for depth, count in depth_counts.items():
        print(f"    Level {depth}: {count} proteins")
    
    # Get embeddings
    X = hog_proteins[embedding_cols].values
    
    # For each depth level, compute intra-group metrics
    results = []
    
    for depth in sorted(hog_proteins['depth'].unique()):
        if depth == 0:
            continue
            
        depth_mask = hog_proteins['depth'] == depth
        n_proteins = depth_mask.sum()
        
        if n_proteins < 5:
            continue
        
        X_depth = X[depth_mask]
        
        # Compute centroid
        centroid = X_depth.mean(axis=0)
        
        # Compute intra-group distances (in ORIGINAL space)
        intra_distances = cdist([centroid], X_depth, metric='cosine')[0]
        
        # Compute pairwise distances for a sample
        if len(X_depth) > 100:
            idx = np.random.choice(len(X_depth), 100, replace=False)
            X_sample = X_depth[idx]
        else:
            X_sample = X_depth
        
        pairwise = pdist(X_sample, metric='cosine')
        
        results.append({
            'depth': depth,
            'n_proteins': n_proteins,
            'mean_dist_to_centroid': float(intra_distances.mean()),
            'std_dist_to_centroid': float(intra_distances.std()),
            'median_pairwise_dist': float(np.median(pairwise)),
            'max_pairwise_dist': float(np.max(pairwise)) if len(pairwise) > 0 else np.nan,
        })
    
    results_df = pd.DataFrame(results)
    
    # Compute correlation: Does deeper hierarchy = tighter clustering?
    if len(results_df) >= 3:
        corr_centroid, p_centroid = spearmanr(results_df['depth'], 
                                               results_df['mean_dist_to_centroid'])
        corr_pairwise, p_pairwise = spearmanr(results_df['depth'], 
                                               results_df['median_pairwise_dist'])
    else:
        corr_centroid, p_centroid = np.nan, np.nan
        corr_pairwise, p_pairwise = np.nan, np.nan
    
    print(f"\n  Clustering by hierarchy depth:")
    print(results_df.to_string(index=False))
    
    print(f"\n  Correlation (depth vs tightness):")
    print(f"    Depth vs Mean-Centroid-Dist: r = {corr_centroid:.3f} (p = {p_centroid:.3f})")
    print(f"    Depth vs Median-Pairwise:    r = {corr_pairwise:.3f} (p = {p_pairwise:.3f})")
    
    # Interpretation
    if corr_centroid < -0.3 and p_centroid < 0.05:
        interpretation = "✅ GOOD: Deeper sub-HOGs cluster more tightly (as expected)"
    elif corr_centroid > 0.3 and p_centroid < 0.05:
        interpretation = "❌ BAD: Deeper sub-HOGs are MORE dispersed (unexpected)"
    else:
        interpretation = "⚪ NEUTRAL: No significant relationship between depth and clustering"
    
    print(f"\n  Interpretation: {interpretation}")
    
    return {
        'roothog_id': roothog_id,
        'n_proteins': len(hog_proteins),
        'max_depth': int(hog_proteins['depth'].max()),
        'depth_results': results,
        'corr_centroid': float(corr_centroid) if not np.isnan(corr_centroid) else None,
        'p_centroid': float(p_centroid) if not np.isnan(p_centroid) else None,
        'corr_pairwise': float(corr_pairwise) if not np.isnan(corr_pairwise) else None,
        'p_pairwise': float(p_pairwise) if not np.isnan(p_pairwise) else None,
        'interpretation': interpretation,
    }


def compute_sequence_embedding_correlation(protein_df, embedding_cols, n_pairs=5000, seed=42):
    """
    THE KEY ANALYSIS: Does embedding distance correlate with sequence divergence?
    
    This directly tests whether ESM2 embeddings capture evolutionary relationships.
    """
    print("\n" + "="*70)
    print("SEQUENCE IDENTITY vs EMBEDDING DISTANCE CORRELATION")
    print("="*70)
    print("This is the most important metric for evaluating if embeddings")
    print("capture evolutionary relationships.\n")
    
    np.random.seed(seed)
    
    # Sample protein pairs from the same root HOG
    hog_counts = protein_df.groupby('roothog_id').size()
    large_hogs = hog_counts[hog_counts >= 10].index.tolist()
    
    embedding_distances = []
    sequence_similarities = []
    
    print(f"Sampling {n_pairs} protein pairs from {len(large_hogs)} large HOGs...")
    
    pairs_collected = 0
    attempts = 0
    max_attempts = n_pairs * 20
    
    while pairs_collected < n_pairs and attempts < max_attempts:
        attempts += 1
        
        # Pick a random HOG
        hog = np.random.choice(large_hogs)
        hog_proteins = protein_df[protein_df['roothog_id'] == hog]
        
        if len(hog_proteins) < 2:
            continue
        
        # Pick two random proteins
        idx = hog_proteins.sample(n=2, random_state=seed + attempts)
        
        seq1 = idx.iloc[0]['Sequence']
        seq2 = idx.iloc[1]['Sequence']
        
        if pd.isna(seq1) or pd.isna(seq2):
            continue
        
        # Compute sequence similarity (normalized Levenshtein)
        seq_sim = levenshtein_ratio(seq1, seq2)
        
        # Compute embedding distance (cosine)
        emb1 = idx.iloc[0][embedding_cols].values
        emb2 = idx.iloc[1][embedding_cols].values
        emb_dist = cdist([emb1], [emb2], metric='cosine')[0, 0]
        
        embedding_distances.append(emb_dist)
        sequence_similarities.append(seq_sim)
        pairs_collected += 1
        
        if pairs_collected % 1000 == 0:
            print(f"  Collected {pairs_collected} pairs...")
    
    embedding_distances = np.array(embedding_distances)
    sequence_similarities = np.array(sequence_similarities)
    
    # Compute correlations
    # Note: High sequence similarity should correlate with LOW embedding distance
    spearman_r, spearman_p = spearmanr(sequence_similarities, embedding_distances)
    pearson_r, pearson_p = pearsonr(sequence_similarities, embedding_distances)
    
    print(f"\nResults ({pairs_collected} pairs analyzed):")
    print(f"  Spearman correlation: r = {spearman_r:.4f} (p = {spearman_p:.2e})")
    print(f"  Pearson correlation:  r = {pearson_r:.4f} (p = {pearson_p:.2e})")
    
    print(f"\n  Sequence similarity range: [{sequence_similarities.min():.3f}, {sequence_similarities.max():.3f}]")
    print(f"  Embedding distance range:  [{embedding_distances.min():.3f}, {embedding_distances.max():.3f}]")
    
    # Interpretation
    # Negative correlation means: higher sequence similarity → lower embedding distance
    # This is what we WANT
    if spearman_r < -0.5:
        interpretation = "✅ STRONG: Embeddings capture sequence similarity well"
    elif spearman_r < -0.3:
        interpretation = "🟡 MODERATE: Embeddings partially capture sequence similarity"
    elif spearman_r < -0.1:
        interpretation = "🟠 WEAK: Embeddings weakly capture sequence similarity"
    else:
        interpretation = "❌ NONE: Embeddings do NOT capture sequence similarity"
    
    print(f"\n  Interpretation: {interpretation}")
    
    return {
        'n_pairs': pairs_collected,
        'spearman_r': float(spearman_r),
        'spearman_p': float(spearman_p),
        'pearson_r': float(pearson_r),
        'pearson_p': float(pearson_p),
        'seq_sim_range': [float(sequence_similarities.min()), float(sequence_similarities.max())],
        'emb_dist_range': [float(embedding_distances.min()), float(embedding_distances.max())],
        'interpretation': interpretation,
        'embedding_distances': embedding_distances.tolist(),
        'sequence_similarities': sequence_similarities.tolist(),
    }


def compare_silhouette_spaces(protein_df, embedding_cols, sample_size=5000, seed=42):
    """
    Compare silhouette scores computed in different spaces:
    1. Original 640D embedding space
    2. UMAP 2D projection
    
    This demonstrates the flaw in the original analysis.
    """
    print("\n" + "="*70)
    print("SILHOUETTE SCORE: ORIGINAL SPACE vs UMAP SPACE")
    print("="*70)
    print("The original analysis computed silhouette in UMAP space,")
    print("which distorts distances. This comparison shows the difference.\n")
    
    from umap import UMAP
    
    np.random.seed(seed)
    
    # Sample proteins
    sample_df = protein_df.sample(n=min(sample_size, len(protein_df)), random_state=seed)
    X = sample_df[embedding_cols].values
    labels = sample_df['roothog_id'].values
    
    # Filter to HOGs with >= 2 proteins
    unique_hogs, counts = np.unique(labels, return_counts=True)
    valid_hogs = unique_hogs[counts >= 2]
    valid_mask = np.isin(labels, valid_hogs)
    X = X[valid_mask]
    labels = labels[valid_mask]
    
    print(f"Sample size: {len(X)} proteins, {len(np.unique(labels))} HOGs")
    
    # Silhouette in original space
    print("\nComputing silhouette in 640D embedding space...")
    sil_original = silhouette_score(X, labels, metric='cosine')
    
    # UMAP projection
    print("Running UMAP projection...")
    umap = UMAP(n_neighbors=15, min_dist=0.1, metric='cosine', random_state=seed)
    X_umap = umap.fit_transform(X)
    
    # Silhouette in UMAP space
    print("Computing silhouette in 2D UMAP space...")
    sil_umap = silhouette_score(X_umap, labels, metric='euclidean')
    
    print(f"\n  Silhouette in 640D space: {sil_original:.4f}")
    print(f"  Silhouette in 2D UMAP:    {sil_umap:.4f}")
    print(f"  Difference:               {sil_umap - sil_original:+.4f}")
    
    if abs(sil_umap - sil_original) > 0.1:
        interpretation = "⚠️ WARNING: UMAP significantly changes silhouette scores"
    else:
        interpretation = "✓ OK: UMAP roughly preserves silhouette scores"
    
    print(f"\n  {interpretation}")
    
    return {
        'sample_size': len(X),
        'n_hogs': len(np.unique(labels)),
        'silhouette_original': float(sil_original),
        'silhouette_umap': float(sil_umap),
        'difference': float(sil_umap - sil_original),
        'interpretation': interpretation,
    }


def plot_results(hierarchy_results, seq_emb_results, output_dir):
    """Generate visualization plots."""
    
    # Plot 1: Hierarchy depth vs clustering tightness
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    for i, result in enumerate(hierarchy_results[:4]):
        if i >= 4:
            break
        ax = axes[i // 2, i % 2]
        
        depth_data = pd.DataFrame(result['depth_results'])
        if len(depth_data) > 0:
            ax.errorbar(depth_data['depth'], depth_data['mean_dist_to_centroid'],
                       yerr=depth_data['std_dist_to_centroid'], 
                       marker='o', capsize=5, alpha=0.7)
            ax.set_xlabel('HOG Hierarchy Depth')
            ax.set_ylabel('Mean Distance to Centroid')
            ax.set_title(f"Root HOG {result['roothog_id']}\n(r={result['corr_centroid']:.2f})")
            
            # Add trend line
            z = np.polyfit(depth_data['depth'], depth_data['mean_dist_to_centroid'], 1)
            p = np.poly1d(z)
            ax.plot(depth_data['depth'], p(depth_data['depth']), "r--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'hierarchy_depth_clustering.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'hierarchy_depth_clustering.png'}")
    plt.close()
    
    # Plot 2: Sequence similarity vs embedding distance
    fig, ax = plt.subplots(figsize=(10, 8))
    
    seq_sim = np.array(seq_emb_results['sequence_similarities'])
    emb_dist = np.array(seq_emb_results['embedding_distances'])
    
    # Hex bin for density
    hb = ax.hexbin(seq_sim, emb_dist, gridsize=30, cmap='YlOrRd', mincnt=1)
    plt.colorbar(hb, ax=ax, label='Count')
    
    ax.set_xlabel('Sequence Similarity (Levenshtein ratio)', fontsize=12)
    ax.set_ylabel('Embedding Distance (Cosine)', fontsize=12)
    ax.set_title(f'Sequence Similarity vs Embedding Distance\n'
                 f'Spearman r = {seq_emb_results["spearman_r"]:.3f}', fontsize=14)
    
    # Add trend line
    z = np.polyfit(seq_sim, emb_dist, 1)
    p = np.poly1d(z)
    x_line = np.linspace(seq_sim.min(), seq_sim.max(), 100)
    ax.plot(x_line, p(x_line), 'b-', linewidth=2, label='Linear fit')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sequence_embedding_correlation.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_dir / 'sequence_embedding_correlation.png'}")
    plt.close()


def main():
    """Run the corrected deep HOG analysis."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load data
    protein_df, embedding_cols, hog_cache = load_data()
    
    # Add depth info to HOG cache
    hog_cache['depth'] = hog_cache['hog_id'].apply(get_hog_depth)
    
    # Find HOGs with deepest hierarchies and sufficient proteins
    roothog_stats = hog_cache.groupby('roothog_id').agg({
        'depth': 'max',
        'EntryID': 'count'
    }).rename(columns={'EntryID': 'n_proteins', 'depth': 'max_depth'})
    
    # Select top HOGs with deep hierarchies
    candidates = roothog_stats[(roothog_stats['max_depth'] >= 8) & 
                               (roothog_stats['n_proteins'] >= 50)]
    candidates = candidates.sort_values('max_depth', ascending=False).head(10)
    
    print("\n" + "="*70)
    print("SELECTED HOGs FOR DEEP ANALYSIS")
    print("="*70)
    print(candidates)
    
    results = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'methodology': 'Corrected analysis with original embedding space metrics',
        'hierarchy_analysis': [],
        'sequence_embedding_correlation': None,
        'silhouette_comparison': None,
    }
    
    # 1. Analyze hierarchy depth clustering for top HOGs
    for roothog_id in candidates.index[:5]:
        hierarchy_result = analyze_hierarchy_depth_clustering(
            protein_df, embedding_cols, hog_cache, int(roothog_id)
        )
        results['hierarchy_analysis'].append(hierarchy_result)
    
    # 2. Compute sequence-embedding correlation
    results['sequence_embedding_correlation'] = compute_sequence_embedding_correlation(
        protein_df, embedding_cols, n_pairs=5000
    )
    
    # 3. Compare silhouette in different spaces
    results['silhouette_comparison'] = compare_silhouette_spaces(
        protein_df, embedding_cols, sample_size=5000
    )
    
    # Generate plots
    print("\n" + "="*70)
    print("GENERATING VISUALIZATIONS")
    print("="*70)
    plot_results(results['hierarchy_analysis'], 
                 results['sequence_embedding_correlation'], 
                 OUTPUT_DIR)
    
    # Save results
    # Convert non-serializable items
    results_json = results.copy()
    with open(OUTPUT_DIR / 'deep_analysis_results.json', 'w') as f:
        json.dump(results_json, f, indent=2, default=str)
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    
    # Final summary
    print("\n" + "="*70)
    print("SUMMARY OF FINDINGS")
    print("="*70)
    
    print("\n1. HIERARCHY DEPTH CLUSTERING:")
    for result in results['hierarchy_analysis']:
        print(f"   HOG {result['roothog_id']}: {result['interpretation']}")
    
    print(f"\n2. SEQUENCE-EMBEDDING CORRELATION:")
    print(f"   {results['sequence_embedding_correlation']['interpretation']}")
    print(f"   Spearman r = {results['sequence_embedding_correlation']['spearman_r']:.4f}")
    
    print(f"\n3. SILHOUETTE SPACE COMPARISON:")
    print(f"   {results['silhouette_comparison']['interpretation']}")
    print(f"   Original (640D): {results['silhouette_comparison']['silhouette_original']:.4f}")
    print(f"   UMAP (2D):       {results['silhouette_comparison']['silhouette_umap']:.4f}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
