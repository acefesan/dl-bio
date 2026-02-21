#!/usr/bin/env python3
"""UMAP clustering of ESM2 protein embeddings.

Loads the CAFA3 merged dataset and performs UMAP dimensionality reduction
on protein embeddings to visualize clustering patterns by taxa and HOG.

Usage:
    python umap_embeddings.py                      # Default settings
    python umap_embeddings.py --sample-size 1000   # More samples per taxon
    python umap_embeddings.py --min-taxa 1000      # Only species with >= 1000 proteins
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from umap import UMAP

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "assets/proteins/datasets/cafa3_merged/cafa3_with_embeddings.feather"
OUTPUT_DIR = BASE_PATH / "assets/proteins/analysis"


def load_dataset(data_path: Path) -> pd.DataFrame:
    """Load the CAFA3 merged dataset."""
    print(f"Loading dataset from {data_path.name}...")
    df = pd.read_feather(data_path)
    print(f"  Loaded {len(df):,} rows, {df['EntryID'].nunique():,} unique proteins")
    return df


def get_unique_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Get unique proteins with embeddings (deduplicated from GO term expansion)."""
    embedding_cols = [c for c in df.columns if c.startswith('ME:')]
    protein_cols = ['EntryID', 'Length', 'taxonomyID', 'scientific_name',
                    'oma_id', 'hog_id', 'roothog_id'] + embedding_cols

    # Deduplicate and filter to proteins with embeddings
    protein_df = df[protein_cols].drop_duplicates(subset='EntryID')
    protein_df = protein_df[protein_df[embedding_cols[0]].notna()]

    print(f"  Unique proteins with embeddings: {len(protein_df):,}")
    return protein_df


def get_abundant_taxa(protein_df: pd.DataFrame, min_count: int = 500) -> dict:
    """Get taxa with at least min_count proteins."""
    taxa_counts = protein_df.groupby('taxonomyID').size()
    taxa_counts = taxa_counts.sort_values(ascending=False)

    abundant = taxa_counts[taxa_counts >= min_count]

    # Get scientific names for each taxon
    taxid_to_name = protein_df.groupby('taxonomyID')['scientific_name'].first().to_dict()

    # Create mapping: taxonomyID -> short name
    taxa_map = {}
    for taxid, count in abundant.items():
        full_name = taxid_to_name.get(taxid)
        if full_name and isinstance(full_name, str):
            parts = full_name.split()
            if len(parts) >= 2:
                short_name = f"{parts[0][0]}. {parts[1]}"
            else:
                short_name = full_name[:15]
        else:
            short_name = f"Tax_{taxid}"
        taxa_map[taxid] = short_name

    print(f"\nAbundant taxa (>= {min_count} proteins): {len(taxa_map)}")
    for taxid, name in list(taxa_map.items())[:10]:
        count = abundant.loc[taxid]
        print(f"  {name}: {count:,}")
    if len(taxa_map) > 10:
        print(f"  ... and {len(taxa_map) - 10} more")

    return taxa_map


def get_abundant_hogs(protein_df: pd.DataFrame, min_count: int = 100) -> pd.DataFrame:
    """Get root HOGs with at least min_count proteins."""
    hog_counts = protein_df['roothog_id'].value_counts()
    abundant = hog_counts[hog_counts >= min_count]

    # Exclude HOG 0 (no assignment)
    abundant = abundant[abundant.index != 0]

    print(f"\nAbundant root HOGs (>= {min_count} proteins): {len(abundant)}")
    print(abundant.head(10))

    return abundant


def sample_proteins(
    protein_df: pd.DataFrame,
    taxa_map: dict,
    sample_size: int = 500,
    seed: int = 42
) -> pd.DataFrame:
    """Sample proteins from each taxon."""
    np.random.seed(seed)

    sampled_dfs = []
    for taxid, short_name in taxa_map.items():
        taxon_df = protein_df[protein_df['taxonomyID'] == taxid]
        n_sample = min(sample_size, len(taxon_df))
        sampled = taxon_df.sample(n=n_sample, random_state=seed)
        sampled_dfs.append(sampled)

    sample_df = pd.concat(sampled_dfs, ignore_index=True)
    sample_df['taxon_name'] = sample_df['taxonomyID'].map(taxa_map)

    print(f"\nSampled {len(sample_df):,} proteins from {len(taxa_map)} taxa")
    return sample_df


def run_umap(X: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    """Run UMAP dimensionality reduction."""
    print("\nRunning UMAP...")
    umap_model = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        metric='cosine',
        random_state=42,
        verbose=True
    )
    return umap_model.fit_transform(X)


def plot_by_taxa(sample_df: pd.DataFrame, output_path: Path):
    """Create UMAP plot colored by taxa."""
    fig, ax = plt.subplots(figsize=(14, 10))

    # Sort taxa by count for legend
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
    """Create UMAP plot colored by top root HOGs."""
    fig, ax = plt.subplots(figsize=(14, 10))

    top_hog_ids = abundant_hogs.index.tolist()[:15]  # Top 15 HOGs

    # Background: other/no HOG
    other_mask = ~sample_df['roothog_id'].isin(top_hog_ids)
    other_df = sample_df[other_mask]
    ax.scatter(
        other_df['umap_x'], other_df['umap_y'],
        label=f'Other/No HOG (n={len(other_df)})',
        alpha=0.2, s=10, c='lightgray'
    )

    # Top HOGs
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
    """Create grid showing HOG distribution across taxa."""
    # Get top taxa and HOGs
    top_taxa = sample_df['taxon_name'].value_counts().head(12).index.tolist()
    top_hogs = abundant_hogs.head(12).index.tolist()

    # Create cross-tabulation
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


def save_hog_frequencies(protein_df: pd.DataFrame, output_path: Path):
    """Save HOG frequency table."""
    # Filter out NaN roothog_id and 0 (no assignment)
    valid_df = protein_df[protein_df['roothog_id'].notna() & (protein_df['roothog_id'] != 0)]

    hog_counts = valid_df.groupby('roothog_id').agg(
        n_proteins=('EntryID', 'nunique'),
        n_taxa=('taxonomyID', 'nunique'),
        taxa_list=('taxonomyID', lambda x: ', '.join(map(str, sorted(x.unique())[:5])))
    ).sort_values('n_proteins', ascending=False)

    hog_counts.index = hog_counts.index.astype(int)
    hog_counts.index.name = 'roothog_id'

    hog_counts.to_csv(output_path)
    print(f"Saved HOG frequencies: {output_path}")
    print(f"  Total root HOGs: {len(hog_counts):,}")
    print(hog_counts.head(15))


def main():
    parser = argparse.ArgumentParser(description="UMAP clustering of protein embeddings")
    parser.add_argument('--sample-size', type=int, default=500,
                        help='Proteins to sample per taxon (default: 500)')
    parser.add_argument('--min-taxa', type=int, default=500,
                        help='Minimum proteins for taxon inclusion (default: 500)')
    parser.add_argument('--min-hog', type=int, default=100,
                        help='Minimum proteins for HOG inclusion (default: 100)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for plots')
    args = parser.parse_args()

    # Setup output directory
    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load and prepare data
    df = load_dataset(DATA_PATH)
    protein_df = get_unique_proteins(df)

    # Get abundant taxa and HOGs
    taxa_map = get_abundant_taxa(protein_df, min_count=args.min_taxa)
    abundant_hogs = get_abundant_hogs(protein_df, min_count=args.min_hog)

    # Save HOG frequencies
    save_hog_frequencies(protein_df, output_dir / 'hog_frequencies.csv')

    # Sample proteins
    sample_df = sample_proteins(protein_df, taxa_map, sample_size=args.sample_size)

    # Extract embeddings and run UMAP
    embedding_cols = [c for c in sample_df.columns if c.startswith('ME:')]
    X = sample_df[embedding_cols].values
    print(f"\nEmbedding matrix: {X.shape}")

    umap_coords = run_umap(X)
    sample_df['umap_x'] = umap_coords[:, 0]
    sample_df['umap_y'] = umap_coords[:, 1]

    # Generate plots
    print("\nGenerating plots...")
    plot_by_taxa(sample_df, output_dir / 'umap_by_taxa.png')
    plot_by_hog(sample_df, abundant_hogs, output_dir / 'umap_by_hog.png')
    plot_taxa_hog_grid(sample_df, abundant_hogs, output_dir / 'taxa_hog_heatmap.png')

    # Save UMAP coordinates
    coords_df = sample_df[['EntryID', 'taxonomyID', 'taxon_name', 'roothog_id', 'umap_x', 'umap_y']]
    coords_df.to_csv(output_dir / 'umap_coordinates.csv', index=False)
    print(f"Saved UMAP coordinates: {output_dir / 'umap_coordinates.csv'}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Taxa analyzed: {len(taxa_map)}")
    print(f"Proteins sampled: {len(sample_df):,}")
    print(f"Abundant HOGs (>= {args.min_hog}): {len(abundant_hogs)}")
    print(f"Output directory: {output_dir}")
    print("\nDone!")


if __name__ == "__main__":
    main()
