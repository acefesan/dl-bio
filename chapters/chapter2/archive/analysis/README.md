# UMAP Embedding Analysis

Analysis of ESM2 protein embeddings from the CAFA3 dataset using UMAP dimensionality reduction.

## Dataset

- **Source**: CAFA3 merged dataset with ESM2 embeddings
- **Proteins**: 141,659 unique proteins with embeddings
- **Embedding model**: `facebook/esm2_t33_650M_UR50D` (1280 dimensions)
- **Taxa analyzed**: 22 species with ≥500 proteins each

## Files

| File | Description |
|------|-------------|
| `umap_by_taxa.png` | UMAP projection colored by species (22 taxa, 500 samples each) |
| `umap_by_hog.png` | UMAP projection colored by top 15 root HOGs |
| `taxa_hog_heatmap.png` | Cross-tabulation of HOG distribution across taxa |
| `umap_stability_seeds.png` | UMAP stability across 5 random seeds |
| `hog_frequencies.csv` | Root HOG frequencies (36,148 HOGs) |
| `umap_coordinates.csv` | UMAP coordinates for sampled proteins |
| `umap_seed_*.csv` | UMAP coordinates for each random seed (42, 123, 456, 789, 1024) |
| `umap_stability_summary.json` | Procrustes disparity metrics |

## UMAP Stability Analysis

UMAP was run with 5 different random seeds to assess cluster stability.

### Procrustes Disparity Results

| Seed Pair | Disparity | Interpretation |
|-----------|-----------|----------------|
| 42 vs 123 | 0.4157 | Most similar |
| 42 vs 1024 | 0.4534 | Similar |
| 456 vs 1024 | 0.4584 | Similar |
| 123 vs 1024 | 0.4761 | Moderate |
| 789 vs 1024 | 0.5716 | Moderate |
| 42 vs 789 | 0.6218 | Different |
| 42 vs 456 | 0.6489 | Different |
| 123 vs 456 | 0.6677 | Different |
| 123 vs 789 | 0.6820 | Different |
| 456 vs 789 | 0.7252 | Most different |

**Interpretation**: Disparity ranges from 0 (identical) to 1 (completely different).
- Values 0.4-0.5: Good stability, clusters well preserved
- Values 0.5-0.7: Moderate stability, global structure preserved but positions shift
- Values >0.7: Poor stability, significant structural changes

### Findings

1. **Moderate instability** (disparity 0.42-0.73) is typical for UMAP with high-dimensional data
2. **Global structure** (main clusters) is preserved across seeds
3. **Local positions** shift significantly between runs
4. Some seed pairs (42/123, 42/1024) produce more similar results

### Recommendations

To improve stability:
- Increase `n_neighbors` (e.g., 30-50) for more global structure
- Apply PCA pre-reduction (e.g., to 50-100 dimensions) before UMAP
- Use ensemble averaging of multiple UMAP runs
- Consider t-SNE for comparison (different stability characteristics)

## HOG Analysis

- **Total root HOGs**: 36,148
- **Abundant HOGs (≥100 proteins)**: 18
- **Top HOG (801468)**: 2,230 proteins across 77 taxa (likely ribosomal or housekeeping)

## Species Analyzed

| Species | Proteins |
|---------|----------|
| H. sapiens | 25,037 |
| A. thaliana | 14,456 |
| M. musculus | 14,310 |
| D. rerio | 12,601 |
| D. melanogaster | 11,876 |
| R. norvegicus | 8,767 |
| S. cerevisiae | 5,467 |
| T. brucei | 5,179 |
| C. elegans | 4,882 |
| S. pombe | 4,631 |
| ... and 12 more | |

## Usage

View images in the markdown server or directly:

```bash
# View in browser
xdg-open umap_by_taxa.png

# Regenerate analysis
python chapters/umap_embeddings.py --sample-size 500 --min-taxa 500
```

## Generated

- **Date**: 2026-02-03
- **Script**: `chapters/umap_embeddings.py`
