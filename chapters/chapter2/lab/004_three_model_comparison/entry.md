# 004 — Three-Model Comparison (150M vs 650M vs 3B)

**Date:** 2026-03-06
**Models:** ESM2 150M, 650M, 3B
**Status:** complete

## Hypothesis
Comparing ESM2 models across 3 scales (150M, 650M, 3B) reveals whether embedding
quality for taxonomic/HOG clustering scales with model size. The 3B anomaly (entry 002)
suggested larger isn't always better — does 150M follow the same trend or confirm
that 650M is a sweet spot?

## Setup
All models use the same CAFA3 dataset (142k proteins) and sampling strategy (22 taxa ×
500 proteins = 11k). Multi-seed UMAP (seeds 42, 123, 456, 789) for 150M and 3B; single
seed for 650M (run predated multi-seed support).

| Model | Embedding dim | batch_size | max_seq_length |
|-------|--------------|------------|----------------|
| 150M (`esm2_t30_150M_UR50D`) | 640 | 8 | 3000 |
| 650M (`esm2_t33_650M_UR50D`) | 1280 | 32 | 5000 |
| 3B (`esm2_t36_3B_UR50D`) | 2560 | 4 | 3000 |

## Results

### Clustering Metrics (cosine distance, original embedding space)

| Model | Seed | KMeans sil | HOG sil |
|-------|------|-----------|---------|
| 150M  | 42   | **0.0498** | -0.0026 |
| 150M  | 123  | **0.0491** | -0.0004 |
| 150M  | 456  | **0.0509** | 0.0111  |
| 150M  | 789  | **0.0519** | 0.0010  |
| 650M  | 42   | 0.0396    | **0.0313** |
| 3B    | 42   | 0.0003    | -0.2131 |
| 3B    | 123  | 0.0024    | -0.2153 |
| 3B    | 456  | 0.0016    | -0.1988 |
| 3B    | 789  | -0.0008   | -0.2203 |

### Average metrics by model

| Model | Avg KMeans sil | Avg HOG sil |
|-------|---------------|-------------|
| 150M  | **0.0504**    | 0.0023      |
| 650M  | 0.0396        | **0.0313**  |
| 3B    | 0.0009        | -0.2116     |

### Key findings

1. **KMeans clustering: 150M > 650M >> 3B.** The smallest model produces the most
   KMeans-separable clusters. This is consistent across all 4 seeds for 150M
   (sil 0.049–0.052) and stable for 3B (sil ≈ 0).

2. **HOG clustering: 650M > 150M >> 3B.** The 650M model best captures evolutionary
   group structure. 150M is near-zero (proteins in the same HOG are neither clustered
   nor anti-clustered). 3B is strongly negative (HOG members are farther apart than
   random — anti-clustered).

3. **Inverse scaling for clustering quality.** Bigger ESM2 models produce embeddings
   that are *worse* for unsupervised clustering by taxonomy and HOG. This doesn't mean
   the larger models are worse at biology — they may capture finer-grained sequence
   features that orthogonalize coarse taxonomic signal.

## Figures

### UMAP colored by taxa (seed 42 for each model)
- `figures/150m/umap_by_taxa.png`
- `figures/650m/umap_by_taxa.png`
- `figures/3b/umap_by_taxa.png`

### UMAP colored by root HOG
- `figures/150m/umap_by_hog.png`
- `figures/650m/umap_by_hog.png`
- `figures/3b/umap_by_hog.png`

### Subtree HOG overview (top 4 root HOGs highlighted)
- `figures/150m/subtree_overview.png`
- `figures/650m/subtree_overview.png`
- `figures/3b/subtree_overview.png`

### Root HOG 801468 — level-1 sub-HOG coloring
- `figures/150m/subtree_hog_801468_level1.png`
- `figures/650m/subtree_hog_801468_level1.png`
- `figures/3b/subtree_hog_801468_level1.png`

## Interpretation

The inverse scaling pattern is the central finding. Three possible explanations:

1. **Embedding anisotropy**: larger models may concentrate embeddings in a narrow cone,
   making cosine distances less discriminative. The negative 3B HOG silhouette is a
   strong signal that something is geometrically wrong with the embedding space for
   distance-based methods.

2. **Feature granularity**: 150M/650M may default to coarser sequence features (family-
   level motifs), while 3B captures residue-level details that don't correlate with
   taxonomy. The 3B model has more capacity to learn fine distinctions between sequences
   within the same HOG, pushing them apart in embedding space.

3. **Truncation confound**: 650M used max_seq_length=5000 while 150M and 3B used 3000.
   This changes the protein set slightly (253 sequences removed at 5000, more at 3000).
   However, 150M and 3B share the same truncation and show opposite clustering quality,
   so truncation alone doesn't explain the trend.

## Next steps
- **Anisotropy check**: compute singular value spectrum and cosine distance distribution
  for each model's embeddings — test the anisotropy hypothesis directly
- **Whitened embeddings**: if embeddings are anisotropic, try PCA/whitening before
  clustering — this may recover structure in the 3B space
- **Sub-HOG depth sweep**: run level-2 and level-3 sub-HOG coloring to see if deeper
  HOG levels show clustering even in the 3B space
- **Functional annotation clustering**: try GO term labels instead of HOG — do the
  models capture function better than evolutionary group?
