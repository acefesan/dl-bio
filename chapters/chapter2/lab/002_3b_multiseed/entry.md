# 002 — ESM2 3B Multi-Seed UMAP

**Date:** 2026-03-04
**Model:** ESM2 3B (`facebook/esm2_t36_3B_UR50D`, 2560D)
**Status:** complete

## Hypothesis
A larger model (3B vs 650M) might capture stronger biological signal, producing tighter
clusters. Multi-seed UMAP (4 seeds) rules out initialization artifacts.

## Setup
- Same CAFA3 dataset, max_seq_length=3000 (to fit GPU memory), batch_size=4
- Same sampling: 22 taxa × 500 proteins = 11k total
- Seeds tested: 42, 123, 456, 789
- GPU: RTX 5090 (32GB VRAM); embedding took several hours with checkpointing

## Results

| Seed | KMeans sil (2560D) | HOG sil (2560D) |
|------|---------------------|-----------------|
| 42   | 0.0003              | -0.2131         |
| 123  | 0.0024              | -0.2153         |
| 456  | 0.0016              | -0.1988         |
| 789  | -0.0008             | -0.2203         |

Stunning reversal: 3B is *worse* than 650M. KMeans silhouette near zero, HOG silhouette
strongly negative. Consistent across all 4 seeds — not a UMAP initialization artifact.

## Figures
- `figures/umap_by_taxa_seed42.png` — seed 42 UMAP (no visible taxon clusters)
- `figures/umap_by_taxa_seed123.png` — seed 123 UMAP (same pattern)
- `figures/umap_by_taxa_seed456.png` — seed 456 UMAP (same pattern)
- `figures/umap_by_hog_seed42.png` — HOG coloring seed 42

## Interpretation
The 3B model's embedding space is not organized by taxonomy or HOG in the way 650M is.
Possible explanations:
1. **Over-parameterization / anisotropy**: 3B embeddings may live in a degenerate subspace
   where the distances used by UMAP/KMeans don't reflect biological similarity
2. **Truncation effect**: max_seq_length=3000 vs 5000 for 650M changes the protein set
   (253 sequences skipped for 150M vs more for 3B)
3. **Training objective scale**: larger ESM2 may have learned more fine-grained sequence
   features that don't correspond to coarse taxonomic/HOG groupings

HOG silhouette is negative, suggesting proteins in the same HOG are *farther apart* than
average in 2560D space — this is the opposite of what we expect and points to anisotropy
or some systematic issue with the 3B embedding geometry.

## Next steps
- Run 150M model — does the trend hold (smaller=better for this clustering task)?
- Check embedding geometry: are 3B embeddings highly anisotropic? (check singular value
  distribution, cosine distance distribution)
- Go deeper in HOG tree (entry 003): root HOGs may be too coarse; sub-HOG coloring might
  reveal structure even if root HOGs don't
