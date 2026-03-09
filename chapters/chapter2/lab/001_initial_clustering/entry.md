# 001 — Initial ESM2 650M Clustering & UMAP

**Date:** 2026-02-16
**Model:** ESM2 650M (`facebook/esm2_t33_650M_UR50D`, 1280D)
**Status:** complete

## Hypothesis
ESM2 embeddings capture evolutionary and functional signal — proteins from the same species
or homologous group (HOG) should cluster together in embedding space.

## Setup
- 142k CAFA3 protein sequences, max_seq_length=5000, batch_size=32
- Sampled 500 proteins from each of 22 taxa with ≥500 proteins → 11k total
- KMeans (k=20) in original 1280D space (cosine metric)
- UMAP (n_neighbors=15, min_dist=0.1, cosine) for visualization
- Colored by taxa and by root HOG ID

## Results

| Metric | Value |
|--------|-------|
| KMeans silhouette (1280D, cosine) | 0.0396 |
| HOG silhouette (1280D, cosine) | 0.0313 |

Weak but non-trivial positive silhouette scores. UMAP shows visible taxa clusters — species
do group loosely. HOG coloring shows less obvious structure; most proteins in the sample
have no assigned root HOG (roothog_id == 0), limiting interpretability.

## Figures
- `figures/umap_by_taxa.png` — UMAP colored by taxon (22 species, 500 each)
- `figures/umap_by_hog.png` — UMAP colored by top 15 root HOGs; majority are gray (no HOG)
- `figures/taxa_hog_heatmap.png` — HOG distribution across taxa

## Interpretation
650M embeddings have weak but real taxonomic signal. The HOG coloring is limited because
the coloring is by root HOG, which is very coarse, and most proteins have no root HOG in
the sample. The fact that many proteins lack a root HOG means we're only seeing a fraction
of the evolutionary structure in the HOG plot.

## Next steps
- Try a larger model (3B) — does more parameters capture stronger signal?
- Try smaller model (150M) for comparison
- Go deeper in the HOG tree: instead of root HOGs, look at sub-HOGs within major families
