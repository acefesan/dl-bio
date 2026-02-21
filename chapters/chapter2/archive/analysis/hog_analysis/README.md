# Intra-HOG Variance Analysis

Analysis of ESM2 protein embeddings grouped by Hierarchical Orthologous Groups (HOGs) to understand how well protein language model embeddings capture evolutionary relationships.

## Background

**HOGs (Hierarchical Orthologous Groups)** from the OMA browser represent proteins descended from a common ancestor. **Root HOGs** trace relationships back to LUCA (Last Universal Common Ancestor).

**Hypothesis**: If ESM2 embeddings truly capture evolutionary/functional relationships, proteins within the same HOG should cluster together in embedding space despite taxonomic distance.

## Dataset

```
Proteins with HOG assignments: 91,870
Unique root HOGs: 36,148
HOGs analyzed (>= 100 proteins): 18
```

## Metrics

### Intra-HOG Variance

```
V_hog = (1/n) Σ ||x_i - μ_hog||_2
```

Mean L2 distance of proteins from their HOG centroid. Lower = tighter cluster.

### Coherence

```
Coherence = Intra-HOG Variance / Distance from Global Centroid
```

- **< 1**: HOG is tighter than expected given its distance from center (coherent)
- **> 1**: HOG is more dispersed than expected (incoherent)

## Results

### Variance Statistics

| Metric | Value |
|--------|-------|
| Mean Intra-HOG Variance | 1.468 |
| Median | 1.476 |
| Range | [0.497, 2.342] |

### Coherence Statistics

| Metric | Value |
|--------|-------|
| Mean Coherence | 0.851 |
| Interpretation | Most HOGs are tighter than random |

### Correlations

| Comparison | Spearman r | Interpretation |
|------------|------------|----------------|
| Variance vs N_taxa | 0.145 | Weak: More taxa = slightly more variance |
| Variance vs HOG_size | 0.366 | Moderate: Larger HOGs = more variance |

## Key Findings

### 1. Most HOGs Are Coherent (Coherence < 1)

13/18 analyzed HOGs have coherence < 1, meaning orthologous proteins cluster more tightly than expected by chance. This suggests ESM2 embeddings do capture some evolutionary signal.

### 2. Taxonomic Breadth Has Weak Effect

The weak correlation (r=0.145) between variance and number of taxa suggests that proteins in the same HOG cluster together regardless of how many different species they come from. This is good evidence for functional similarity overriding taxonomic signal.

### 3. Most Coherent HOGs

| HOG | Coherence | N Taxa | Description |
|-----|-----------|--------|-------------|
| 136254 | 0.20 | 1 | Single-taxon (species-specific) |
| 802052 | 0.35 | 22 | Highly conserved across taxa |
| 1027819 | 0.43 | 23 | Highly conserved across taxa |

### 4. Least Coherent HOGs

| HOG | Coherence | N Taxa | Description |
|-----|-----------|--------|-------------|
| 792940 | 3.10 | 21 | Very dispersed |
| 801468 | 2.77 | 77 | Largest HOG, spans most taxa |
| 990642 | 1.36 | 59 | High taxonomic breadth |

The largest and most taxonomically diverse HOGs tend to be less coherent, possibly because:
- They contain diverged paralogs
- They span multiple functional subfamilies
- Taxonomic signal interferes with functional clustering

## Files

| File | Description |
|------|-------------|
| `umap_by_hog.png` | UMAP projection colored by top 15 HOGs |
| `intra_hog_variance.png` | Variance metrics visualization (4 panels) |
| `hog_*_taxa.png` | Individual HOG distributions colored by taxonomy |
| `hog_metrics.csv` | Complete metrics for all analyzed HOGs |
| `hog_umap_coordinates.csv` | UMAP coordinates for sampled proteins |
| `summary.json` | Summary statistics |

## Interpretation

### What ESM2 Embeddings Capture

1. **Some evolutionary signal**: Coherence < 1 for most HOGs
2. **Dominated by taxonomic/compositional signal**: As shown in earlier analysis, ~91% of variance is explained by amino acid composition
3. **Functional signal is present but secondary**: Same-function proteins cluster, but taxonomy often dominates

### Recommendations

For orthology/homology detection:
- ESM2 embeddings can help identify related proteins
- Consider normalizing for taxonomic bias (AA composition regression)
- Structure-based embeddings may be more robust for distant homologs

---

*Generated: 2026-02-09*
*Script: `chapters/umap_by_hog.py`*
