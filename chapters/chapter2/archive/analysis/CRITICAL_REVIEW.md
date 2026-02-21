# Critical Review: HOG/Protein Embedding Analysis

**Date:** 2026-02-14  
**Reviewer:** Subagent Analysis  
**Status:** ⚠️ SIGNIFICANT METHODOLOGICAL CONCERNS

---

## Executive Summary

After thorough review, this analysis has **serious methodological issues** that undermine its conclusions. The silhouette scores are misleading, the hierarchy analysis is too shallow, and key metrics are computed incorrectly. This document provides a detailed critique and actionable recommendations.

---

## 🔴 Critical Issues Found

### 1. The "Perfect ARI" Problem (BUG)

**What was found:** All pairwise Adjusted Rand Index values = 1.0

```
Pairwise ARI (Mean: 1.000)
```

**Why this is wrong:** The analysis computes ARI between samples using **HOG labels** (which are deterministic database IDs), not actual cluster assignments. When the same protein appears in two samples, it obviously has the same HOG label in both! This is circular reasoning.

**What ARI should measure:** Whether the SAME proteins get assigned to the SAME clusters when you re-run the analysis with different random seeds.

**The fix:** Compute ARI on the KMeans cluster assignments (or UMAP-derived clusters), not on the original HOG labels.

### 2. Negative KMeans Silhouette Scores

**What was found:**
```
KMeans Silhouette: -0.277 ± 0.074
```

**Why this is damning:** A negative silhouette score means points are, on average, **closer to clusters they DON'T belong to** than to their own clusters. This indicates:

1. UMAP's 2D projection destroys cluster structure
2. The clusters don't actually exist in the data
3. KMeans with k=20 is inappropriate for this data

**Interpretation:** The KMeans clustering in UMAP space is **worse than random assignment**.

### 3. Weak HOG Silhouette (0.30)

**What was found:**
```
HOG Silhouette: 0.302 ± 0.026
```

**Why this is concerning:**

| Silhouette Score | Interpretation |
|------------------|----------------|
| > 0.7 | Strong clustering |
| 0.5 - 0.7 | Reasonable structure |
| 0.25 - 0.5 | Weak, overlapping clusters |
| < 0.25 | No substantial structure |

A score of 0.30 indicates **weak, overlapping clusters**. Proteins within the same HOG are NOT tightly grouped in embedding space.

### 4. Shallow Hierarchy Analysis

**What was done:** Analyzed only 18 root HOGs with ≥100 proteins

**What was missed:** The HOG hierarchy goes 14+ levels deep!

```
Root HOG 801468: 2,236 proteins
Depth distribution:
  Level 1:  541 proteins
  Level 2:  433 proteins  
  Level 3:  366 proteins
  ...
  Level 14:   1 protein
```

**The missed opportunity:** The analysis treats all proteins in a root HOG as equivalent. But sub-HOGs represent finer evolutionary distinctions:
- `HOG:E0801468` = ancient ancestral protein (LUCA level)
- `HOG:E0801468.10h.66e.89g...` = recent divergence (more similar)

Proteins at deeper hierarchy levels should cluster MORE tightly in embedding space. **This was never tested.**

### 5. UMAP Silhouette Fallacy

**The problem:** Silhouette scores are computed on UMAP's 2D projection, not the original 640-dimensional ESM2 embeddings.

**Why this matters:** UMAP is a non-linear dimensionality reduction optimized for visualization. It:
- Distorts distances
- Creates artificial clusters
- Destroys global structure

**The right approach:** Compute silhouette in the ORIGINAL embedding space (640D), or use intrinsic clustering metrics.

---

## 🟡 Questionable Assumptions

### 1. "Root HOG = Evolutionary Group"

**Assumption:** Proteins sharing a root HOG are "orthologous" and should cluster together.

**Problem:** Root HOGs trace back to LUCA (Last Universal Common Ancestor) - 4 billion years ago! Proteins can:
- Retain similar function but diverge in sequence
- Acquire new functions while retaining sequence
- Have convergent evolution (different origin, similar structure)

**Reality check:** LUCA-level groupings are too broad. A protein kinase from bacteria and humans share an ancestor but have enormously different sequences.

### 2. "Embedding Similarity = Evolutionary Relatedness"

**Assumption:** ESM2 embeddings capture evolutionary relationships.

**Problem:** ESM2 was trained on masked language modeling of protein sequences. It captures:
- Sequence patterns (amino acid co-occurrence)
- Structural features (secondary structure)
- Local context (domain signatures)

It was NOT trained to capture:
- Evolutionary distance
- Phylogenetic relationships
- HOG-level groupings

**Better question:** Do ESM2 embeddings correlate with sequence identity? This was never tested!

### 3. "Taxa and HOG Are Independent Signals"

**Finding:**
```
Taxa Silhouette: -0.124
HOG Silhouette: +0.033
ARI (Taxa vs HOG): 0.003
```

**Interpretation:** Neither taxa NOR HOGs form coherent clusters. The near-zero ARI shows they capture different (poor) signals.

**The real finding:** ESM2 embeddings don't cleanly separate ANYTHING at these granularities.

---

## 🟢 What Was Actually Learned

Despite the issues, some genuine insights emerge:

1. **Sampling stability is high** - Different random samples give consistent (though weak) silhouette scores
2. **Taxonomic signal exists** - Taxa centroids separate mammals from fungi/protists in UMAP
3. **Large HOGs are heterogeneous** - High intra-HOG variance shows these groups are functionally diverse
4. **Coherence < 1 is good news** - Most HOGs are "tighter" than their distance from the global mean

---

## 📊 Recommended Metrics (Instead of Silhouette)

| Metric | What It Measures | Good For |
|--------|------------------|----------|
| **Davies-Bouldin Index** | Ratio of within-cluster to between-cluster distances | Comparing clustering algorithms |
| **Cophenetic Correlation** | How well dendrogram preserves pairwise distances | Hierarchical clustering |
| **Normalized Mutual Information** | Information shared between clustering and labels | Comparing to ground truth |
| **DBCV** | Density-based cluster validation | Non-spherical clusters |
| **Sequence Identity Correlation** | Does embedding distance ∝ sequence divergence? | THE key biological question |

### The Missing Analysis: Embedding Distance vs Sequence Identity

```python
# What should have been done:
from scipy.stats import spearmanr

# For pairs of proteins in the same HOG:
embedding_distances = [...]  # Cosine distance in embedding space
sequence_identities = [...]  # % sequence identity

correlation = spearmanr(embedding_distances, sequence_identities)
# This tells you if embeddings capture evolutionary divergence!
```

---

## 🔬 HOGs for Deeper Analysis

The current analysis only looked at root HOGs (LUCA level). Here are protein families with **rich evolutionary trees** worth exploring:

| Root HOG | Max Depth | N Proteins | Description |
|----------|-----------|------------|-------------|
| 801468 | 14 levels | 2,236 | Large, deeply structured family |
| 802573 | 15 levels | ~200 | Very deep hierarchy |
| 1027273 | 15 levels | ~100 | Maximum depth found |
| 1028112 | 11 levels | 123 | Medium size, deep tree |
| 1027819 | 11 levels | 121 | Similar structure |

### Recommended Analysis: Hierarchy-Depth Clustering

For HOG 801468, proteins should cluster as:
```
Level 1 proteins  → Loose cluster (ancient divergence)
Level 5 proteins  → Tighter cluster (moderate divergence)  
Level 10 proteins → Tight cluster (recent divergence)
```

**Test:** Does silhouette score INCREASE as you go deeper in the hierarchy?

---

## 📖 Beginner-Friendly Explanation

### What is a HOG (Hierarchical Orthologous Group)?

Imagine a family tree, but for proteins instead of people:

```
                    LUCA (4 billion years ago)
                         ↓
              ┌─────────────────────┐
              ↓                     ↓
         Bacteria              Eukaryotes
              ↓                     ↓
         ┌────┴────┐         ┌─────┴─────┐
         ↓         ↓         ↓           ↓
      E.coli    Bacillus   Fungi      Animals
                              ↓           ↓
                           Yeast      ┌───┴───┐
                                      ↓       ↓
                                    Mouse   Human
```

A **HOG** groups proteins that share a common ancestor. The **root HOG** goes all the way back to LUCA (the Last Universal Common Ancestor). Sub-HOGs capture more recent branching points.

### Why Should Clustering Matter?

The promise of protein embeddings (like ESM2) is that they capture biological meaning:
- Similar sequences → Similar embeddings
- Similar function → Similar embeddings
- Related evolution → Similar embeddings

**The test:** If proteins in the same HOG are truly "related," they should cluster together in embedding space.

### What Do These Results Actually Mean?

**The sobering truth:** With a silhouette score of 0.30, proteins in the same HOG are only **weakly grouped** in ESM2 embedding space.

This could mean:
1. **ESM2 captures function, not phylogeny** - Two proteins can do the same job with very different sequences
2. **HOGs are too broad** - LUCA-level groupings span 4 billion years of evolution
3. **Embeddings need better analysis** - Maybe linear methods (like silhouette) miss non-linear relationships

---

## ✅ Action Items

### Immediate Fixes

1. **Fix the ARI calculation** - Compare cluster assignments, not HOG labels
2. **Compute silhouette in original space** - Use 640D embeddings, not UMAP projection
3. **Add sequence identity correlation** - The most biologically meaningful metric

### Extended Analysis

4. **Analyze sub-HOG hierarchy** - Does clustering improve at deeper levels?
5. **Test multiple k values** - What's the optimal number of clusters?
6. **Use density-based clustering** - HDBSCAN may find natural groupings

### Better Visualizations

7. **Annotate UMAP by protein function** - GO terms, domain families
8. **Color by hierarchy depth** - Not just root HOG
9. **Add hierarchical tree plots** - Show the actual HOG structure

---

## Conclusion

This analysis asks the right question: **Do protein embeddings capture evolutionary relationships?**

But the current methodology can't answer it. The silhouette scores are misleading (computed on UMAP, not embeddings), the ARI is buggy (comparing labels, not clusters), and the hierarchy analysis is shallow (only root HOGs).

The answer is probably "partially" - ESM2 captures sequence patterns that correlate with evolution, but not perfectly. A proper analysis would:

1. Correlate embedding distance with sequence identity
2. Test if deeper sub-HOGs cluster more tightly
3. Use intrinsic metrics in the original embedding space

**The good news:** The data and infrastructure are in place. These fixes are implementable.

---

*Report generated by critical analysis subagent. Recommendations should be validated before implementation.*
