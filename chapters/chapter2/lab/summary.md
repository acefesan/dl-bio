# Chapter 2 Lab — Experiment Summary

Exploring how well ESM2 protein language model embeddings capture evolutionary
and functional structure (taxonomy, HOG phylogenetic groups).

## Results Table

| #   | Name                    | Model         | Key finding                                              | Date       |
|-----|-------------------------|---------------|----------------------------------------------------------|------------|
| 001 | Initial clustering      | ESM2 650M     | KMeans sil 0.04, visible taxa clusters in UMAP           | 2026-02-16 |
| 002 | 3B multi-seed UMAP      | ESM2 3B       | Near-zero KMeans sil, HOG sil -0.21; consistent 4 seeds  | 2026-03-04 |
| 003 | Subtree HOG coloring    | ESM2 650M     | Level-1 sub-HOG view of top 4 root HOGs in UMAP          | 2026-03-06 |
| 004 | Three-model comparison  | 150M/650M/3B  | **Inverse scaling**: 150M best KMeans (0.05), 650M best HOG (0.03), 3B worst (-0.21) | 2026-03-06 |
| 005 | Pretraining model scaling | 7.7M→650M scratch | No overfitting at any scale on 10K seqs; 650M performs worst; effective rank drops more in deeper models | 2026-03-07 |
| 006 | Subtree depth sweep | 150M/650M/3B | Sub-HOG levels 1–3, multi-seed (27 runs, 135 figures); awaiting visual review | 2026-03-07 |
| 007 | GO prediction        | 150M/650M/3B | **Inverse scaling in supervised task**: 650M best (auPRC 0.24), 3B worst (0.09) | 2026-03-07 |
| 008 | 1.47B pretraining    | 1.47B scratch | Inverse scaling extends to 1.47B; best eval 2.68 vs 2.56 for small models; peaks at step 2500 like 650M | 2026-03-08 |
| 009 | GO grid search       | 150M/650M/3B  | **Not classifier capacity**: best 3B (0.114) < worst 150M (0.120); depth hurts 3B | 2026-03-08 |
| 010 | 200-epoch divergence | 7.7M→1.47B scratch | 650M **diverged** (eval 2.68→4.38) due to lr/batch interaction; 7.7M–148M stable at 2.58 | 2026-03-09 |
| 011 | Uniform clustering | 150M/650M/3B | Inverse scaling confirmed on common protein set (141K); UMAP sil negative for all models; 150M best (0.052), 3B ~0 | 2026-03-15 |

## Key Result: Inverse Scaling

| Model | Avg KMeans sil (orig) | Avg KMeans sil (UMAP) | Entry |
|-------|-----------------------|-----------------------|-------|
| 150M  | **0.0518**            | **-0.067**            | 011   |
| 650M  | 0.0370                | -0.091                | 011   |
| 3B    | 0.0001                | -0.174                | 011   |

Entry 011 (uniform protein set, 4 seeds): larger ESM2 models produce worse
KMeans clusters. UMAP silhouette is negative for all models — KMeans clusters
from high-D do not form compact regions in 2D. 3B is worst on both metrics.

## Key Result: Pretraining Scaling (entry 005)

Training ESM2 from scratch on 10K UniRef sequences (MLM, 84 epochs):

| Model | Params | Best eval loss | Gen gap |
|-------|--------|---------------|---------|
| 6L-320H | 7.7M | 2.5625 | +0.015 |
| 12L-320H | 15.1M | 2.5611 | +0.015 |
| 20L-480H | 55.9M | 2.5625 | +0.015 |
| 30L-640H | 148.5M | 2.5712 | +0.013 |
| 33L-1280H | 651.7M | 2.6275 | +0.003 |
| 33L-1920H | 1,465M | 2.6755 | -0.022 |

No overfitting at any scale (7.7M → 1.47B). Larger models perform equal or worse.
Effective rank of attention layers drops more steeply in deeper models
(0.63 for 6L → 0.48 for 33L), confirming low-rank simplicity bias.

## Key Result: Supervised GO Prediction (entry 007)

| Model | Params | Test Loss | Test auPRC | Test auROC |
|-------|--------|-----------|------------|------------|
| 150M  | 537K   | 0.0742    | 0.2135     | 0.8212     |
| 650M  | 865K   | **0.0718**| **0.2389** | **0.8376** |
| 3B    | 1.5M   | 0.0934    | 0.0905     | 0.6279     |

Inverse scaling extends to supervised prediction — not just a geometry artifact.
650M > 150M here (unlike clustering), suggesting supervised methods extract more
from richer representations while unsupervised methods suffer dimensionality curse.

## Open Questions

1. **Anisotropy hypothesis**: are 3B embeddings concentrated in a narrow cone,
   making cosine distances non-discriminative? (check SVD spectrum)
2. **Whitening rescue**: can PCA/whitening recover structure in the 3B space?
3. **Sub-HOG depth**: does level-2/3 sub-HOG coloring reveal finer structure? *(entry 006 generated — awaiting review)*
4. ~~**Functional clustering**: do GO term labels cluster better than HOG labels?~~ *(answered by entry 007 — supervised GO prediction works, 650M best)*
5. ~~**Feature granularity vs geometry**: is the inverse scaling a geometry artifact?~~ *(entry 007 shows it's not — supervised prediction also shows inverse scaling)*
6. **Force overfitting**: remove weight decay or shrink dataset to <1K to see if
   memorization finally occurs
7. **Fair epoch comparison**: rerun all from scratch for 200 epochs with lr scaled by batch size *(entry 010: first attempt diverged due to resume bug + lr/batch mismatch)*
8. **Dataset cartography**: track per-sequence confidence/variability across checkpoints
9. **3B rescue**: can PCA/whitening on 3B embeddings before GO prediction improve performance?
10. **Longer training**: does 3B catch up with more steps (currently only 300)?
11. ~~**Classifier capacity**: is the 3B problem just an undersized MLP?~~ *(entry 009: NO — best 3B with 2.9M params still worse than worst 150M with 121K params)*
12. **Why depth hurts 3B**: deeper networks actively degrade 3B performance — amplifying noise in high-dim space?
