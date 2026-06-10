# UMAP And Dimensionality Reduction

## Summary

A single-cell dataset is a matrix of ~1M cells × ~60k genes. You cannot look at 60k dimensions. Dimensionality reduction is the chain of steps that turns that matrix into a 2D scatter you *can* look at — and UMAP is the last link in that chain, the one that produces the picture. This page explains what each link does, what UMAP actually computes, and — critically — which features of a UMAP picture carry information and which are artifacts you must not over-read.

For where these steps sit in the full workflow, see [scRNA visualization and analysis](scrna-visualization-and-analysis.md). For the specific cross-atlas embedding this lab built, see [001 joint UMAP](../labs/001-joint-umap.md).

## The Chain, Not Just UMAP

"Making a UMAP" is shorthand for a pipeline. UMAP is step 5; steps 1–4 decide what UMAP even sees.

| Step | What it does | Why |
|---|---|---|
| 1. **Normalize** (`normalize_total` + `log1p`) | put every cell on a comparable count scale, compress the dynamic range | a cell with 50k counts and one with 2k counts should not separate just because of sequencing depth |
| 2. **Select features** (HVG, or a fixed gene panel) | keep the ~2k–3k genes that vary informatively | 60k genes is mostly noise and housekeeping; feeding all of them buries the signal |
| 3. **Linear reduction** (PCA / TruncatedSVD → ~50 components) | collapse 2k genes into ~50 orthogonal axes that capture co-variation | denoises, and makes the neighbor search in step 4 tractable |
| 4. **(Optional) batch correction** (Harmony, scVI, centering) | remove technical axes (donor, assay, atlas) from those ~50 components | otherwise cells cluster by *who sequenced them*, not by what they are |
| 5. **UMAP** | turn the ~50-D corrected space into 2 coordinates for plotting | the only step whose output you actually look at |

The single most common misconception is that UMAP "finds" the structure. It doesn't. The structure is decided by steps 1–4 — the normalization, the gene set, and especially whether you corrected for batch. UMAP just *draws* whatever neighbor relationships those steps produced. A bad UMAP is almost always a steps-1–4 problem.

## PCA / SVD vs t-SNE vs UMAP

Three reducers show up constantly; they are not interchangeable.

| Method | Kind | Preserves | Use it for |
|---|---|---|---|
| **PCA** | linear, deterministic | global variance directions; distances along top PCs are meaningful | the ~50-D denoised input to everything downstream; a fast sanity check |
| **TruncatedSVD** | linear, deterministic | same as PCA but **without mean-centering**, so it runs directly on a sparse matrix | exactly when your normalized matrix is sparse and you don't want to densify it (this lab's case) |
| **t-SNE** | nonlinear, stochastic | local neighborhoods, tightly | structured data where you want crisp local clusters; slow at large N |
| **UMAP** | nonlinear, stochastic | local neighborhoods, plus *some* global layout | the genre-default 2D picture; faster than t-SNE at million-cell scale |

PCA and SVD differ only by mean-centering. PCA centers the columns first; SVD does not. On a sparse single-cell matrix, centering destroys sparsity (every zero becomes `-mean`), which is why the practical choice on sparse data is TruncatedSVD — it gives you PCA-like components without ever materializing a dense matrix. That is the reason [the joint UMAP](../labs/001-joint-umap.md) uses `TruncatedSVD`, not `PCA`.

## What UMAP Actually Computes

UMAP runs in two phases.

**Phase 1 — build a fuzzy neighbor graph (in the ~50-D input space).**
For each cell it finds its `n_neighbors` nearest neighbors under a chosen `metric` (Euclidean on PCA, often **cosine** on integrated latents). It then assigns each edge a weight that falls off with distance, calibrated *per cell* so that every cell is connected to its nearest neighbor with full weight. The result is a weighted graph: who is near whom, locally, with strength.

**Phase 2 — lay that graph out in 2D.**
It drops all the points into 2D and runs stochastic gradient descent with two opposing forces:
- **attraction** pulls together pairs of cells that share a strong edge in the graph,
- **repulsion** pushes apart random pairs (negative sampling) so everything doesn't collapse to a point.

It iterates until the 2D layout's edge structure best matches the high-D graph. `min_dist` enters here: it sets how tightly connected points are allowed to pack (small `min_dist` → dense, clumpy clusters; large → looser, more spread out). Because phase 2 is randomized, the picture **depends on the random seed** — two runs differ in rotation, reflection, and exact blob placement even on identical input.

### The two knobs that change the picture

| Parameter | Low value | High value |
|---|---|---|
| `n_neighbors` | emphasizes fine local structure, more fragmented islands | emphasizes global structure, smoother continents |
| `min_dist` | tight, dense clumps (good for seeing discrete types) | spread-out points (good for seeing continua/trajectories) |

The lab's joint UMAP uses `n_neighbors=15`, `min_dist=0.5`, `metric="cosine"` — Scanpy-typical defaults. Tabula Sapiens' *native* published UMAP used `metric="euclidean"` on its `X_scvi` latent with internal `a≈0.583, b≈1.334`. Different inputs and metrics mean the two pictures are not directly comparable coordinate-for-coordinate — only structurally.

## What A UMAP Picture Means — And Doesn't

This is the part to internalize before reading any embedding figure in this lab.

**Carries information:**
- *Which cells are near which.* Tight local neighborhoods are real — cells in the same little blob genuinely had similar transcriptomes in the input space.
- *Whether a labeled group is coherent.* If `cell_type = microglia` lights up one compact region, that label is internally consistent.
- *Whether a sparse gene is localized.* If ADORA-high cells concentrate in one zone rather than scattering, that's a real co-localization signal (this is the whole point of the ADORA-high overlays).

**Does NOT carry information:**
- *Distance between clusters.* Two blobs far apart are not "more different" than two blobs close together. Inter-cluster distance is essentially arbitrary.
- *Cluster size / area.* A big blob is not a more important or more variable population. Area is a layout artifact of `min_dist` and point count.
- *Density.* UMAP does not preserve how dense regions were in the original space.
- *Absolute position / orientation.* Flip the seed and the whole thing rotates. There is no meaningful "up" or "left."

**The practical rule:** use UMAP to *find* questions ("these ADORA3 cells cluster together — what are they?"), then answer them with quantitative tools — dotplots, marker-gene rankings, the pseudobulk table. Never report a conclusion whose only evidence is a UMAP distance or a blob's size. See the trap table in [scRNA visualization and analysis](scrna-visualization-and-analysis.md#the-recurring-traps).

## The Sparse-Gene Overdraw Trap

ADORA receptors are expressed in <1% of cells. Color a million-cell UMAP by raw `ADORA1` and you get a near-white plot: the 99% zero-cells are drawn last and paint over the rare positive cells. Every ADORA UMAP in this lab fixes this by drawing zero/low cells as a faint grey context layer and only the above-threshold cells in saturated color. That is a *rendering* decision, not a biological one — the threshold is a visualization floor, documented per figure in [ARTIFACTS.md](../../lab/001_adora_expression/ARTIFACTS.md) and [001 ADORA interpretation](../labs/001-adora-interpretation.md).

## Lab 001 Figures That Use This Page

| Figure | What the reduction chain produced |
|---|---|
| `tabula_sapiens_embeddings_by_tissue.png` | the *same* cells in six embeddings (PCA, scVI, UMAP variants) — a direct view of how steps 3–4 change the picture |
| `tabula_sapiens_adora_high_umap_four_color.png` | the overdraw-fix UMAP, one color per receptor |
| `shared_umap_tabula_hbca/` and `_harmony/` | the cross-atlas joint UMAP — see [001 joint UMAP](../labs/001-joint-umap.md) for exactly how it was built |

## Sources

- [UMAP: Uniform Manifold Approximation and Projection (McInnes, Healy, Melville 2018)](https://arxiv.org/abs/1802.03426) — the original method paper.
- [Understanding UMAP (Coenen & Pearce, Google PAIR)](https://pair-code.github.io/understanding-umap/) — interactive demonstration of why distance and density are not preserved.
- [scanpy `tl.umap` reference](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.umap.html) — the parameters this lab sets.

Related pages: [scRNA visualization and analysis](scrna-visualization-and-analysis.md), [single-cell RNA-seq measurement](single-cell-rna-seq-measurement.md), [001 joint UMAP](../labs/001-joint-umap.md), [001 ADORA interpretation](../labs/001-adora-interpretation.md).
