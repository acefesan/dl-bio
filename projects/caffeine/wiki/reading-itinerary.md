# Tonight's Reading Itinerary

A single front-to-back path through everything Lab 001 has done, ordered the way the work actually happened: **get the data → learn the embedding → see the joint embedding we built → read every figure → know the open questions.** Read top to bottom; each stop links forward to the next. Roughly 60–90 minutes end to end.

For the biology-first ramp (caffeine → receptors → signaling) read [learning-path.md](learning-path.md) instead. This itinerary is the *data-and-methods* path.

## 1. What We Downloaded And Why (≈15 min)

The project is an integration problem over public atlases, so it starts with getting the right files onto disk and the false starts that taught us how.

1. [public data landscape](concepts/public-data-landscape.md) — the atlases in play and which are references vs perturbations.
2. [001 CellxGene Census API](labs/001-cellxgene-census-api.md) — the API we pull from.
3. [TileDB-SOMA storage](concepts/tiledb-soma-storage.md) — *why the obvious "few genes, many cells" query is the worst case*, the thing that cost us two overnight runs.
4. [Census source H5ADs](concepts/census-source-h5ads.md) — the escape hatch we switched to: download whole per-dataset H5ADs instead of querying SOMA.
5. [001 fetch stall post-mortem](labs/001-fetch-stall-postmortem.md) — what actually went wrong, in order.
6. [001 data flow](labs/001-data-flow.md) and [001 H5AD and AnnData cache](labs/001-h5ad-anndata-cache.md) — the files that ended up on disk (45 GB Tabula, 4.4 GB HBCA non-neuronal) and how we read them memory-safely with h5py.

**You'll be able to say:** what Tabula Sapiens and HBCA are, why we have them as local H5ADs rather than live SOMA queries, and what's still not downloaded (HBCA neurons, 30 GB).

## 2. UMAP Fundamentals (≈15 min)

Before any picture, the machinery that makes the picture.

1. [single-cell RNA-seq measurement](concepts/single-cell-rna-seq-measurement.md) — what the numbers in the matrix even are.
2. [UMAP and dimensionality reduction](concepts/umap-and-dimensionality-reduction.md) — the normalize → select genes → PCA/SVD → batch-correct → UMAP chain, what UMAP actually computes (two-phase neighbor graph + force layout), the two knobs (`n_neighbors`, `min_dist`), and **which features of a UMAP are real and which are artifacts you must not over-read.**
3. [scRNA visualization and analysis](concepts/scrna-visualization-and-analysis.md) — the full plot-type catalogue (dotplot, stacked violin, matrixplot, heatmap, density UMAP) and the standard pipeline each figure came from.

**You'll be able to say:** why cluster *distance* and *size* in a UMAP carry no information, why a sparse gene like ADORA needs the overdraw fix, and what a dotplot separates that a UMAP can't.

## 3. The Joint UMAP We Built (≈10 min)

The cross-atlas embedding is its own method — read it as method, not just a figure.

1. [001 joint UMAP](labs/001-joint-umap.md) — exactly how `make_shared_umap_tabula_hbca.py` puts Tabula and HBCA on one set of axes: shared-gene intersection (56,999 → 2,504), stratified 40k+40k sampling, identical CPM-10k+log1p normalization, TruncatedSVD to 50 components, the `center` vs **Harmony** batch correction (Harmony is the one to trust), and the final cosine UMAP. Includes how to read the dataset-colored QC plot and the ADORA-high overlay, plus the four output-directory variants.

**You'll be able to say:** why you can't just stack two published UMAPs, what batch correction is removing, and what the Harmony joint plot does and doesn't prove.

## 4. Every Figure, And How It Was Made (≈25 min)

The guided gallery. This is the long, rewarding part.

1. [ARTIFACTS.md](../lab/001_adora_expression/ARTIFACTS.md) — the walkthrough of *every* file in `figures/` and `cache/`, each with **what it shows / how to read it / what it's good for / caveats.** Eight sections, from source data through the story-telling Q1 figures, the original Tabula figures, the HBCA brain figures, the Tabula×HBCA comparison set, and the Q1 summary tables.
2. [001 ADORA interpretation](labs/001-adora-interpretation.md) — the companion that reads the current outputs: thresholds, dotplot reading, the sparse-UMAP overlays, and the missing-brain-coverage caveat.
3. [entry.md](../lab/001_adora_expression/entry.md) — the lab notebook entry with the receptor-by-receptor Interpretation section (the actual Q1 findings).

**You'll be able to say:** for any PNG in the figures folder, what question it answers and where it can mislead you.

## 5. The Questions This Feeds (≈10 min)

Where Q1 points next.

1. [research questions](concepts/research-questions.md) — the full slate of 20 questions.
2. [cell type response model](concepts/cell-type-response-model.md) — the four-layer prediction stack Q1 is layer one of.
3. [caffeine sensitivity genetics](concepts/caffeine-sensitivity-genetics.md) — the pharmacokinetic (CYP1A2/AHR) vs pharmacodynamic (ADORA2A) split behind several questions.

**You'll be able to say:** how the cell-type-by-receptor map gates the chromatin (Q9), NFAT/vascular (Q10), and CellOracle (Q13) questions that come after it.

---

**The one gap to keep in mind throughout:** HBCA neurons (2.5M cells, 30 GB) are not downloaded, so the strongest ADORA1/2A prior — neuronal expression in striatum, cortex, hippocampus, cerebellum — is not yet testable. Everything brain-related in the current figures is the *non-neuronal* (glial/vascular) half of the story.
