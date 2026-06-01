# Workstreams

Master view across threads of investigation. Each thread has an ID (T1, T2, ...). Book-following work lives in `chapters/`; standalone projects in `projects/`. Lab entries live in `{chapter,project}/lab/NNN_name/`.

---

## Area A — ESM2 protein embeddings (chapter 2)

**Core question:** How well do ESM2 embeddings capture evolutionary (HOG) and functional (GO) structure, and why do larger models do worse on broad clustering?

### T1 — Inverse scaling (broad clustering)
**Status:** ✅ confirmed. Larger ESM2 → worse KMeans/HOG silhouette on CAFA3.
**Entries:** 001, 002, 004, 011
**Key:** 011 nails it on a uniform 141K protein set across 4 seeds. 150M best (sil 0.052), 3B worst (~0).
**Open:** is it geometry (anisotropy/cone) or representation? → see T6.

### T2 — Does inverse scaling extend to supervised tasks?
**Status:** ✅ confirmed. 3B worse on GO prediction too.
**Entries:** 007, 009
**Key:** 007 shows 650M best (auPRC 0.24), 3B worst (0.09). 009 rules out classifier capacity — best 3B MLP (2.9M params) still < worst 150M MLP (121K).
**Open:** why depth hurts 3B specifically — amplifying noise in high-dim space?

### T3 — Pretraining from scratch (mechanism probe)
**Status:** 🔄 partial. No overfitting at any scale on 10K seqs.
**Entries:** 005, 008, 010, 012_1m_scaling
**Key:** 005 tests 7.7M→650M, all converge to ~2.56 eval. 008 extends to 1.47B (2.68). 010 tried 200-epoch rerun, 650M diverged from lr/batch mismatch.
**Open:** force overfitting (<1K seqs, no weight decay); fair 200-epoch comparison with properly scaled lr.

### T4 — Sub-HOG / paralog resolution
**Status:** 🔄 partial. Sub-HOG structure visible in larger models.
**Entries:** 003, 006, 012_hog_explorer
**Key:** 006 generated 135 figures across depth sweep 1-3, awaiting visual review. 012_hog_explorer found L3 best trade-off (53 classes, CV=0.46, 19% coverage) on HOG 801468.
**Open:** systematic silhouette scoring at each sub-HOG depth.

### T5 — Mammalian focused dataset
**Status:** 🟢 active. UMAP shows all model scales separate 8 HOGs cleanly; larger models resolve paralog subtypes.
**Entries:** 013, 014
**Key:** 30 mammalian species × 8 root HOGs × 1,576 proteins. ESM2 embeddings for 4 scales. Interactive UMAP with click-to-compare HOG paths at `chapters/chapter2/lab/014_hog_umap/figures/umap_interactive.html`.
**Open:**
- Compute silhouette for 8 HOGs across all 4 models to quantify 014's visual finding
- Test if ESM2 distances recover OMA sub-HOG hierarchy
- Within-HOG distance vs species divergence time (ortholog phylogenetic signal)

### T6 — Scaling rescue via whitening
**Status:** ⚠️ orphan — result exists, not yet documented as a lab entry.
**Artifacts:** `chapters/chapter2/07_whitening_analysis.py`, `chapters/chapter2/results/whitening/whitening_results.json`
**Key (from AGENTS.md, 2026-04-11):** whitening recovers HOG silhouette to 0.0325 at 256 dims vs 0.0150 in original 3B space. Plain PCA without whitening does poorly on HOG.
**Open:** promote to a lab entry (015?); try whitened 3B as input to the GO classifier (T2).

### Cross-thread questions (Area A)

1. **Geometry vs representation:** does whitening (T6) rescue the inverse scaling in T1/T2? If yes → geometry. If no → the signal isn't there.
2. **Scale transfer:** do the 8-HOG findings (T5) hold when lifted to the broad 141K set, or is the clean separation an artifact of the curated HOG selection?
3. **Sub-HOG depth (T4) vs paralog resolution (T5):** are they probing the same thing on different data?

---

## Area B — Caffeine epigenome (standalone project)

Lives under [`projects/caffeine/`](projects/caffeine/). Not a book chapter.

**Core question:** Can public single-cell atlases + caffeine GWAS/EWAS + the Findley/Boye HUVEC ATAC-seq datasets be integrated computationally to predict caffeine's epigenomic effects across every human cell type?

### T7 — Caffeine epigenomic fingerprint across human cell types
**Status:** 📄 proposal only — no code yet.
**Proposal:** [projects/caffeine/PROPOSAL.md](projects/caffeine/PROPOSAL.md)
**Scope:** purely computational; no new experiments. Integrate ENCODE, Roadmap, GEO (Findley 2019 / Boye 2024 HUVEC ATAC-seq), HCA, GTEx, CistromeDB, JASPAR, 4D Nucleome, LINCS L1000, GWAS/EWAS catalogs.
**Headline insight driving the design:** direct caffeine-treatment epigenomic data is extremely rare (essentially just Findley/Boye in HUVECs). Strategy must *predict* caffeine response in untested cell types via: (1) adenosine receptor expression mapping, (2) chromatin accessibility at ADORA loci, (3) chromVAR/TOBIAS TF activity, (4) CellOracle in-silico perturbation.
**20 candidate questions** in the proposal, ranked beginner→publication-worthy. Suggested starting points:
- Q1: map ADORA1/2A/2B/3 expression across all HCA cell types (Scanpy, 1-week scope)
- Q5: CYP1A2 expression across liver cell subtypes (single-cell, scoped)
- Q10: reanalyze Findley/Boye HUVEC datasets focused on NFAT motif activity
**Gaps if pursued (from proposal §10):** no genome-wide ChIP-seq under caffeine, no Hi-C under caffeine, no multi-omic caffeine data, no GRN inference specific to caffeine — all are opportunities.
**Resource reality check:** full pipeline spec in proposal is 5–10 TB storage, 40–80h compute. A meaningful first deliverable (Q1 or Q5) is far smaller — decide which before starting.
**Open (before starting):**
- Pick the first question from the 20 — probably Q1, Q5, or Q10 based on data availability + tractability
- Decide whether this competes for bandwidth with Area A or runs in parallel

---

## Pointers

- `chapters/chapter2/lab/summary.md` — per-entry results table + cumulative tables (Area A)
- `projects/caffeine/PROPOSAL.md` — Area B proposal
- `docs/architecture.md` — pipelines, repo layout, data flow
- `docs/running.md` — how to run each pipeline
- `AGENTS.md` — agent-facing pointers + conventions
