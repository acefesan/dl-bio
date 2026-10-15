# Dataset candidates for the adenosine / caffeine vasodilation question

Written 2026-10-15 by the session asked to scout (request from dl-bio).
Question: which vascular cell types and organs (brain first, heart second) express
ADORA2A/ADORA2B and the vasodilation machinery (17-gene set), and does that match
what the literature says caffeine does to vessels?

**Verification tags.** `[verified]` = I read the paper's data-availability text or the
record this session. `[from memory]` = accession/size recalled from background
knowledge, not checked; confirm before relying on it. Nothing was downloaded.

Ranked by how directly each would improve on the Census for this question.

---

## Tier 1: purpose-built brain-vascular atlases (replace the Census brain block)

### 1. Winkler et al. 2022, Science: normal + AVM human brain vasculature
- **What:** scRNA-seq of vessel-enriched brain tissue, plus bulk RNA-seq of AVMs. Human.
- **Size:** 181,388 cells total; normal cortex 5 individuals / 74,535 cells; AVM 5 patients / 106,853 cells; 39 bulk AVMs. `[verified]`
- **Resolves:** endothelial (arterial / capillary / venous), pericyte, smooth muscle, fibroblast, perivascular macrophage; labels are the study's own, so fine arterial-venous zonation is available. Use the normal cortex only.
- **Confounders recorded:** surgical tissue from living patients (no PMI, no agonal state, which is the main gain), demographics in Table S1. Caveat: the "normal" cortex is resected tissue from neurosurgery patients, not healthy volunteers. `[verified]`
- **Access:** processed data via UCSC cell browser <https://adult-brain-vasc.cells.ucsc.edu> (open); raw reads dbGaP phs002624.v2.p1 (controlled); code <https://github.com/cnk113/vascular-analysis>; Zenodo <https://zenodo.org/records/5826205> holds only the code zip (12.9 MB). `[verified]`
- **Why better:** one lab, one protocol, native vascular labels with zonation, no post-mortem hypoxia confound for ADORA2B.

### 2. Garcia et al. 2022, Nature: single-cell dissection of the human brain vasculature (BVE-seq)
- **What:** snRNA-seq with a vascular-enrichment protocol (BBB nuclei enrichment). Human.
- **Size:** 16,681 cerebrovascular nuclei, 11 subtypes (from abstract); 17 fresh temporal-lobe resections (epilepsy), 8 HD + 8 controls post-mortem, 12 ROSMAP controls across 7 regions. `[verified]`
- **Resolves:** endothelial, mural (pericyte / SMC), three perivascular fibroblast subtypes; in-silico sorted vascular nuclei from post-mortem ROSMAP cortex across 7 regions.
- **Confounders recorded:** fresh-vs-post-mortem is an explicit split, so you can *test* the PMI confound; PMI used as a covariate in the HD analysis (per-sample values probably in supplementary tables, not read). `[verified]`
- **Access:** GEO **GSE173731**, count matrices on <http://compbio.mit.edu/scBBB/> (open). `[verified]`
- **Why better:** the only candidate that lets us compare fresh surgical vs post-mortem vessels directly and see whether ADORA2B/hypoxia genes move with PMI.

### 3. Yang et al. 2022, Nature: human brain vascular atlas (VINE-seq)
- **What:** snRNA-seq of vessel-isolated nuclei, hippocampus and cortex, post-mortem. Human.
- **Size:** 143,793 nuclei; 25 samples; 17 subjects (cognitively normal and AD). (from search summary of abstract)
- **Resolves:** arterial / capillary / venous endothelium, two pericyte subtypes, SMC, perivascular vs meningeal fibroblasts.
- **Confounders recorded:** diagnosis, age, sex, region; PMI in supplementary metadata (not checked). Post-mortem only, so the agonal/hypoxia confound is present.
- **Access:** GEO **GSE163577** ("The human brain vasculature at single-cell resolution", Homo sapiens; checked against GEO 2026-10-15 by dl-bio), and the lab's data portal. Open.
- **Why better:** deeper vascular sampling and richer pericyte subtyping than the Census labels; pairs with #1 and #2 as a three-study replication set, ideally harmonised with a batch-aware method.

### 4. Mouse brain vascular atlas: Vanlandewijck et al. 2018 (Betsholtz) and Zhao/Zeisel-derived
- **What:** scRNA-seq of mouse brain vessels and mural cells, with arteriovenous zonation. Mouse.
- **Size:** about 3,500 vascular and perivascular cells in the original. `[from memory]`
- **Resolves:** endothelial zonation (arterial → capillary → venous), arteriolar SMC, pericytes, vSMC; a well-known continuous-zonation model.
- **Confounders:** none of the human ones (controlled lab animals, fresh tissue); no human PMI/agonal-state issue.
- **Access:** <http://betsholtzlab.org/VascularSingleCells/database.html> and GEO **GSE98816** (mouse brain vascular scRNA-seq; checked 2026-10-15). GSE99235 exists but is the mouse **lung** vascular set, not brain. Open.
- **Why better:** a clean ground truth for *which* vascular segment carries Adora2a/2b, free of post-mortem effects. Caffeine/adenosine can also be tested in vivo in mice.

---

## Tier 2: bulk RNA with donor covariates (controls the confound the Census can't)

### 5. GTEx v8/v10 bulk RNA-seq (and snRNA-seq)
- **What:** bulk RNA-seq across about 50 tissues, post-mortem. Human. Also a GTEx snRNA-seq release (8 tissues, 2022).
- **Size:** about 950 donors, about 17,000 samples (v8). `[from memory]`
- **Resolves:** tissue level only: coronary artery, aorta, tibial artery, brain regions (cortex, hippocampus, caudate, etc.), heart left ventricle / atrial appendage. No cell types, but vessel-rich arteries act as an SMC/endothelium proxy.
- **Confounders recorded:** ischemic time, death circumstance on the Hardy scale (0–4, ventilator vs sudden vs slow), PMI, age, sex, BMI, RIN, cause of death. This is the field set the Census lacks.
- **Access:** open expression (GTEx portal, <https://gtexportal.org>); individual-level genotype and reads are controlled (dbGaP phs000424). TPM matrix is about 1–2 GB; per-gene extraction via the API is KB.
- **Why better:** lets us regress ADORA2B, NOS3 and the rest on Hardy scale and ischemic time, quantifying exactly how much of the "signal" is hypoxia.

### 6. Human Protein Atlas: RNA consensus + single-cell type + brain regions
- **What:** consensus RNA from HPA and GTEx plus single-cell-type clusters; protein (IHC) for tissues. Human.
- **Resolves:** brain regional RNA and cell-type clusters (endothelial, pericyte etc.) and tissue-level IHC for ADORA2A/2B. Gene pages exist for both: ADORA2A `ENSG00000128271`, ADORA2B `ENSG00000170425`. `[verified that pages exist; expression values not read]`
- **Confounders:** summarised, so few; useful as a sanity check rather than a primary analysis.
- **Access:** open, downloadable TSVs (<https://www.proteinatlas.org/about/download>), small.
- **Why better:** adds a *protein-level* line of evidence, which RNA-only Census data cannot.

---

## Tier 3: perturbation data (what caffeine / adenosine actually does)

### 7. Boye et al. 2024, eLife: caffeine response regulatory variants in vascular endothelial cells
- **What:** MPRA plus caffeine treatment of human vascular endothelial cells (HUVEC), 665 variants with allelic effects, 29 caffeine-response variants. Human, in vitro.
- **Access:** GEO **GSE221514** (FASTQ + counts), open. `[verified from search result]`
- **Why:** an MPRA is a reporter assay, not transcriptomics, so it will not give ADORA expression. It does show caffeine acting on endothelial regulatory elements. Check the authors' earlier paired RNA-seq (Findley et al. 2019, caffeine-treated endothelial cells; accession not found, `[unverified]`).

### 8. Gap: adenosine-agonist / A2A-A2B-agonist (NECA, CGS21680, BAY 60-6583) or caffeine-treated vascular cells
- Searched GEO-style terms and found no dataset that I could confirm. **Recommended action:** search GEO directly for "NECA", "CGS21680", "caffeine HUVEC", "adenosine HBMEC", filtered to Homo sapiens + expression profiling by high-throughput sequencing; also the LINCS L1000 / CMap (<https://clue.io>), which holds adenosine-receptor ligands in HUVEC and smooth-muscle-like lines. `[from memory]`

---

## Tier 4: spatial and physiology

### 9. Spatial brain atlases with vessels
- Human DLPFC Visium (spatialLIBD, 12 samples, 3 donors; Maynard 2021), Allen / BICAN MERFISH of human and mouse brain, 10x Xenium brain panels. `[from memory]`
- Value: validates that vascular expression is positioned at vessels, not ambient RNA. Panels may not include all 17 genes, so check panel contents before investing.
- Access: open, spatialLIBD via Bioconductor; Allen via <https://portal.brain-map.org>. Size ranges from MB to hundreds of GB.

### 10. Caffeine and cerebral blood flow (ASL MRI)
- Published results to cross-check: caffeine lowers global CBF about 20–30% (e.g. a 26% decrease with ASL reported in a 2021 ISMRM abstract; a 3 T ASL study of 9 volunteers at 200 mg). `[verified from search snippets]`
- Searched for an OpenNeuro caffeine ASL dataset and could not confirm one; **search OpenNeuro directly** for "caffeine" and "ASL". Do not rely on a specific accession from me.

### 11. Pharmacology
- Receptor affinities: IUPHAR/BPS Guide to Pharmacology (<https://www.guidetopharmacology.org>) and ChEMBL for caffeine at A1/A2A/A2B/A3 (Ki values). `[from memory]` Small, open.
- PET occupancy: A1 and A2A occupancy studies after caffeine exist in the literature (e.g. with [11C]TMSX / SCH442416); these are papers, not datasets. `[from memory]`

---

## Suggested order of work
1. Pull #1, #2, #3 processed matrices (small, open for #2 and #1's browser); score the 17-gene set per vascular cell type per study, report per-study agreement. Garcia's fresh-vs-post-mortem split gives the PMI test.
2. Pull the GTEx per-gene table for the 17 genes and regress on Hardy scale / ischemic time.
3. Cross-check with the mouse atlas (#4) and HPA (#6).
4. Treat #7–#8 as the caffeine-perturbation leg; needs a proper GEO search.

## Open items I could not settle
- GTEx sizes: from memory, unverified. (Yang GSE163577 and Vanlandewijck GSE98816 were checked against GEO by dl-bio on 2026-10-15.)
- Whether any public caffeine-treated vascular RNA-seq exists (#8) and any OpenNeuro caffeine ASL set (#10).
- Per-sample PMI for Winkler and Garcia (supplementary tables not read).
