# Blogpost Outline — Reference & Analytics Draft

**Purpose.** This is a *reference list matching the outline of the blogpost*, not the blogpost itself. You write the prose and own the understanding; this file is the index that tells you, for each section, **which source to read** (wiki page, external paper) and **which analytics/figure backs the claim**. Every source is here so you can read and verify it before you write the section.

**Legend:**
- 📖 **Read** — a wiki page or external source to understand the section.
- 📊 **Analytics/figure** — a computed result or figure in this repo you can cite or embed.
- ⚠️ **Gap** — something the outline calls for that we don't have yet (write-around or future work).

Paths are relative to this file (`projects/caffeine/lab/001_adora_expression/`).

---

## Hook — *caffeine enters your body; neurons, epithelial cells react. How do we paint this story?*

📖
- [wiki/concepts/caffeine-molecular-targets.md](../../wiki/concepts/caffeine-molecular-targets.md) — what caffeine actually binds.
- [wiki/concepts/adenosine-receptors.md](../../wiki/concepts/adenosine-receptors.md) — the receptor family the story centers on.
- [wiki/concepts/cell-type-response-model.md](../../wiki/concepts/cell-type-response-model.md) — the "different cells react differently" thesis, framed.
- (color, optional) [wiki/concepts/caffeine-in-plants.md](../../wiki/concepts/caffeine-in-plants.md), [wiki/concepts/caffeine-cultural-history.md](../../wiki/concepts/caffeine-cultural-history.md).

📊
- The whole-project framing figure to set up "paint this story": [lab/001_adora_expression/figures/ranked_top20_per_receptor.png](figures/ranked_top20_per_receptor.png) — the punchline preview (which cells light up per receptor).

---

## 1) What are ADORA receptors

### 1a) The idea behind the G-protein coupled receptor *(+ pathway sketch)*

📖 (read in this order — this is the signaling cascade end to end)
- [wiki/concepts/pharmacology-vocabulary.md](../../wiki/concepts/pharmacology-vocabulary.md) — ligand/agonist/antagonist vocabulary first.
- [wiki/concepts/gpcr.md](../../wiki/concepts/gpcr.md) — receptor architecture: ligand outside → signal inside.
- [wiki/concepts/g-protein.md](../../wiki/concepts/g-protein.md) — the GDP/GTP molecular switch.
- [wiki/concepts/g-protein-switching.md](../../wiki/concepts/g-protein-switching.md) — *(optional deep dive)* how the GPCR flips the switch (GEF mechanism).
- [wiki/concepts/g-protein-coupling.md](../../wiki/concepts/g-protein-coupling.md) — decoding Gi/o, Gs, Golf, Gq.
- [wiki/concepts/camp-signaling.md](../../wiki/concepts/camp-signaling.md) — cAMP cascade.
- [wiki/concepts/kinase.md](../../wiki/concepts/kinase.md) — PKA → CREB step.
- [wiki/concepts/signaling-to-transcription.md](../../wiki/concepts/signaling-to-transcription.md) — how the cascade becomes gene regulation.

📊 *pathway sketch* — ⚠️ no diagram asset exists yet. **Gap:** you'll draw/source the GPCR→Gα→cAMP→PKA→CREB schematic. The text for every arrow is in `camp-signaling.md` and `g-protein-coupling.md`.

### 1b) The 4 types — what each does differently

📖
- [wiki/concepts/adenosine-receptors.md](../../wiki/concepts/adenosine-receptors.md) — A1/A2A/A2B/A3 overview.
- [wiki/concepts/camp-signaling.md](../../wiki/concepts/camp-signaling.md) — the **receptor-by-receptor** coupling schematic (A1/A3 = Gi↓cAMP; A2A/A2B = Gs↑cAMP). This is the "what each does differently" core.
- [wiki/concepts/secondary-caffeine-targets.md](../../wiki/concepts/secondary-caffeine-targets.md) — higher-dose / non-ADORA targets, to bound the claim.
- [wiki/concepts/caffeine-molecular-targets.md](../../wiki/concepts/caffeine-molecular-targets.md) — caffeine as antagonist across the four.

📊
- Per-receptor expected vs observed cell-type biology: [lab/001_adora_expression/entry.md](entry.md) (Hypothesis + Interpretation sections) — your prose source for "A1 = glia/neurons, A2A = striatum/immune, A2B = epithelium/astrocytes, A3 = myeloid/mast".

---

## 2) How can we measure these receptors

### 2a) Single-cell transcriptomics — array of methods & caveats

📖
- [wiki/concepts/single-cell-rna-seq-measurement.md](../../wiki/concepts/single-cell-rna-seq-measurement.md) — what scRNA-seq actually measures (transcripts, not protein), dropout, sparsity.
- [wiki/concepts/scrna-visualization-and-analysis.md](../../wiki/concepts/scrna-visualization-and-analysis.md) — the standard pipeline + plot catalogue.
- [wiki/concepts/umap-and-dimensionality-reduction.md](../../wiki/concepts/umap-and-dimensionality-reduction.md) — embedding mechanics & what not to over-read.
- [wiki/concepts/public-data-landscape.md](../../wiki/concepts/public-data-landscape.md) — the atlases we drew on.
- Data plumbing (for the "how we got the data" aside): [wiki/concepts/census-source-h5ads.md](../../wiki/concepts/census-source-h5ads.md), [wiki/concepts/tiledb-soma-storage.md](../../wiki/concepts/tiledb-soma-storage.md), [wiki/labs/001-data-flow.md](../../wiki/labs/001-data-flow.md).

📊 **Array of methods — assays present across our datasets** (recomputed 2026-06-13, **neurons included**; total ~4.5M cells across the three real atlases):

| Assay | EFO | Tabula Sapiens | HBCA non-neuronal | HBCA neurons |
|---|---|---:|---:|---:|
| 10x 3′ v3 | EFO:0009922 | 1,025,717 (90.3%) | 888,263 (100%) | 2,480,956 (100%) |
| 10x 5′ v2 | EFO:0009900 | 67,331 (5.9%) | — | — |
| Smart-seq2 | EFO:0008931 | 41,501 (3.7%) | — | — |
| Smart-seq3 | EFO:0022488 | 1,669 (0.15%) | — | — |

Per-dataset totals: Tabula 1,136,218 · HBCA non-neuronal 888,263 · HBCA neurons 2,480,956 cells. (A 100-cell `stratified_100` probe also exists — Smart-seq2 + BD Rhapsody — ignore for analysis.)

*Caveat to write about:* the data is overwhelmingly **droplet 3′-tag (10x 3′ v3)** — and now even more so: both HBCA halves are 100% 10x 3′ v3, so the entire ~3.4M-cell brain is single-assay and the only assay heterogeneity lives in Tabula. The ~3.8% plate-based **Smart-seq2/3** (full-length) Tabula cells read sparse genes like ADORA differently. The check for this is [lab/001_adora_expression/figures/assay_stratified_dotplot.png](figures/assay_stratified_dotplot.png) (Tabula-only by necessity).

### 2b) Proteomics?

📖
- [wiki/concepts/single-cell-rna-seq-measurement.md](../../wiki/concepts/single-cell-rna-seq-measurement.md) — the "RNA is evidence of transcript, not receptor on the membrane" limit; this is the *motivation* for asking about proteomics.
- [wiki/concepts/scrna-visualization-and-analysis.md](../../wiki/concepts/scrna-visualization-and-analysis.md) — mentions **CITE-seq / totalVI** (surface-protein + RNA) as the multimodal option.

⚠️ **Gap.** There is no proteomics page or proteomics dataset in the project. Honest framing for the blogpost: scRNA gives transcript abundance; GPCR surface protein can diverge from mRNA; orthogonal validation would need antibody/CITE-seq or mass-spec proteomics, which we have not pulled. **Decide:** write this as a "limitation / what's next" subsection, or research it (could spin up a sources pass).

---

## 3) Results

### 3a) Distribution by tissue

📖
- [lab/001_adora_expression/ARTIFACTS.md](ARTIFACTS.md) §3 — how the tissue dotplot was built and how to read it.
- [wiki/labs/001-adora-interpretation.md](../../wiki/labs/001-adora-interpretation.md) — reading thresholds.

📊
- [lab/001_adora_expression/figures/tabula_sapiens_adora_dotplot_tissue.png](figures/tabula_sapiens_adora_dotplot_tissue.png) (+ `.csv`) — ADORA × 28 broad tissues.
- [lab/001_adora_expression/figures/tabula_sapiens_embeddings_by_tissue.png](figures/tabula_sapiens_embeddings_by_tissue.png) — atlas overview by tissue.
- Drill-down example: [lab/001_adora_expression/figures/tabula_sapiens_tongue_adora_cell_type_breakdown.png](figures/tabula_sapiens_tongue_adora_cell_type_breakdown.png) — "tongue → basal cells" tissue→cell-type move.

### 3b) Distribution by cell-type

📖
- [lab/001_adora_expression/entry.md](entry.md) — Interpretation section (the actual per-receptor findings + cross-receptor co-expression).
- [lab/001_adora_expression/ARTIFACTS.md](ARTIFACTS.md) §2, §4, §5 — the Q1 figures, HBCA brain figures, and the joint Tabula×HBCA comparison.
- [wiki/labs/001-joint-umap.md](../../wiki/labs/001-joint-umap.md) — the cross-atlas embedding method.

📊
- **Headline:** [lab/001_adora_expression/figures/ranked_top20_per_receptor.png](figures/ranked_top20_per_receptor.png) — top 20 cell types per receptor (Tabula + HBCA).
- [lab/001_adora_expression/figures/tabula_sapiens_adora_dotplot_cell_type.png](figures/tabula_sapiens_adora_dotplot_cell_type.png) — cell-type dotplot.
- [lab/001_adora_expression/figures/hbca_adora_dotplot_cell_type.png](figures/hbca_adora_dotplot_cell_type.png) — brain non-neuronal.
- Tables to quote numbers from: `cache/pseudobulk_by_cell_type.feather` (191 cell types × per-gene mean/pct), `cache/cross_receptor_overlap.feather` (multi-receptor cells — the microglia quadruple-positive finding).
- Sanity checks: [lab/001_adora_expression/figures/donor_stratified_dotplot.png](figures/donor_stratified_dotplot.png), [assay_stratified_dotplot.png](figures/assay_stratified_dotplot.png).

> **Coverage caveat to state plainly:** the brain figures currently shown are the **non-neuronal** half (glia/vascular). The neuronal half (HBCA neurons, 2,480,956 cells) is now **downloaded** (`cache/human_brain_cell_atlas/hbca_all_neurons_8e10f1c4.h5ad`) but **not yet folded into the Q1 figures/pseudobulk** — that regeneration is the next compute step before this section can claim the A1/A2A neuronal story.

---

## Open gaps this outline exposes (for your "what's next" / honesty section)
- ⚠️ **Pathway diagram (1a)** — no asset; needs drawing.
- ⚠️ **Proteomics (2b)** — no page, no data; decide limitation vs research.
- ⏳ **Brain neurons (3b)** — downloaded (2.48M cells); assay count updated, but Q1 figures/pseudobulk not yet regenerated to include them (currently glia-only).
- Future questions queued: GTEx bulk cross-check, Q5 (CYP1A2 liver), Q9 (chromatin at ADORA loci). See [wiki/concepts/research-questions.md](../../wiki/concepts/research-questions.md).
