> **Scratchpad — module-score run, both atlases.** Not prose. Data +
> interpretation notes for post 2, answering the open TODO from
> `pathway-scratchpad-ca4.md` ("no baseline comparison cell type").

## Method

AddModuleScore-style (Tirosh et al. 2016, *Science* — same method behind
Seurat's `AddModuleScore` / Scanpy's `score_genes`), simplified to run at
the **cell-type/pseudobulk level** rather than per-cell, for speed. See
docstring in `score_a1_module_pathway.py` for the exact algorithm (25
expression bins, ~20 matched control genes per panel gene). Same 61-gene
panel as the CA4 co-expression pull (`pathway-scratchpad-ca4.md`).

Benchmark context (from PMC11420841, read this session): AddModuleScore
performs reasonably above ~50 genes (we're at 61 — thin margin, not a big
cushion) but is one of several methods shown to be fooled by
condition-specific expression imbalances in that paper's synthetic
benchmarks. Worth a grain of salt on any single surprising result below.

## HBCA neurons — full ranked table

`figures/a1_module_score_hbca_neurons.csv` (21 superclusters, all ranked,
no filtering needed — smallest group here is >10k cells)

**Top 5:**

| Group | n cells | module score |
|---|---:|---:|
| Hippocampal CA4 | 10,654 | **1.94** |
| Hippocampal CA1-3 | 74,979 | **1.92** |
| Eccentric medium spiny neuron | 40,144 | 1.35 |
| Upper-layer intratelencephalic | 455,006 | 0.99 |
| Midbrain-derived inhibitory | 126,782 | 0.98 |

**Bottom 3:**

| Group | n cells | module score |
|---|---:|---:|
| Thalamic excitatory | 85,546 | 0.02 |
| Cerebellar inhibitory | 14,411 | -0.10 |
| Mammillary body | 16,602 | -1.45 |

**Headline: CA4 answers the open TODO.** Relative to an expression-matched
random background, CA4 (and CA1-3 right behind it) are genuinely elevated
for the full A1 pathway panel, not just "has the machinery like every
neuron does" — this is the real signal the raw co-expression pull
(previous scratchpad) couldn't distinguish on its own.

**Flag before writing anything: CA4 and CA1-3 are a near-tie (1.94 vs.
1.92), and both are dramatically ahead of everything else.** Worth
deciding whether the post's "representative cell type" framing should
be CA4 alone or "the hippocampus (CA1-4)" as a region-level story —
picking CA4 specifically over CA1-3 by 0.02 is not obviously
defensible as a meaningful distinction.

**Real discrepancy to resolve, not yet explained:** Thalamic excitatory
was our **#2 ranked cell type for ADORA1 alone** (`entry.md`: mean 2.41
among expressing cells, 60.2%) but scores near **zero (0.02)** on the
full pathway panel — and Mammillary body scores the single worst
(-1.45) despite reasonable ADORA1 presence. Two non-exclusive
possibilities, neither checked yet:
- These cell types express ADORA1 itself but comparatively little of
  the downstream machinery (a real biological asymmetry — worth
  actually looking at which specific panel genes are low there).
- An artifact of the simplified pseudobulk-level scoring (vs. true
  per-cell scoring) or the control-gene sampling for this particular
  region — can't rule this out without re-running with different
  random seeds / more control genes per bin.

## Tabula Sapiens — top 15 of 161 cell types

`figures/a1_module_score_tabula_sapiens.csv` (full ranked table, all 161
Tabula cell types)

**Do not trust the top of this ranking as-is.** Overall scores are much
lower than HBCA (top score 0.21 vs. HBCA's 1.94 — different scale,
expected since Tabula X is log-normalized not raw counts, another
reminder of the cross-atlas scale caveat from the RNA-seq interpretation
deck). More importantly: several top "hits" have **single- or
double-digit cell counts** (pancreatic PP cell n=2, pancreatic A cell
n=48, retinal bipolar neuron n=61) — exactly the small-sample regime
PMC11420841 warns produces unreliable scores (their own scPS method
needs 50+ cells/group to be trusted; no reason to assume AddModuleScore
is more robust here).

## Open TODOs

- [ ] **Filter Tabula results to a minimum cell count (e.g. n >= 500)
  before drawing any conclusion from it** — rerun the ranking with that
  filter applied rather than treating the current top-15 as meaningful.
- [ ] **Investigate the thalamic/mammillary discrepancy** — pull the
  per-gene breakdown for those two groups from the raw pseudobulk data
  (same approach as the CA4 gene table) to see which specific branches
  are driving the low score, before deciding whether it's biological or
  a scoring artifact.
- [ ] Decide CA4 vs. CA1-3 framing given the near-tie.
- [ ] Consider re-running with Seurat's actual defaults (24 bins / 100
  controls per gene) instead of the lighter 25/20 used here, as a
  sensitivity check on the CA4/CA1-3 vs. everything-else gap — if the
  gap survives with full defaults, that's stronger evidence it's real
  and not a control-sampling artifact.
