> **Scratchpad — per-branch module scoring, both atlases.** Not prose.
> Answers "does each cell type express *which* pathway branch," which the
> single combined score in `module-score-scratchpad.md` could not.

## Method

Same AddModuleScore-style approach (`score_a1_branches.py`), but the
61-gene panel is split into its 10 documented branches and each is scored
separately, against control genes matched to *that branch's* genes.
Also dumps the full per-gene × per-cell-type mean table.

Outputs:
- `figures/a1_branch_scores_hbca_neurons.csv` / `..._tabula_sapiens.csv`
- `figures/a1_pergene_means_hbca_neurons.csv` / `..._tabula_sapiens.csv`

## ⚠️ Read the table DOWN columns, not ACROSS rows

**Branch scores are raw expression-difference units, and branches differ
enormously in baseline expression magnitude.** In CA4: `PLCB1` mean is
85.7 while `ADORA1` mean is 1.8. So the PLC/PKC branch scoring 13.38 vs
the receptor scoring 1.04 does **not** mean PLC/PKC is "13× more
enriched" — the units aren't comparable. Comparing cell types *within* a
branch (down a column) is valid; comparing branches *within* a cell type
(across a row) is not, unless the scores are normalized first (not done).

## ⚠️ The cAMP arm scores negative almost everywhere — treat as suspect

`Gi_alpha`, `adenylyl_cyclase`, and `PKA` are negative in ~20 of 21
superclusters. A branch being negative in essentially *every* cell type
is more likely a systematic method artifact than biology: with 25
quantile bins over 60,606 genes (~2,400 genes/bin), control genes drawn
from a bin can sit well above their matched panel gene in expression,
biasing scores down. **Do not interpret the sign of these three branches
as "the cAMP arm is depleted"** without re-running with finer bins /
Seurat's full defaults (24 bins, 100 controls per gene).

## HBCA neurons — what actually holds up

Comparing down columns, the hippocampus leads the branches with real signal:

| Branch | Top cell types |
|---|---|
| `Gbg_GIRK` | **CA4 4.04**, dentate gyrus 3.97, CA1-3 3.93, midbrain-derived inhibitory 2.83 |
| `Gbg_Ca_channels` | **CA4 11.49**, amygdala excitatory 10.91, thalamic excitatory 9.38 |
| `PLC_PKC_Ca` | **CA4 13.38**, CA1-3 10.98, medium spiny 6.92 |
| `receptor` (ADORA1 alone) | **CA4 1.04**, deep-layer near-projecting 0.81, thalamic excitatory 0.71 |

**CA4 tops four of the four branches that carry signal.** That's a
stronger result than the combined score gave — it's not one branch
dragging an average up.

Note the combined column changed vs. the earlier run (CA4 1.51 here vs.
1.94 before, and Amygdala excitatory jumped to #2). That's expected, not
a bug: this "combined" is the mean of *branch* scores (each branch
weighted equally regardless of gene count), whereas the earlier one
weighted each *gene* equally. Worth picking one convention deliberately.

## Resolved: the Thalamic excitatory anomaly (was an open TODO)

The per-gene table explains it concretely. Raw mean expression:

| Gene | CA4 | Thalamic exc. | note |
|---|---:|---:|---|
| `ADORA1` | 1.8 | 1.5 | receptor present in both |
| `CACNA1A` / `CACNA1B` | 24.7 / 26.6 | 19.8 / 20.7 | Ca²⁺ channels present in both |
| `KCNJ6` (GIRK2) | 22.8 | **3.2** | GIRK2 largely absent |
| `PLCB1` | 85.7 | **21.8** | 4× lower |
| `PRKCA` | 63.1 | **0.6** | effectively absent |
| `PRKCE` | 61.0 | 10.1 | 6× lower |

**So it's real, not an artifact: thalamic excitatory neurons carry the
receptor and the presynaptic calcium-channel machinery, but are largely
missing GIRK2 and the PLC/PKC effector arm that hippocampal neurons are
loaded with.** Receptor presence ≠ same downstream wiring. This is a
genuinely interesting result for the post — arguably more interesting
than CA4 winning, because it's a concrete case where "expresses ADORA1"
and "can run the A1 pathway the way the hippocampus does" come apart.

Same logic explains **Mammillary body** (worst combined, −0.83):
`PLCB1` = 1.1 (vs. CA4's 85.7) and `ADORA1` = 0.2 — depleted across the
board, not just one branch.

Bonus observation: **Medium spiny neurons have essentially no GIRK**
(`KCNJ3` 0.6, `KCNJ6` 0.7) despite high `PLCB1` (77.8) — the A2A cell
type is wired differently from the A1 hippocampal cells.

## Tabula Sapiens — no meaningful signal, as expected

Filtered to groups with n ≥ 500 (105 of 161 cell types survive; the
earlier unfiltered top-15 was small-n noise, as suspected). After
filtering, **every combined score is 0.03–0.08 — flat, essentially
nothing.** Top entries are myeloid/immune (mononuclear phagocyte,
neutrophil, NK T cell, monocyte, macrophage, mast cell) but at
magnitudes indistinguishable from noise.

Reasonable read: the A1 pathway signal is a brain story; peripheral
tissue shows no comparable enrichment. Caveat: Tabula X is
log-normalized while HBCA is raw counts, so the *magnitudes* are not
cross-comparable with the HBCA table anyway — the claim here is only
"flat within Tabula," not "smaller than HBCA."

## Open TODOs

- [ ] Re-run with Seurat's full defaults (24 bins / 100 controls per
  gene) as a sensitivity check, specifically to test whether the
  all-negative cAMP-arm columns are an artifact.
- [ ] Decide the weighting convention for any "combined" number
  (gene-equal vs. branch-equal) and use it consistently.
- [ ] Normalize branch scores (e.g. z-score within branch across cell
  types) if the post wants to compare branches within a cell type.
- [ ] CA4 vs CA1-3 framing question still open — but note CA4 now wins
  each individual branch, which strengthens the case for CA4 over the
  earlier 0.02 near-tie.
