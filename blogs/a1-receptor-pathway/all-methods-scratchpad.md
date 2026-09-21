> **Scratchpad — all seven scGSA methods from Wang & Thakar 2024, run as the
> authors run them, on the minimal 8-gene A1 panel vs. the 61-gene panel.**
> Not prose. Run notes, fidelity accounting and results for post 2, section 5.

## Why this run exists

Every previous number in these scratchpads came from **one** method
(AddModuleScore) in a **simplified pseudobulk variant** — scoring each
cell-type group once from group-level means rather than scoring each cell.
That is not the implementation the paper benchmarks, so the paper's
reassurance about AddModuleScore's low false-positive rate never actually
applied to our numbers. This run replaces that with the real per-cell
implementations of all seven methods.

## The two panels

**minimal8** — receptor, neuronal transducer, and the two effector arms that
have hippocampus-specific literature behind them, plus one canonically
Gi-inhibited adenylyl cyclase:

`ADORA1`, `GNAO1`, `GNAI1`, `KCNJ3`, `KCNJ6`, `CACNA1A`, `CACNA1B`, `ADCY5`

Chosen for specificity to A1, not for coverage. Dropped on purpose:
`GNAI2/3` (ubiquitous peripheral paralogs), `KCNJ5` (cardiac GIRK4, ~2% in
CA4 anyway), `KCNJ9` (accessory subunit), the PKA subunits and `CREB1` (too
far downstream, shared with every Gs/Gi receptor in the cell), and the whole
PLC/PKC, PI3K/MAPK and NO/cGMP branches (no neuronal A1 evidence — they come
from cardiovascular/immune reviews).

Note the cost of that choice: `PLCB1`, `PRKCA` and `PRKCE` were the genes that
did the discriminating work in the thalamic-vs-CA4 comparison
(`branch-score-scratchpad.md`), and they are exactly what a specificity-first
panel has to drop, because PLCβ is canonically Gq-coupled.

**full61** — the exhaustive panel from `pathway-scratchpad-ca4.md`, unchanged,
carried along so the two can be compared directly.

## Environment

The machine migration left no R, no `uv` and no project venv, so the whole
stack was rebuilt. Conda env at `~/.local/share/mamba/envs/scgsa`
(micromamba, ~4 GB — safe to delete, see "Reproducing" below).

| Package | Ours | Paper |
|---|---|---|
| Seurat | 5.5.1 | 4.2.0 |
| AUCell | 1.28.0 | 1.18.1 |
| UCell | 2.10.1 | 2.0.1 |
| GSVA (ssGSEA) | 2.0.0 | 1.44.5 |
| irGSEA | GitHub `chuiqin/irGSEA` | 2.1.5 |
| scPS | authors' `scPS_2024.R` | — |
| JASMINE | authors' `JASMINE_V1_11October2021.r` | — |

## Fidelity: we use the authors' own invocations

The paper's repo (`github.com/Thakar-Lab/scPS`) ships the exact scripts used
for the comparison. We reproduce their calls rather than guessing:

- **Five methods via irGSEA**, exactly as `Comparison analysis/Rscript irGSEA.R`:
  `irGSEA.score(assay="RNA", slot="data", seeds=42, custom=TRUE, ncores=8,`
  `method=c("AUCell","JASMINE","ssgsea","UCell","scSE"), kcdf="Gaussian",`
  `ucell.MaxRank=NULL, aucell.MaxRank=NULL, JASMINE.method="oddsratio")`
- **AddModuleScore** with bare defaults (nbin 24, ctrl 100), as their
  `Rscript AddmoduleScore.R`.
- **scPS** sourced from their `scPS_2024.R` (PCA on scaled data,
  `weight.by.var=FALSE`, 50% cumulative-variance PC cutoff, score = weighted
  PC aggregate × gene set mean expression).

A second, independent implementation (`score_direct.R`) calls AUCell, UCell,
GSVA, JASMINE and the SCSE formula directly, as a cross-check on the wrapper.

## Deviations from the paper — read before trusting anything below

1. **200 cells per cell-type group**, sampled with seed 0; groups with
   <500 cells dropped. Rank-based methods score each cell against the whole
   gene axis, which is not tractable across 2.5M cells. 200/group is the
   paper's own main sample size (their scripts live under `Sample size 200/`)
   and the point where they report all methods stabilising.
2. **No zero imputation.** The paper runs scImpute and reports that every
   method improves under it. We ran non-imputed data only. This is exactly
   the interpretive gap flagged in the cognitive-debt list in `post.md`.
3. **scPS `npcs` clamped from 10 to 7** on the 8-gene panel. The authors
   hardcode `npcs = 10`, which cannot run on a panel of 8 genes. See the
   stability section — this is not a cosmetic change.
4. **Newer package versions** than the paper's (table above).
5. **Both atlases normalised identically from raw counts** (HBCA X is counts;
   Tabula counts pulled from `raw/X` rather than its log-normalised `X`), then
   `NormalizeData` in Seurat. This removes the cross-atlas scale caveat that
   made earlier HBCA/Tabula magnitudes non-comparable.

## HEADLINE: the panel choice flips which cell type wins

CA4's rank among 21 HBCA superclusters, same data, same methods, only the
gene panel differs (1 = top):

| Method | full61 | minimal8 |
|---|---:|---:|
| scPS | 5 | **1** |
| AUCell | 9 | **1** |
| UCell | 5 | **1** |
| ssGSEA | 10 | **1** |
| AddModuleScore | 5 | **2** |
| JASMINE | 7 | **2** |
| SCSE | 6 | **2** |

On the 61-gene panel CA4 is mid-pack and **Hippocampal CA1-3 wins instead**
(rank 1 under AUCell, AddModuleScore, JASMINE and UCell). On the 8-gene panel
CA4 is first or second under all seven methods, and CA1-3 falls to 3rd–8th.

**This inverts the earlier finding.** The old pseudobulk run reported CA4 1.94
vs. CA1-3 1.92 on the full panel — a near-tie with CA4 nominally ahead. Run
per-cell with real implementations, the full panel does not support CA4 at
all. The CA4 story survives only on the A1-specific panel.

Read together with the branch scratchpad, the likely explanation is the one
already suspected there: the 61-gene panel is dominated by ubiquitous
signalling genes (PLCβ, PKC, calmodulin, MAPK/AKT), and `PLCB1`/`PRKCA`/
`PRKCE` are abundant in hippocampus generally. The full panel was partly
measuring "how much generic signalling machinery does this neuron carry,"
which CA1-3 wins on population grounds. Stripping to the 8 A1-specific genes
is what isolates the receptor's own wiring.

## scPS is seed-unstable on an 8-gene panel

Seurat's `RunPCA` defaults to truncated SVD (irlba). Asking it for 7 of 8
possible components triggers `did not converge--results might be invalid!` on
both atlases. Consequences, over three seeds (`scps_sensitivity_hbca.csv`):

| SVD | seed 42 | seed 1 | seed 7 |
|---|---|---|---|
| `approx=TRUE` (Seurat default) | CA4 | CA4 | **Amygdala excitatory** |
| `approx=FALSE` (exact) | CA4 | CA4 | CA4 |

With the default solver the winning cell type changes with the random seed.
With exact SVD the ranking is identical across seeds. The authors' own code
has `# approx = FALSE` sitting commented out on that line — they evidently hit
this too. **Any scPS number on a small panel should use `approx=FALSE`.**
The main run used the default, and happened to agree with the stable answer
at seed 42, but that is luck rather than evidence.

Related: the authors' own example filters gene sets to **50–60 genes**
(`predBind_maxMin <- list(max = 60, min = 50)`). scPS is not designed for
panels the size of ours.

## The CA4 lead is distributional, not a few outlier cells

Per-cell scores, minimal8, CA4 (n=200) vs. all other 4,000 cells pooled:

| Method | CA4 rank by *median* | Mann-Whitney p | effect size (CLES) |
|---|---:|---:|---:|
| UCell | 1 | 4.6e-43 | 0.787 |
| AUCell | 1 | 1.9e-37 | 0.766 |
| ssGSEA | 1 | 5.1e-35 | 0.757 |
| scPS | 1 | 8.4e-31 | 0.740 |
| JASMINE | 1 | 2.9e-23 | 0.706 |
| AddModuleScore | 2 | 2.9e-27 | 0.725 |
| SCSE | 2 | 2.6e-26 | 0.721 |

CLES ≈ 0.75 means a random CA4 cell outscores a random non-CA4 cell about
three times in four. Real, consistent, but a long way from separable — worth
stating honestly in the post rather than as "CA4 is the A1 cell type."

The two methods where CA4 is second both put **Upper rhombic lip** first, and
both are the magnitude-based methods (AddModuleScore subtracts binned
controls; SCSE is a fraction of that cell's total counts). Worth one check
before writing: whether Upper rhombic lip simply has low total counts, which
would inflate SCSE's ratio.

## Implementation cross-check: the two paths agree exactly

`score_all_methods.R` (irGSEA wrapper, the authors' route) and
`score_direct.R` (AUCell / UCell / GSVA / authors' JASMINE / SCSE formula
called directly) were run independently on the same HBCA sample. Spearman rho
of group ranks was **+1.000 for all five shared methods on both panels**, with
the same top-ranked cell type every time. The wrapper adds no hidden
transformation, and the numbers below are not an artifact of one code path.

## Cross-method agreement (Spearman rho of group ranks, HBCA)

Between panels, per method — i.e. how much the panel choice moves each
method's ranking:

| Method | rho(minimal8, full61) |
|---|---:|
| ssGSEA | 0.821 |
| AUCell | 0.774 |
| JASMINE | 0.736 |
| SCSE | 0.673 |
| UCell | 0.638 |
| scPS | 0.595 |
| **AddModuleScore** | **0.392** |

AddModuleScore — the method all our previous numbers came from — is the one
whose ranking is *least* stable to the panel change.

Within the minimal8 panel, the methods split into two families: the
rank-based ones agree tightly (AUCell/JASMINE/ssGSEA rho 0.89–0.96) and the
magnitude-based ones agree with each other (AddModuleScore/SCSE rho 0.95) but
poorly across the divide (SCSE vs. ssGSEA rho 0.47). On the full61 panel the
methods agree much more (mostly 0.82–0.96). Small panels are where method
choice starts to matter — which is the paper's point, arriving from the
opposite direction to the one we expected.

## Tabula Sapiens: the panels are measuring different things entirely

105 cell types (n >= 500 cells), 200 cells sampled each.

**full61 → myeloid/immune.** Top hits across all seven methods are classical/
intermediate monocytes, neutrophils and mural cells. Cross-method agreement is
high (rho 0.81–0.97).

**minimal8 → the vascular contractile lineage.** Completely different, and
internally coherent. Ranks out of 105, consistent across methods:

| Cell type | scPS | AUCell | UCell | SCSE |
|---|---:|---:|---:|---:|
| mural cell | 1 | 1 | 1 | 1 |
| smooth muscle cell | 4 | 3 | 4 | 3 |
| pericyte | 5 | 4 | 3 | 4 |
| myofibroblast cell | 6 | 5 | 6 | 5 |
| vascular assoc. smooth muscle | 15 | 12 | 17 | 11 |
| *endothelial cell* | *50* | *50* | *46* | *56* |

Six of the top six are the same contractile mural family, and **endothelium is
not enriched** — which is the control that makes this look like signal rather
than a generic "vascular tissue" effect. This is a plausible A1 story
(adenosine acting on vascular smooth muscle and pericytes to regulate tone),
but note we have not checked it against the vascular A1 literature at all yet.

**The two panels are essentially uncorrelated here.** Spearman rho between
minimal8 and full61 group ranks: scPS +0.07, AddModuleScore +0.09,
AUCell −0.06, JASMINE −0.10, UCell −0.00, SCSE −0.05, ssGSEA +0.32. On HBCA
the same comparison ran 0.39–0.82. In peripheral tissue the 61-gene panel and
the 8-gene panel rank cell types independently of each other.

That is the clearest evidence yet for the specificity argument. The full
panel's peripheral hits are monocytes and neutrophils — exactly what you would
expect from a panel carrying `NFKB1`/`RELA`, the PI3K/MAPK arm and the
PLC/PKC arm, which are immune signalling machinery. The full panel was
substantially reporting "which cells have lots of generic signalling
proteins."

**Brain vs. periphery, now a legitimate comparison** (both atlases normalised
identically from raw counts — this is what the old cross-atlas caveat
blocked). Top HBCA score ÷ top Tabula score, minimal8: UCell 10.1x,
AUCell 5.8x, scPS 5.4x, SCSE 5.1x. On full61 the same ratios collapse to
1.0–2.0x. The A1-specific panel preserves the brain/periphery distinction that
the exhaustive panel washes out.

**This revises the earlier read.** `branch-score-scratchpad.md` concluded
Tabula was "flat, essentially nothing... no comparable enrichment." With
per-cell scoring on an A1-specific panel there is coherent structure — it is
weaker than brain, and it is vascular rather than neuronal.

ssGSEA is the method that disagrees with the others on Tabula minimal8
(rho 0.31–0.67 against the rest, vs. 0.85–0.97 among the others) — consistent
with the paper's finding that ssGSEA is the weakest of the seven.

## Outputs

In `projects/caffeine/lab/001_adora_expression/figures/`:

- `a1_allmethods_group_scores_hbca.csv` — 7 methods x 2 panels x 21 groups
- `a1_allmethods_group_scores_tabula.csv` — 7 x 2 x 105 groups
- `a1_allmethods_percell_hbca.csv` — per-cell scores, 4,200 cells
- `a1_allmethods_percell_tabula.csv` — per-cell scores, 21,000 cells
- `a1_scps_sensitivity_hbca.csv` — scPS seed/SVD stability
- `a1_allmethods_direct_hbca.csv` — independent-implementation cross-check

Scripts in the lab dir (`scgsa_*`): `scgsa_setup_env.sh` builds the R env and
fetches the authors' source; `scgsa_extract_subsample.py` writes the
subsampled counts; `scgsa_score_all_methods.R` and `scgsa_score_direct.R` are
the two implementations; `scgsa_scps_sensitivity.R` is the stability check.

## Open TODOs

- [ ] **Answer the cognitive-debt questions in `post.md` before using any of
      this in prose.** Nothing here resolves them; the zero-count question in
      particular is now more load-bearing, not less, because no imputation was
      run.
- [ ] Re-run scPS everywhere with `approx=FALSE` and treat those as the
      canonical scPS numbers.
- [ ] Check whether Upper rhombic lip's SCSE/AddModuleScore lead is a
      total-count artifact.
- [ ] Decide the panel story for the post. The honest version is now: the
      61-gene panel does **not** support CA4, the 8-gene panel does, and the
      difference is interpretable — but that needs saying out loud rather than
      quietly switching panels.
- [ ] Consider the with/without-imputation comparison the paper runs, since
      that is the single biggest untested variable here.
- [ ] **Check the mural/pericyte/smooth-muscle result against the vascular A1
      literature.** It is a clean, method-independent signal with endothelium
      as a working negative control, but we have read nothing on vascular A1
      yet — it could be well-established, or it could be an artifact of the
      panel sharing genes with vascular contractile machinery (`CACNA1A/B`,
      `ADCY5` and `KCNJ3` are not neuron-exclusive).
- [ ] Only HBCA got the independent-implementation cross-check. Tabula was
      run through the irGSEA path only — cheap to repeat with
      `scgsa_score_direct.R tabula` if any Tabula number becomes load-bearing.
- [ ] The subsample is a single draw at seed 0. Nothing here has been repeated
      across draws, so the stability of these rankings to *which* 200 cells
      were picked is untested (distinct from the scPS seed issue above).
