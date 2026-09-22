> **Scratchpad — the A1 minimal panel across the whole local CELLxGENE Census:
> 96.6M primary human cells plus three non-human primates.** Not prose. Run
> notes and results for post 2, sections 1 and 5.

## What this run is, and why it isn't scored

Local Census build 2025-11-08 at `/mnt/bulk/dl_bio/cellxgene_census/`.

**No scGSA scoring here, deliberately.** With ~100M human cells, per-cell
rank-based scoring (AUCell, UCell, JASMINE, ssGSEA) is not tractable, and the
method-comparison work concluded an 8-gene panel does not need scoring anyway
— scoring exists to compress gene sets too large to read, and this one isn't.
So this computes **exact per-gene statistics over every cell**: `n_cells`,
`pct_expressing`, `mean_all`, `mean_expressing`, grouped by `cell_type` and by
`tissue_general`. That is a stronger result than a scored subsample, not a
weaker one — no sampling, no control-gene draw, no seed.

Only `is_primary_data == True` cells are counted, since the Census duplicates
cells across overlapping datasets (158,982,719 human obs → 96,591,226 unique).

Efficiency note: the X read touches only the 8 panel gene columns, so it
returns their nonzeros (212M human entries) rather than the full
61,497-gene axis. Human raw X is 538 GB on disk and reads at ~115 MiB/s over
the FUSE/NTFS mount, so the human pass takes roughly 1.5 hours.

## Coverage

| Experiment | primary cells | cell types | panel genes |
|---|---:|---:|---:|
| human | 96,591,226 | 870 | 8/8 |
| macaque | 2,929,014 | 54 | 8/8 |
| marmoset | 1,712,738 | 32 | 8/8 |
| chimp | 158,099 | 24 | 8/8 |
| spatial human | 3,042,411 | 207 | 8/8 |
| spatial mouse | 2,645,566 | 63 | 8/8 |
| **mouse (dissociated)** | **absent** | — | — |

`census_data/mus_musculus` did not land in the sync. The group metadata still
lists it, so both a directory listing and the SOMA API's `.keys()` show mouse
— it only fails on open (`DoesNotExistError`). Being re-synced.

## Headline: ADORA1 is a deep-layer cortical projection neuron gene

Top human ADORA1 cell types by `pct_expressing` (min 5,000 cells), of 870:

| Cell type | n cells | % expressing | mean |
|---|---:|---:|---:|
| near-projecting glutamatergic cortical neuron | 61,801 | 66.9 | 2.09 |
| L5 extratelencephalic projecting cortical neuron | 16,466 | 63.5 | 1.70 |
| L5/6 near-projecting glutamatergic neuron | 39,073 | 62.2 | 1.46 |
| corticothalamic-projecting cortical neuron | 87,791 | 57.1 | 1.20 |
| indirect pathway medium spiny neuron | 11,091 | 55.0 | 1.10 |
| direct pathway medium spiny neuron | 11,276 | 54.3 | 1.07 |
| L2/3-6 intratelencephalic projecting neuron | 2,364,373 | 53.4 | 1.11 |
| L6 intratelencephalic projecting neuron | 107,182 | 52.3 | 1.14 |

By tissue, brain dominates outright: 18.7% expressing across 28,967,109 brain
cells, mean 0.335. Next tissues are far behind — spinal cord 12.2%, CNS 7.0%,
eye 4.7%, everything else under 4%.

**Both medium spiny neuron classes carry ADORA1 at ~54-55%.** MSNs are the
canonical *A2A* cell type in the story so far, so the fact that they are also
near the top for A1 is worth a line in the post — the "A1 in cortex/
hippocampus, A2A in striatum" split is cleaner in the telling than in the data.

## The hippocampus question: not contradicted, but not reproduced either

Census hippocampal entries:

| Cell type | n cells | % expressing | mean |
|---|---:|---:|---:|
| hippocampal pyramidal neuron | 12,345 | 49.1 | 0.93 |
| hippocampal granule cell | 34,220 | 25.2 | 0.33 |
| hippocampal interneuron | 9,976 | 9.0 | 0.10 |
| hippocampal astrocyte | 27,092 | 6.8 | 0.08 |

Hippocampal pyramidal neurons rank 10 by `pct_expressing` among cell types
with >= 5,000 cells, behind the deep-layer cortical types — whereas the HBCA
work put Hippocampal CA4 top.

**This is a granularity difference, not a disagreement.** The Census uses the
Cell Ontology `cell_type` axis, which has no CA4 term at all; HBCA's
`supercluster_term` does. Every CA subfield is collapsed into "hippocampal
pyramidal neuron" here, which mixes CA1/CA2/CA3/CA4. HBCA is itself a
CELLxGENE dataset, so these are overlapping, not independent, observations.
Do not present this as replication or as refutation. If the post wants a
CA4 claim it has to come from HBCA's own labels.

## The finding that most affects the post: the receptor is the variable part

Percent of cells expressing each panel gene, in the top ADORA1 cell type of
each species:

| | ADORA1 | GNAO1 | GNAI1 | KCNJ3 | KCNJ6 | CACNA1A | CACNA1B | ADCY5 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| human | 66.9 | 95.8 | 71.1 | 96.9 | 93.7 | 92.6 | 98.6 | 34.8 |
| chimp | 66.6 | 95.4 | 69.6 | 90.2 | *1.4* | 96.5 | 89.5 | 19.8 |
| macaque | 75.8 | 99.1 | 82.0 | 91.9 | 83.4 | 85.9 | 98.6 | 65.0 |
| marmoset | 68.4 | 89.3 | 56.6 | 94.3 | 91.3 | 98.9 | 94.4 | 67.4 |

The downstream machinery is close to universal in these neurons — GIRK
channels, Gi/o, and both presynaptic calcium channels sit at 85-99%. **ADORA1
itself is the limiting, variable component at 67-76%, and `ADCY5` is the other
variable one.**

This is the cleanest argument yet for the minimal-panel logic, and it partly
undercuts pathway scoring as a whole: if seven of eight genes are on in nearly
every cortical neuron, a panel score is mostly measuring those, and the one
gene that actually varies is the receptor. It also supports simply reporting
ADORA1 and treating the rest as presence/absence context.

*Chimp `KCNJ6` = 1.4% is an artifact, not biology.* It is near-absent
census-wide in chimp (0.78% of all 158K cells, max 5.2% in any cell type)
versus 11.3% / 14.4% / 50.9% census-wide in human / macaque / marmoset. The
chimp census is a single dataset on one reference build; read this as a
genome-annotation gap for `KCNJ6`.

## Cross-species: conserved across ~40M years

Top ADORA1 cell type per species, by `pct_expressing` (min 500 cells):

- **human** — near-projecting glutamatergic cortical neuron, 66.9%
- **chimp** — L5/6 near-projecting, primary motor cortex, 66.6%
- **macaque** — L6 intratelencephalic projecting, primary motor cortex, 75.8%
- **marmoset** — L6 corticothalamic-projecting cortical neuron, 68.4%

All four are deep-layer (L5/L6) cortical projection neurons, at strikingly
similar prevalence. Marmoset is a New World monkey, so this spans roughly 40
million years of primate divergence. Caveat: the three non-human primate
censuses are cortex-focused, so "top cell type is cortical" is partly
constrained by what tissue was sampled — the *prevalence* similarity is the
more defensible claim than the ranking.

## Spatial is not useful for this post

Both spatial collections are spot-based (Visium V1 + Slide-seqV2), with
`in_tissue`/`array_row`/`array_col` coordinates and `suspension_type = "na"`.
A spot spans several cells, which is why the largest human spatial `cell_type`
is `unknown` (886,855 of 3.04M).

**Neither collection contains brain tissue.** Human spatial is kidney, heart,
lung, lymph node, breast (37% cancer); mouse spatial is 76% kidney plus
embryo, six tissues total. Mouse spatial ADORA1 is low everywhere (lung 7.4%,
colon 2.1%, kidney 0.4%).

The one angle: ADORA1 in kidney was in the original Q1 hypothesis, and these
are kidney-heavy, so spatial could locate A1 within kidney (cortex vs.
medulla). Secondary to the brain story.

## Gotcha worth remembering

Mouse orthologs are `Adora1`, `Kcnj3` — not uppercase. An exact
`feature_name.isin(PANEL)` match silently returned **0/8 genes** for mouse and
would have looked like missing data. `census_a1_panel.py` now matches
case-insensitively and records the canonical symbol. This will matter again
when the dissociated mouse census lands.

## Outputs

In `projects/caffeine/lab/001_adora_expression/figures/`:
`a1_census_human_by_cell_type.csv`, `a1_census_human_by_tissue.csv`,
`a1_census_{macaque,marmoset,chimp}_by_cell_type.csv`,
`a1_census_spatial_{human,mouse}_by_tissue.csv`.
Script: `census_a1_panel.py` (`CENSUS_URI` / `CENSUS_OUT` env vars).

## Open TODOs

- [ ] Re-run once dissociated mouse lands — that is the species the A1
      electrophysiology literature actually used, so it is the one that can
      test whether expression matches the measured physiology.
- [ ] Decide how the post handles CA4 given the Census has no CA4 term. The
      honest framing is that the CA4 claim is HBCA-label-specific.
- [ ] The `L4/5 IT` entry has the highest `mean_all` (4.53) but only 37.5%
      expressing and `mean_expressing` 12.08 — far above every other type.
      Check whether that is one donor or one dataset before quoting it.
- [ ] Quantify the HBCA/Census overlap rather than hand-waving it; HBCA cells
      are inside this 96.6M, so the two analyses are not independent.
- [ ] Nothing here is disease-stratified. `disease` was read but not used;
      the Census is not all healthy tissue.
