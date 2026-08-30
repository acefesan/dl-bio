> **Scratchpad — groundwork for post 2, section 1 ("where is A1
> expressed?").** Not prose. Just pointers to what already exists in the
> repo, plus a pulled-out A1-only table so we're not squinting at
> 4-receptor figures when we only care about one.

## Existing figures (all 4 receptors — A1 is one line/dot among four)

| Figure | What it shows | Path |
|---|---|---|
| `ranked_top20_per_receptor.png` | Top 20 cell types per receptor, 3-atlas combined (Tabula + HBCA non-neuronal + HBCA neurons) — the main one | `projects/caffeine/lab/001_adora_expression/figures/ranked_top20_per_receptor.png` |
| `tabula_sapiens_adora_dotplot_cell_type.png` | Cell-type dotplot, Tabula Sapiens (peripheral tissues only, no brain) | `.../figures/tabula_sapiens_adora_dotplot_cell_type.png` |
| `hbca_adora_dotplot_supercluster_term.png` | Cell-type dotplot, HBCA neurons, grouped by supercluster | `.../figures/hbca_adora_dotplot_supercluster_term.png` |
| `hbca_adora_dotplot_cell_type.png` | Cell-type dotplot, HBCA non-neuronal (glia/vascular/immune) | `.../figures/hbca_adora_dotplot_cell_type.png` |
| `tabula_sapiens_adora_dotplot_tissue.png` | By-tissue view (28 broad tissues), Tabula Sapiens only — no brain in this one | `.../figures/tabula_sapiens_adora_dotplot_tissue.png` |

Full walkthrough of all of these: `projects/caffeine/lab/001_adora_expression/ARTIFACTS.md`.

## Master data table (all 4 receptors, all 3 atlases, per cell type)

`projects/caffeine/lab/001_adora_expression/cache/pseudobulk_by_cell_type.feather`
— 212 rows (one per cell type per source atlas), columns are
`{RECEPTOR}_mean`, `{RECEPTOR}_pct_expressing`, `{RECEPTOR}_mean_nonzero`
for each of ADORA1/2A/2B/3. This is the actual source of the numbers in
`entry.md`'s Interpretation section and of `ranked_top20_per_receptor.png`.
It's a binary feather file (not human-browsable directly) — load with
`pandas.read_feather(...)`.

## A1-only table, pulled out and sorted (new export, this session)

Exported to: `projects/caffeine/lab/001_adora_expression/figures/adora1_top_cell_types.csv`
(166 cell types with nonzero A1 signal, all 3 atlases, sorted by
`ADORA1_mean_nonzero` descending). Top 15:

| Source atlas | Cell type | n cells | mean (expressing cells) | % expressing |
|---|---|---:|---:|---:|
| HBCA neurons | Hippocampal CA4 | 10,654 | 2.58 | 68.4% |
| HBCA neurons | Thalamic excitatory | 85,546 | 2.41 | 60.2% |
| HBCA neurons | Deep-layer corticothalamic and 6b | 78,396 | 2.03 | 56.2% |
| HBCA neurons | Deep-layer near-projecting | 18,856 | 2.01 | 62.3% |
| HBCA neurons | Deep-layer intratelencephalic | 228,467 | 1.86 | 53.2% |
| HBCA neurons | Hippocampal CA1-3 | 74,979 | 1.86 | 52.0% |
| HBCA neurons | Amygdala excitatory | 109,452 | 1.78 | 43.6% |
| HBCA neurons | Miscellaneous | 25,071 | 1.77 | 34.2% |
| HBCA neurons | Lower rhombic lip | 52,650 | 1.73 | 52.1% |
| HBCA neurons | Upper-layer intratelencephalic | 455,006 | 1.54 | 35.7% |
| HBCA neurons | Hippocampal dentate gyrus | 67,533 | 1.48 | 34.3% |
| HBCA neurons | Medium spiny neuron | 152,189 | 1.48 | 31.9% |
| HBCA neurons | Eccentric medium spiny neuron | 40,144 | 1.46 | 30.8% |
| HBCA non-neuronal | oligodendrocyte precursor cell | 105,734 | 1.33 | 23.4% |
| HBCA non-neuronal | oligodendrocyte | 494,966 | 1.31 | 20.2% |

(Two Tabula Sapiens rows — `glial cell`, `platelet` — also appear near the
top of the raw sort but with tiny `n_cells` (20 and 494) and sub-1%
`pct_expressing`; almost certainly small-sample noise, same call already
made in `entry.md`. Excluded from the table above; still in the full CSV
if we want to double check.)

## Open / not done yet

- [ ] No by-tissue-only A1 table pulled out yet — the tissue dotplot
  (`tabula_sapiens_adora_dotplot_tissue.png`) has no brain, so "by tissue"
  for A1 is really just "by broad Tabula tissue label," which undersells
  A1's story (it's mostly a brain receptor per the table above). Worth
  deciding whether "by tissue" in post 2 should really mean "by brain
  region" (HBCA has `ROIGroupCoarse`/`dissection` columns — see the
  region-provenance note in ARTIFACTS.md / the recap from earlier in this
  chat) instead of the Tabula tissue labels.
- [ ] Haven't decided which of these cell types is "the representative
  one" to walk the pathway through first (post 2 outline, section 3) —
  Hippocampal CA4 is the top hit by the numbers above, but worth deciding
  deliberately rather than just taking rank #1.
