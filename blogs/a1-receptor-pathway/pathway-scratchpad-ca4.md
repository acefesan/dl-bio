> **Scratchpad — groundwork for post 2, sections 2-5 (pathway walkthrough +
> data check), starting with the representative cell type: Hippocampal
> CA4.** Not prose. Data pull + interpretation notes for us to draw on
> when actually writing.

## Scope decision (made this session)

Went with the **maximal** reading of "exhaustive" — every A1 signaling
branch documented in the literature *anywhere in the body*, not just the
branches with direct hippocampal/neuronal evidence. See the two
hippocampus-confirmed papers vs. the two broader reviews in the reading
list (life-tracker, 2026-09-02) — the reviews are where the extra
branches below (PLC/PKC, PI3K/MAPK, NO/cGMP) come from, and that evidence
is mostly cardiovascular/immune tissue, not brain. Flag this when writing
the post — "exhaustive" here means "every branch anyone has published for
A1," not "every branch known to happen in a neuron."

## The gene panel, by branch

| Branch | Genes | Literature basis |
|---|---|---|
| Receptor | `ADORA1` | — |
| Gi/o alpha subunits | `GNAI1`, `GNAI2`, `GNAI3`, `GNAO1` | canonical |
| Adenylyl cyclase (Gi-inhibited isoforms) | `ADCY1`, `ADCY5`, `ADCY6`, `ADCY8` | canonical; ADCY1/8 are the Ca/CaM-sensitive, neuronally-emphasized ones |
| PKA subunits | `PRKAR1A/1B/2A/2B` (regulatory), `PRKACA/B` (catalytic) | canonical |
| cAMP downstream | `CREB1` | canonical |
| Gβγ → GIRK channels | `KCNJ3` (GIRK1), `KCNJ6` (GIRK2), `KCNJ9` (GIRK3), `KCNJ5` (GIRK4) | hippocampus-specific, CA1 (PMC4416607) |
| Gβγ → presynaptic Ca²⁺ channels | `CACNA1B` (N-type/Cav2.2), `CACNA1A` (P/Q-type/Cav2.1) | hippocampus-specific, CA3-CA1 (ScienceDirect) |
| PLC/PKC/Ca²⁺ branch | `PLCB1-4`, `PRKCA/B/D/E`, `CALM1-3` | cardiovascular/immune reviews, not neuron-specific |
| NF-κB (downstream of PLC/PKC) | `NFKB1`, `RELA` | same reviews |
| PI3K/MAPK branch | `PIK3CA/B/D`, `PIK3R1`, `PREX1`, `RAC1`, `MAPK1/3` (ERK1/2), `MAP2K1/2` (MEK1/2), `MAPK8/9/10` (JNK1-3), `MAPK11/14` (p38), `AKT1/2/3`, `VEGFA` | cardiovascular/immune reviews |
| NO/cGMP branch | `NOS1/2/3`, `GUCY1A1`, `GUCY1B1`, `PRKG1/2` | cardiovascular reviews |

61 genes total. All 61 matched a `var/feature_name` symbol in the HBCA
neuron H5AD on the first try — no aliasing issues.

## Result: co-expression in Hippocampal CA4 (n=10,654 cells)

Full table: `figures/a1_pathway_hippocampal_ca4_coexpression.csv`
Script (reproducible): `extract_a1_pathway_ca4_coexpression.py`

**The headline finding is a methodological one, not a biological one:**
almost every gene across almost every branch shows 80-99% of CA4 cells
expressing it (GIRK channels ~97-99%, calcium channels ~99%, most
PI3K/MAPK genes >80%, most PLC/PKC genes >90%). A handful of genes are the
exception — low-detection outliers like `NOS2`/`NOS3` (~1%), `PLCB2/3`
(~3-10%), `PREX1` (~9%), `KCNJ5` (~2%) — these read as either genuinely
absent isoforms or dropout, can't tell which from this alone.

**Why the high numbers don't prove much on their own:** these are mostly
core, ubiquitous cell-signaling genes (PKA, PKC, calmodulin, MAPK/ERK) that
essentially every neuron in the brain likely expresses at some level, A1
pathway or not. **We have no baseline/comparison cell type in this pull** —
without contrasting CA4 against, say, a cell type with weak A1 signal, we
can't tell whether CA4's machinery expression is elevated, average, or
just "what every neuron has." This is the single biggest open TODO before
this data point can support any claim in the post.

## Open TODOs

- [ ] **Add a comparison/baseline cell type** (or several) to the same
  pull — ideally one with strong A1 signal (already have several: thalamic
  excitatory, deep-layer cortical) and one with weak/no A1 signal, to see
  if the "machinery" genes are actually differential or just universal.
  Without this, the co-expression table above is descriptive, not
  evidence of anything CA4-specific.
- [ ] Decide whether to keep the full 61-gene exhaustive panel in the post
  itself, or lead with just the two hippocampus-confirmed branches (GIRK,
  presynaptic Ca²⁺ channels) and treat the other three branches
  (PLC/PKC, PI3K/MAPK, NO/cGMP) as an appendix/aside given they're not
  neuron-specific evidence.
- [ ] The low-detection outliers (NOS2/3, PLCB2/3, PREX1, KCNJ5) — worth
  a sentence on dropout vs. genuine absence before using any of them as a
  "this branch isn't relevant here" claim (see the RNA-seq interpretation
  deck's dropout card — same caution applies).
- [ ] No claim yet about whether the branches that *are* well-expressed
  (GIRK, N/P-Q-type Ca²⁺ channels) are the ones actually driving CA4's
  physiology specifically, vs. hippocampal CA1/CA3 literature simply being
  assumed to transfer over. CA4 (dentate hilus / mossy cells) has much
  thinner dedicated A1 electrophysiology literature than CA1/CA3 — flagged
  already in chat, worth a line in the post itself.
