# A1 Receptor: Distribution and Signaling Pathway

**Status:** outline only — no prose written. **Blocked on purpose**, see
"What's blocking" below.

**Post file:** [`post.md`](./post.md)

Second post in the caffeine/ADORA series (follows
[`../caffeine-adenosine-receptors/`](../caffeine-adenosine-receptors/)).
Deep-dives the A1 receptor specifically: where it's expressed, then a
full step-by-step walk through its signaling pathway, cell type by cell
type, checked against what our data in
[`projects/caffeine`](../../projects/caffeine/) actually shows.

## What's blocking

Section 0 of [`post.md`](./post.md) is a **cognitive-debt list** — three
questions that have to be answerable in my own words before any of this gets
written:

1. How do you properly interpret a zero count for a gene?
2. How does that relate to the paper's treatment of zero imputation as a
   benchmark variable rather than a preprocessing step?
3. Why are bigger gene sets *shadier*, per the paper's scenarios 3 and 4?

Reference material for these is `lqae124.pdf` (Wang & Thakar 2024) — sitting
in the main checkout, not committed. None of the analysis below resolves these
questions; if anything the zero-count one is now more load-bearing, because
nothing was run with imputation.

## The two gene panels

**minimal8** — `ADORA1`, `GNAO1`, `GNAI1`, `KCNJ3`, `KCNJ6`, `CACNA1A`,
`CACNA1B`, `ADCY5`. Chosen for specificity to A1.

**full61** — the exhaustive panel from
[`pathway-scratchpad-ca4.md`](./pathway-scratchpad-ca4.md), every A1 branch
published anywhere in the body.

## Scratchpads, in the order they were written

| File | What it settled |
|---|---|
| [`distribution-scratchpad.md`](./distribution-scratchpad.md) | Where A1 is expressed (section 1 groundwork) |
| [`pathway-scratchpad-ca4.md`](./pathway-scratchpad-ca4.md) | The 61-gene panel by branch; raw CA4 co-expression |
| [`module-score-scratchpad.md`](./module-score-scratchpad.md) | First module scoring, both atlases |
| [`branch-score-scratchpad.md`](./branch-score-scratchpad.md) | Per-branch scoring; the thalamic-vs-CA4 dissociation |
| [`all-methods-scratchpad.md`](./all-methods-scratchpad.md) | All seven scGSA methods, run properly per-cell |
| [`census-scratchpad.md`](./census-scratchpad.md) | Whole-Census: 96.6M human cells + three primates |

## Where the analysis currently stands

**The panel choice, not the method choice, decides the answer.** Run per-cell
with all seven methods from Wang & Thakar 2024: on **full61** hippocampal CA4
is mid-pack (ranks 5–10) and CA1-3 leads; on **minimal8** CA4 is first or
second under all seven. Two independent implementations agree exactly
(rho +1.000). The likely cause is that the 61-gene panel is dominated by
ubiquitous signalling genes and was partly measuring generic signalling
abundance. Details in [`all-methods-scratchpad.md`](./all-methods-scratchpad.md).

**This superseded the earlier read.** Every number before that run came from
one method (AddModuleScore) in a simplified pseudobulk variant that the paper
never benchmarked.

**At Census scale the framing shifts again.** Across all 96.6M primary human
cells, ADORA1 reads as a *deep-layer cortical projection neuron* gene, and the
same holds in chimp, macaque and marmoset. Most importantly: the downstream
machinery sits at 85–99% wherever ADORA1 is expressed while the receptor
itself is at 67–76% — **the receptor is the only variable part of the panel**,
which argues for reporting A1 gene-by-gene rather than as any panel score.
Details in [`census-scratchpad.md`](./census-scratchpad.md).

**The CA4 claim is HBCA-label-specific.** The Census `cell_type` ontology has
no CA4 term and collapses all CA subfields into "hippocampal pyramidal
neuron", and HBCA sits inside the Census anyway — so the Census neither
replicates nor refutes CA4. The post has to say this rather than quietly
picking whichever framing is convenient.

## Data and code

Results and scripts live in
[`projects/caffeine/lab/001_adora_expression/`](../../projects/caffeine/lab/001_adora_expression/)
— see its [`ARTIFACTS.md`](../../projects/caffeine/lab/001_adora_expression/ARTIFACTS.md)
sections 8 and 9 for a per-file walkthrough, and
[`entry.md`](../../projects/caffeine/lab/001_adora_expression/entry.md) for the
lab-level interpretation.

Reproducing the method comparison needs an R toolchain that is not part of
this repo's normal `uv` setup: run `scgsa_setup_env.sh` in that directory.
The Census work needs `tiledbsoma`, which has no Python 3.14 wheels — use 3.12.

## Next steps

- [ ] Answer the three cognitive-debt questions. Everything else is secondary.
- [ ] Re-run once dissociated mouse lands in the local Census — mouse is
      absent from the current sync, and it is the species nearly all the A1
      electrophysiology was measured in.
- [ ] Re-run scPS with `approx = FALSE` and treat those as canonical; the
      default truncated SVD is seed-unstable on an 8-gene panel.
- [ ] Check the Census mural/pericyte/smooth-muscle result against vascular A1
      literature before treating it as a finding.
- [ ] Decide the CA4-vs-cortex framing for section 1, given the two analyses
      answer at different granularities.

## Housekeeping

All of the above is on branch **`worktree-cognitive-debt-outline`**, pushed but
**not merged to `main`** — a fresh checkout of `main` shows none of it.
