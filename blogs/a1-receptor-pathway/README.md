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

Reference material is the committed [`lqae124.pdf`](./lqae124.pdf) (Wang &
Thakar 2024). The [reading companion](./scgsa-site/) has worked explanations
of all three questions, but the author has not yet answered them in their own
words in the post. The zero-count question remains especially consequential
because the A1 scoring runs did not use imputation.

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

**Panel choice strongly changes the answer.** Run per-cell
with all seven methods from Wang & Thakar 2024: on **full61** hippocampal CA4
is mid-pack (ranks 5–10) and CA1-3 leads; on **minimal8** CA4 is first or
second under all seven. Two independent implementations agree exactly
(rho +1.000). The likely cause is that the 61-gene panel is dominated by
ubiquitous signalling genes and was partly measuring generic signalling
abundance. This interpretation is provisional: HBCA superclusters differ
4.2-fold in genes detected per cell, and an expression-matched random-set
control has not yet been run. Details in
[`all-methods-scratchpad.md`](./all-methods-scratchpad.md) and the companion's
[A1 analysis page](./scgsa-site/ours.html).

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

**Direction (2026-10-14): narrow to adenosine-driven vasodilation.** It's a
short, well-characterised chain (A2A/A2B → cAMP in smooth muscle; → eNOS/NO in
endothelium) with a measurable outcome (vessel diameter, blood flow; caffeine
lowers cerebral blood flow ~20–30%) in identifiable cell types. Signalling
diagrams and gene sets for A1, A2A and vasodilation are in section 4 of
[`scgsa-site/pathways.html`](./scgsa-site/pathways.html).

Open questions to answer once the pipeline runs on SOMA data are collected in
[`soma-questions.md`](./soma-questions.md).

- [ ] Answer the three cognitive-debt questions. Everything else is secondary.
- [ ] Re-run once dissociated mouse lands in the local Census — mouse is
      absent from the current sync, and it is the species nearly all the A1
      electrophysiology was measured in.
- [ ] Re-run scPS sign-stably: exact SVD, but with variance explained
      normalised over the first 10 PCs (what the authors' default
      `approx = TRUE` computes) and a fixed PC sign convention. **Not**
      simply `approx = FALSE` — that also redefines the 50%-variance cutoff.
- [ ] Control for genes detected per cell (the paper's scenario 3 confound):
      HBCA superclusters differ 4.2x, and CA1-3/CA4 detect the most. Score
      expression-matched random 8-gene sets in the same cells; if CA4 tops
      those too, the A1 lead is the confound. See the site's "Our A1 data".
- [ ] Check the Census mural/pericyte/smooth-muscle result against vascular A1
      literature before treating it as a finding.
- [ ] Decide the CA4-vs-cortex framing for section 1, given the two analyses
      answer at different granularities.

## Reading companion site

[`scgsa-site/`](./scgsa-site/) is a static site for reading the paper
(`lqae124.pdf`) and working through all seven scoring methods from first
principles, with an in-browser simulator of the four scenarios. Served on the
tailnet at <https://acefsan-ubuntu.taila37484.ts.net:3044>. Pedagogical and
disposable; see its README for how it is served and how its JavaScript was
validated against the R packages.
