> **OUTLINE — no prose sections written yet.** Research and analysis now exist
> for several items below; the checklist describes what still has to be
> evaluated or incorporated into the post. See the README for current results.

# A1 Receptor: Distribution and Signaling Pathway (outline)

Second post in the series. Picks up where the motivation post left off —
zooms in on one receptor (A1) and goes deep instead of wide.

## 0. Cognitive debt — BLOCKING, answer before writing any of this

> These are questions for *me*, not TODOs to delegate or to look up once
> and paste in. The post does not get written until I can answer each one
> in my own words. Reference material is in this folder (`lqae124.pdf`,
> Wang & Thakar 2024) — but having read it is not the same as being able
> to answer these.

- [ ] **How do you properly interpret a zero count for a gene?** When
      does a zero mean "this gene is genuinely not expressed in this
      cell" versus "the assay missed it"? What would I need to know about
      the data to tell those two apart — and can I ever tell them apart
      for a *single* cell, or only in aggregate across many cells?
- [ ] **How does that relate to how the paper treats this kind of work?**
      It runs zero imputation (scImpute) and reports that every method
      improved under it. Why does imputing change the answer at all, what
      assumption is imputation making, and do I believe that assumption
      for our data? And what does it signify that the paper treats
      zero-handling as a *benchmark variable* rather than as a
      preprocessing detail?
- [ ] **Why are bigger gene sets shadier, according to the paper?** The
      naive intuition says more genes means more signal and more
      statistical power. Scenarios 3 and 4 show the opposite — large sets
      produced *higher* false-positive rates even with no
      condition-specific signal present at all. What is the actual
      mechanism driving that, and does it apply to our 61-gene panel?

Why these are load-bearing rather than academic: our current data claims
rest on co-expression across a 61-gene panel, and on readings where
low-detection genes (`NOS2/3`, `PLCB2/3`, `PREX1`, `KCNJ5`) were treated
as possibly-absent. Both of those are exactly the interpretations these
three questions govern. Answer them first, or the post asserts things I
can't defend.

## 1. Where is A1 expressed?

- [x] Pull the A1-specific cell-type distribution from the existing
      artifacts. The export and interpretation are in
      `distribution-scratchpad.md`; Census-scale results and their different
      cell-type resolution are in `census-scratchpad.md`.
- [ ] By-tissue view — which figure/table shows this cleanly for A1
      specifically (as opposed to all four receptors at once)?
- [x] Identify the by-cell-type source table — `adora1_top_cell_types.csv`.
- [ ] Decide how to present the HBCA table alongside the broader Census
      result, which uses different cell-type labels.
- [ ] Decide whether we need a new, A1-only figure (isolating just the one
      receptor from the existing multi-receptor dotplots) or whether the
      existing ones are legible enough to reference directly.

## 2. The full signaling pathway, molecule by molecule

- [ ] Enumerate every molecule/step in the chain, not just the headline
      players — adenosine → A1 receptor → Gi/o (α and βγ arms separately)
      → [α arm] adenylyl cyclase inhibition → lower cAMP → less active PKA
      → fewer phosphorylated targets (HCN channels, CREB, ...) and
      [βγ arm] GIRK channel opening + presynaptic calcium channel
      inhibition → direct hyperpolarization / reduced neurotransmitter
      release.
- [ ] For each step, name the actual gene/protein involved (not just the
      generic role) so we can check whether it's present in our
      single-cell data.
- [ ] Note which parts of the chain are well-established vs. which are
      inferred/generic (same rigor bar as the citation work on post 1).

## 3. Representative cell type — walk the pathway end to end

- [ ] Pick one cell type to fully narrate first (candidate: hippocampal
      CA4 neurons, given it's our strongest A1 signal — confirm this is
      the right pick before committing).
- [ ] Walk every step of the pathway specifically in that cell type,
      citing literature where the mechanism is documented for that cell
      type specifically (not just "neurons in general").

## 4. Repeat across the other A1-relevant cell types

- [ ] Same walk-through for: thalamic excitatory neurons, deep-layer
      cortical neurons, oligodendrocytes, OPCs, astrocytes.
- [ ] Flag where the pathway differs by cell type (e.g. oligodendrocytes/
      OPCs likely don't use the same GIRK/glutamate-release mechanism as
      neurons — the "differentiation via cAMP" story from our cAMP
      discussion is a different downstream consequence entirely).

## 5. Check the pathway against our actual data

- [x] Extract and score downstream pathway genes in HBCA and Tabula; see
      `pathway-scratchpad-ca4.md`, `all-methods-scratchpad.md`, and the lab
      artifact guide. The 8-gene and 61-gene panels give different rankings.
- [ ] Define a defensible mini-methods section for the completed analysis,
      including panel selection, cell sampling, scoring, zeros, and the
      expression-matched random-set control still to be run.
- [ ] Be explicit about what this can and can't prove — co-expression of
      pathway genes is suggestive, not functional proof the pathway is
      active in that cell type (same epistemic caution as the "can our
      data corroborate this" discussions from post 1's research phase).

## Open questions / not yet decided

- [ ] Exact scope boundary for "every single molecule in the pathway" —
      how far downstream do we go before it stops being about A1
      specifically (e.g. do we include general PKA biology, or just the
      A1-specific entry points into it)?
- [ ] Whether this post also revisits/answers the lingering questions from
      post 1 (why coffee still "hits" after good sleep, individual
      sensitivity differences) or keeps those for a later post.
