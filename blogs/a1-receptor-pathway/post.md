> **OUTLINE — nothing below is written yet.** This is a structural skeleton
> only: section headers and bullet points describing what needs to be
> researched, analyzed, and written. Every bullet is a TODO for us, not a
> claim or a draft sentence. Do not read anything here as finished prose.

# A1 Receptor: Distribution and Signaling Pathway (outline)

Second post in the series. Picks up where the motivation post left off —
zooms in on one receptor (A1) and goes deep instead of wide.

## 1. Where is A1 expressed?

- [ ] Pull the A1-specific slice out of our existing artifacts (tissue +
      cell-type distributions) rather than re-running the analysis —
      `projects/caffeine/lab/001_adora_expression/ARTIFACTS.md` and
      `entry.md` already have this for A1 (hippocampal CA4, thalamic
      excitatory, deep-layer cortical neurons; oligodendrocytes, OPCs,
      astrocytes).
- [ ] By-tissue view — which figure/table shows this cleanly for A1
      specifically (as opposed to all four receptors at once)?
- [ ] By-cell-type view — same question.
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

- [ ] For each cell type above, does our single-cell data show expression
      of the *downstream* pathway components too (not just the ADORA1
      gene itself) — e.g. relevant adenylyl cyclase isoforms, PKA
      subunits, GIRK channel genes, CREB?
- [ ] This is new analysis, not something we've already run — scope it as
      its own mini-methods section (which genes to query, which atlas(es),
      what "signs of the chain" operationally means — some minimum
      co-expression threshold?).
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
