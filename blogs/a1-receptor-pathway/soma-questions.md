# Questions to tackle when the pipeline runs on SOMA data

Parked while learning the paper's pipeline in the playground
(`scgsa-site/playground.html`). Answer each one once the same steps run on the
local CELLxGENE Census (SOMA) instead of the paper's naive-B cells. Node
numbers refer to the playground graph.

1. **Node 5, "Keep genes detected in the drawn cells".** The paper drops
   genes that are zero in every drawn cell. How could that filter affect the
   SOMA results? (Asked 2026-10-12.)

2. **10x chemistry and molecules per cell.** Which 10x chemistries
   (3′ v2/v3/v3.1, 5′ v1/v2, …) make up the SOMA data, and what is the
   distribution of molecules (UMIs) per cell, overall and by chemistry? How
   do people estimate how many RNA molecules a cell actually contains (i.e.
   what fraction a 10x run captures)? (Asked 2026-10-12.)

3. **Transcriptome breadth confounds cross-cell-type scores.** Gene-set
   scores compared across cell types or groups that differ in how many genes
   each cell detects behave like the paper's scenario 3, with biology rather
   than a simulation providing the difference. This is the CA4 caveat: HBCA
   neuron groups differ ~4.2× in genes detected per cell, and the A1 score
   ranking follows that gap. Rank-based scores move because the extra genes
   reshuffle each cell's ranks; share-based scores (SCSE) move because the
   cell's total grows. On SOMA: how large is the genes-per-cell gap between
   the cell types we compare, does the ADORA/A1 ranking still hold against
   expression-matched random gene sets and after depth-matched downsampling,
   and which scoring methods are least affected? (Asked 2026-10-12.)

4. **Making cross-cell-type comparisons defensible.** Not considered at
   first, but scores (and detection rates) compared across *different cell
   types* are the analyses most exposed to depth and transcriptome-breadth
   effects — be careful how they're used. Checks to run on SOMA:
   - **Depth-matched downsampling:** randomly thin every cell to the same
     number of molecules, then recompute detection and scores. If CA4 (or
     the top SOMA cell type) still leads, depth isn't the explanation.
   - **Expression-matched random gene sets as a null:** does a random 8-gene
     set with similar expression levels also rank the same cell type first?
   - **Report detection alongside genes-detected-per-cell** for every group,
     so the confound is visible next to each result.
   - **Compare within a cell type across conditions** (e.g. ADORA1 in the
     same neuron type, caffeine-exposed vs not) wherever such data exists —
     the comparison normalisation handles well.
   (Added 2026-10-12.)

5. **Narrow focus: adenosine-driven vasodilation.** (Direction set
   2026-10-14.) Caffeine's alertness runs mostly through A2A in the nucleus
   accumbens (Huang 2005; Lazarus 2011), but mapping that onto behaviour
   goes through whole circuits. Vasodilation is a short, measurable chain:
   A2A/A2B → Gs → cAMP → PKA → K⁺ channels / MLCK → smooth-muscle relaxation,
   plus endothelial A2A/A2B → eNOS → NO → sGC → cGMP → PKG. On SOMA: which
   vascular cell types (endothelial cells, vascular smooth muscle, pericytes;
   brain vs heart vs kidney vs elsewhere) detect ADORA2A/ADORA2B and the rest
   of the vasodilation gene set (diagram 4c on `scgsa-site/pathways.html`),
   with the depth checks from items 3–4. Genes not yet in our Census panel
   (NOS3, GUCY1A1/B1, PRKG1, KCNJ8, ABCC9, KCNMA1, CACNA1C, MYLK, PPP1R12A)
   need a new scan.
