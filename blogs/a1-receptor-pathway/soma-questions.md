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
