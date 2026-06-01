# Pharmacology Vocabulary

## Why This Page Exists

The rest of the wiki assumes you already know words like *ligand*, *receptor*, *agonist*, and *antagonist*. This page defines those foundational terms in plain language so a non-biologist can read [caffeine molecular targets](caffeine-molecular-targets.md) and [adenosine receptors](adenosine-receptors.md) without getting stuck. It owns the "what does this basic word mean" layer; [epigenomics vocabulary](epigenomics-vocabulary.md) owns the assay and chromatin terms.

## How Drugs Talk to Cells

| Term | Plain meaning | In this project |
|---|---|---|
| ADORA | The gene-symbol prefix for human adenosine receptor genes. It expands roughly to "adenosine receptor." | ADORA1, ADORA2A, ADORA2B, and ADORA3 encode the A1, A2A, A2B, and A3 receptor proteins |
| Receptor | A protein, usually on the cell surface, that detects a specific molecule and triggers a response inside the cell. Think of it as a lock. | The four adenosine receptors (ADORA1/2A/2B/3) are the locks caffeine fits |
| Ligand | Anything that binds a receptor. The key that fits the lock. A ligand can turn the lock or just block it. | Adenosine and caffeine are both ligands for adenosine receptors |
| Binding site / pocket | The specific cavity on the receptor where the ligand fits. | Caffeine fits the same pocket adenosine uses |
| Orthosteric site | The receptor's **main** pocket — the one the natural ligand is meant to use (*ortho-* = proper place). Two molecules competing here fight over the same spot. | Adenosine's normal pocket; caffeine competes for this same orthosteric site, which is why it is a *competitive* antagonist |
| Allosteric site | A **different** pocket elsewhere on the same receptor (*allo-* = other). Binding here does not block the main site directly — it reshapes the receptor so the orthosteric site works better or worse. A side dial, not the keyhole. | Not how caffeine works, but allosteric modulators of adenosine receptors exist in research |
| Agonist | A ligand that binds **and activates** the receptor — turns the lock, fires the downstream signal. | Adenosine is the agonist; it switches the receptor on |
| Antagonist | A ligand that binds **but does not activate** — occupies the lock so the real key can't get in. | Caffeine is an antagonist; it blocks adenosine without firing the signal |
| Competitive antagonist | An antagonist that fights for the **same** binding site as the agonist, so enough agonist can out-compete it. | Caffeine is competitive; high adenosine can displace it |
| Affinity | How tightly a ligand grips its receptor. High affinity = binds at low concentrations. | Caffeine has highest affinity for A1 and A2A, which is why those matter most at coffee doses |
| Selective | Acts on essentially one receptor subtype. | A drug hitting only A2A would be A2A-selective |
| Non-selective | Acts across several subtypes at once. | Caffeine is non-selective: it blocks A1, A2A, and more, so its effects are broad |
| Endogenous | Made by the body itself. | Adenosine is endogenous |
| Exogenous | Comes from outside the body. | Caffeine is exogenous (you drink it) |
| Neuromodulator | A signaling molecule that tunes the activity of neurons up or down. | Adenosine acts as a neuromodulator that builds up and promotes sleep |

## What Happens After a Receptor Fires

These words show up in the [adenosine receptors](adenosine-receptors.md) cheat sheet (the "Coupling" column) and feed into [signaling to transcription](signaling-to-transcription.md). For the fuller decoder, see [G-protein coupling](g-protein-coupling.md).

| Term | Plain meaning |
|---|---|
| [GPCR](gpcr.md) | G-protein-coupled receptor — a huge family of receptors (adenosine receptors are members) that relay their signal through an attached [G-protein](g-protein.md) |
| [G-protein](g-protein.md) | A molecular switch inside the cell that the receptor flips when activated; held off by GDP and on by GTP. Different flavors push the signal in different directions — see [G-protein coupling](g-protein-coupling.md) |
| GTP / GDP | The two small molecules that set the [G-protein](g-protein.md) switch. **GTP** (guanosine *tri*phosphate) = the "on" fuel; **GDP** (guanosine *di*phosphate) = the spent "off" form with one fewer phosphate. Swapping GDP→GTP turns the switch on |
| GTPase | An enzyme activity that chops the third phosphate off GTP, turning it back into GDP. The G-protein has this built in — it is the self-timer that switches itself off |
| Hydrolysis | Breaking a chemical bond using water. The GTPase "hydrolyzes" GTP to GDP — that is the off-switch step |
| GEF | Guanine-nucleotide Exchange Factor — anything that helps a [G-protein](g-protein.md) drop its GDP so GTP can replace it. An activated [GPCR](gpcr.md) *is* a GEF; see [G-protein switching mechanics](g-protein-switching.md) |
| Gs / Golf | G-protein flavors that **raise** cAMP. Golf is the striatum/smell-specialized version coupled to A2A; see [G-protein coupling](g-protein-coupling.md) |
| Gi/o | A G-protein flavor that **lowers** cAMP; see [G-protein coupling](g-protein-coupling.md) |
| Adenylyl cyclase | The enzyme that produces cAMP; G-proteins turn it up or down |
| cAMP | A small "second messenger" molecule that carries the signal onward inside the cell. More cAMP or less cAMP is the lever many of these receptors pull; see [cAMP signaling and ADORA cascades](camp-signaling.md) |
| Second messenger | Any internal molecule (like cAMP or calcium) that relays a receptor's signal deeper into the cell |

Because adenosine normally *lowers* activity through Gi/o, and caffeine *blocks* adenosine, the net effect of caffeine is **disinhibition** — it removes a brake rather than pressing an accelerator.

## From Signal to Gene

These connect the pharmacology above to the epigenomics the project actually measures.

| Term | Plain meaning |
|---|---|
| Gene | A stretch of DNA that encodes instructions for one product, usually a protein |
| Gene expression | How much a gene is being "read out" and turned into product in a given cell |
| Transcription | The step of copying a gene's DNA into RNA — the first move in expressing it |
| Transcription factor (TF) | A protein that binds DNA and turns nearby genes up or down. CREB, NFAT, and MEF2 in this project are TFs |
| [Kinase](kinase.md) | An enzyme that attaches a phosphate group onto another protein. PKA (protein kinase A) is the main one downstream of cAMP here — full detail on [kinase](kinase.md) |
| Phosphatase | The opposite of a kinase: it *removes* a phosphate group, resetting the protein. Kinase writes the tag, phosphatase erases it |
| Phosphorylation | Adding a phosphate group to a protein — a common on/off tag that changes its activity. PKA phosphorylating CREB is what lets CREB switch genes on |
| Upstream / downstream | Earlier vs later in a chain of events. Receptor binding is upstream; changes in gene expression are downstream |

## Worked Example: One Sentence, Fully Decoded

The thesis sentence in [caffeine molecular targets](caffeine-molecular-targets.md) is: *"caffeine is best modeled as a non-selective antagonist of adenosine receptors."* Translated:

- **antagonist** — caffeine binds the receptor but does **not** turn it on; it blocks it.
- **of adenosine receptors** — the receptors normally used by adenosine, the body's "you're getting sleepy" signal.
- **non-selective** — caffeine blocks **several** receptor subtypes (A1, A2A, and others), not just one, so its effects are broad.
- **best modeled as** — this is the dominant, most explanatory mechanism **at the doses people actually consume**; other targets (PDE, ryanodine receptors, mTOR) only kick in at much higher concentrations. See [signaling to transcription](signaling-to-transcription.md) for those.

Put together: caffeine works by *blocking* the receptors adenosine uses to make you sleepy, across *several subtypes*, and that mechanism *best explains* its real-world effects.

Related pages: [G-protein coupling](g-protein-coupling.md), [cAMP signaling](camp-signaling.md), [caffeine molecular targets](caffeine-molecular-targets.md), [adenosine receptors](adenosine-receptors.md), [signaling to transcription](signaling-to-transcription.md), [epigenomics vocabulary](epigenomics-vocabulary.md)
