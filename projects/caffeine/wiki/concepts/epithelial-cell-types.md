# Epithelial Cell Types

## Summary

**Epithelial cells** are the cells that cover surfaces and line cavities: skin, airway, gut, glands, ducts, bladder, mouth, and tongue. They form barriers between the body and the outside world, absorb and secrete material, sense local conditions, and repair tissue after injury.

This page exists because Lab 001's Tabula Sapiens ADORA dotplots show a recurring [ADORA2B](adenosine-receptors.md) signal in epithelial and epithelial-like populations, especially the `Tongue` tissue group.

## Epithelial Cells

An epithelial cell is a surface-lining or gland-forming cell.

Plain version:

```text
epithelial cell = a cell that makes a sheet, lining, duct, gland, or barrier
```

Examples:

| Place | Epithelial role |
|---|---|
| Mouth and tongue | Protective surface layer against food, abrasion, microbes |
| Intestine | Barrier plus absorption and secretion |
| Airway | Barrier, mucus handling, cilia, immune interface |
| Bladder | Stretchy urine-facing barrier |
| Salivary/pancreatic ducts | Fluid and enzyme transport surfaces |

Epithelia are not passive wallpaper. They sense stress, injury, oxygen, microbes, and inflammatory signals. That is why receptor and signaling genes can show strong epithelial patterns.

## Basal Cells

**Basal cells** are epithelial cells sitting near the bottom, or base, of an epithelial layer. They often act as progenitor-like maintenance cells that replenish the surface.

Plain version:

```text
basal cell = the lower-layer epithelial cell that helps maintain and replace the tissue surface
```

Why they matter:

- They are close to the tissue scaffold under the epithelium.
- They can divide and produce more differentiated surface cells.
- They are common in stratified epithelia such as skin, airway, oral mucosa, and tongue.
- They often carry wound-repair, barrier, and stress-response programs.

In the Tabula Sapiens tongue subset, basal cells are the biggest ADORA2B contributor:

| Tongue cell type | Cells | ADORA2B percent expressing | Mean ADORA2B |
|---|---:|---:|---:|
| basal cell | 16,263 | 35.5% | 0.216 |

This means the tissue-level `Tongue` ADORA2B dot is not mainly "taste receptor cells express ADORA2B"; it is mostly a basal epithelial signal.

## Stratified Squamous Epithelial Cells

**Stratified** means many layers. **Squamous** means flat, scale-like cells. A stratified squamous epithelium is a layered protective surface built to tolerate friction and abrasion.

Plain version:

```text
stratified squamous epithelial cell = a flat protective surface cell in a many-layered epithelial sheet
```

Common places:

| Tissue | Why this form helps |
|---|---|
| Tongue and oral mucosa | Handles chewing, swallowing, heat, acid, microbes |
| Esophagus | Handles swallowed food abrasion |
| Skin | Handles external physical stress |
| Cervix/vagina | Barrier and mechanical protection |

In the tongue subset, stratified squamous epithelial cells are the second major ADORA2B contributor:

| Tongue cell type | Cells | ADORA2B percent expressing | Mean ADORA2B |
|---|---:|---:|---:|
| stratified squamous epithelial cell | 12,323 | 14.2% | 0.065 |

So the tongue pattern is best read as:

```text
ADORA2B is visible in tongue barrier epithelium,
especially basal cells and stratified squamous epithelial cells.
```

## Why ADORA2B Might Show Up In Barrier Epithelia

ADORA2B is often discussed as a lower-affinity adenosine receptor that becomes more relevant when extracellular adenosine is high, such as stress, hypoxia, inflammation, or tissue injury contexts. Barrier epithelia are exactly the kind of cells that constantly experience local stress and repair signals.

For this project, the conservative interpretation is:

```text
Tongue ADORA2B marks epithelial/barrier contexts worth following up.
```

It does **not** yet mean:

- caffeine strongly affects tongue epithelium in vivo,
- receptor protein is abundant on the cell surface,
- taste perception is mediated through this signal,
- ADORA2B is unique to tongue.

The dotplot already shows ADORA2B in other epithelial-rich tissues and cell types too.

## Tongue Versus Taste Cells

The word "tongue" can tempt a taste-first interpretation. In this atlas, that would be misleading.

The `Tongue` tissue group is dominated by epithelial barrier cells:

| Tongue cell type | Cells |
|---|---:|
| basal cell | 16,263 |
| stratified squamous epithelial cell | 12,323 |
| taste receptor cell | 32 |

Only 32 taste receptor cells are labeled in this subset. The strong tissue-level ADORA2B signal is therefore a **surface epithelium** result, not a robust taste-cell result.

## Practical Reading Rule

When a tissue-level dot is large:

1. Ask which cell types dominate that tissue.
2. Decompose the tissue into `cell_type x gene`.
3. Check whether the signal is coming from a biologically specific minority or a numerically dominant broad compartment.
4. Treat broad terms like `Tongue`, `Lung`, or `Pancreas` as mixtures, not as single biological cell states.

Related pages: [001 ADORA interpretation](../labs/001-adora-interpretation.md), [single-cell RNA-seq measurement](single-cell-rna-seq-measurement.md), [adenosine receptors](adenosine-receptors.md), [unexpected responsive cell types](unexpected-responsive-cell-types.md)
