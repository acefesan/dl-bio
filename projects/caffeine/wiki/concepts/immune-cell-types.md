# Immune Cell Types

## Summary

Immune cells patrol tissues, detect damage or infection, coordinate inflammation, and remember past threats. In Lab 001, immune cell labels matter because ADORA receptors are not only brain receptors; they also appear in macrophages, monocytes, dendritic cells, mast cells, basophils, T cells, and other immune populations.

This page is the plain-language glossary for immune labels that show up in ADORA dotplots.

## Two Big Branches

Immune cells are often grouped into two broad developmental branches:

| Branch | Plain meaning | Examples |
|---|---|---|
| Myeloid | Fast-acting innate immune and tissue-cleanup cells | monocytes, macrophages, dendritic cells, neutrophils, mast cells, basophils |
| Lymphoid | T/B/NK immune-recognition cells | T cells, B cells, plasma cells, natural killer cells |

This split is a simplification, but it is useful when reading single-cell labels.

## Myeloid Cells

**Myeloid cells** are immune cells from the myeloid lineage. They often act as first responders, tissue sentries, phagocytes, inflammatory signalers, and antigen-presenting cells.

Plain version:

```text
myeloid cell = innate immune / tissue-response cell family
```

In ADORA work, myeloid cells are important because adenosine signaling is deeply tied to inflammation, hypoxia, tissue injury, and immune suppression.

## Monocytes

**Monocytes** circulate in blood and can enter tissues during inflammation. Once in tissues, they can become macrophage-like or dendritic-cell-like states depending on local signals.

Plain version:

```text
monocyte = mobile blood immune cell that can enter tissue and become a local responder
```

In the Tabula Sapiens tongue breakdown, monocytes show some ADORA2B and ADORA3 signal, but the number of tongue monocytes is much smaller than the epithelial populations.

## Macrophages

**Macrophages** are tissue-resident or recruited immune cells that engulf material, clear debris, produce inflammatory signals, remodel tissue, and present antigen.

Plain version:

```text
macrophage = tissue cleanup / alarm / repair immune cell
```

Macrophages frequently show ADORA3 signal in the current Tabula Sapiens ADORA dotplots. That fits the broad idea that adenosine receptors participate in immune regulation and inflammatory tissue states.

## Dendritic Cells

**Dendritic cells** are professional antigen-presenting cells. They sample material from tissues, process it into small pieces, and show those pieces to T cells using MHC molecules.

Plain version:

```text
dendritic cell = immune sentry that collects evidence and presents it to T cells
```

"Dendritic" means branch-like. The name comes from their tree-like projections.

## Myeloid Dendritic Cells

**Myeloid dendritic cells** are dendritic cells from the myeloid branch. They are especially important for detecting tissue danger, sampling microbes or damaged-cell material, and activating or steering T-cell responses.

Plain version:

```text
myeloid dendritic cell = myeloid immune sentry that presents antigen to T cells
```

In Lab 001, a myeloid dendritic cell dot can look striking because mean expression can be bright in a small group. Always check `n_cells`: in the Tongue subset there were only 22 myeloid dendritic cells, so the ADORA3-looking signal is a "look closer" clue, not a stable conclusion by itself.

## T Cells

**T cells** are lymphoid immune cells that recognize antigen and coordinate or execute immune responses.

Plain version:

```text
T cell = adaptive immune cell that recognizes specific molecular evidence
```

Common labels:

| Label | Plain meaning |
|---|---|
| CD4-positive T cell | helper/coordinator T-cell family |
| CD8-positive T cell | cytotoxic/killer T-cell family |
| regulatory T cell | suppressive T-cell family that restrains immune activation |
| thymocyte | developing T cell in the thymus |

ADORA2A is often discussed in T-cell immune regulation, but in the current Tabula Sapiens file the ADORA2A signal is sparse.

## B Cells And Plasma Cells

**B cells** are lymphoid cells that can become antibody-producing cells. **Plasma cells** are the antibody-secreting form.

Plain version:

```text
B cell = antibody-lineage immune cell
plasma cell = antibody factory
```

## Granulocytes, Mast Cells, And Basophils

These are innate immune cells with granules full of signaling molecules.

| Cell type | Plain meaning |
---|---|
| neutrophil | fast inflammatory responder, often antibacterial |
| basophil | rare granulocyte involved in allergy/parasite-type responses |
| mast cell | tissue-resident allergy/inflammation sentinel |

ADORA3 is often visible in mast-cell, basophil, macrophage, and related myeloid contexts in Lab 001 outputs.

## Reading Immune Dots In Lab 001

When an immune-cell dot is bright:

1. Check `n_cells`; rare immune subsets can look bright from small counts.
2. Check tissue context; macrophages in lung, spleen, tongue, and fat may not mean the same state.
3. Separate prevalence from intensity; a small bright dot can mean few high cells.
4. Compare related labels: monocyte, macrophage, tissue-resident macrophage, dendritic cell, and myeloid leukocyte can overlap biologically but are not identical.
5. Treat ADORA mRNA as receptor transcript evidence, not proof of immune response to caffeine.

Related pages: [001 ADORA interpretation](../labs/001-adora-interpretation.md), [adenosine receptors](adenosine-receptors.md), [single-cell RNA-seq measurement](single-cell-rna-seq-measurement.md), [cell type response model](cell-type-response-model.md), [unexpected responsive cell types](unexpected-responsive-cell-types.md)
