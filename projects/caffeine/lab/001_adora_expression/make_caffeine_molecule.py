#!/usr/bin/env python
"""Publication-quality 2D skeletal structure of caffeine, drawn with RDKit.

Caffeine is 1,3,7-trimethylxanthine (C8H10N4O2): a fused purine — a
six-membered pyrimidinedione ring (the two C=O groups) fused to a
five-membered imidazole ring — carrying three N-methyl groups.

RDKit parses the canonical SMILES, generates clean 2D coordinates, and
renders the depiction. This is the standard cheminformatics way to draw a
molecule from text, and it gets the fused-ring geometry, both carbonyls,
and the three N-methyls correct in a way ASCII cannot.

Reference identifiers (what a chemist would store):
    SMILES : CN1C=NC2=C1C(=O)N(C)C(=O)N2C
    InChI  : InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3

Output: figures/caffeine_molecule.{png,svg}
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

LAB_DIR = Path(__file__).resolve().parent
FIG_DIR = LAB_DIR / "figures"

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"


def build_mol() -> Chem.Mol:
    mol = Chem.MolFromSmiles(CAFFEINE_SMILES)
    if mol is None:
        raise ValueError("RDKit failed to parse the caffeine SMILES")
    # Sanity checks against the known formula (C8H10N4O2, MW ~194.19).
    formula = rdMolDescriptors.CalcMolFormula(mol)
    assert formula == "C8H10N4O2", f"unexpected formula: {formula}"
    rdDepictor.Compute2DCoords(mol)
    return mol


def render(mol: Chem.Mol, fmt: str, px: int = 900) -> None:
    if fmt == "svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(px, px)
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(px, px)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = False
    opts.bondLineWidth = 2
    opts.padding = 0.12
    opts.legendFontSize = 20
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol, legend="Caffeine  (1,3,7-trimethylxanthine, C8H10N4O2)"
    )
    drawer.FinishDrawing()
    out = FIG_DIR / f"caffeine_molecule.{fmt}"
    data = drawer.GetDrawingText()
    mode = "w" if fmt == "svg" else "wb"
    with open(out, mode) as fh:
        fh.write(data)
    print(f"  -> {out}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    mol = build_mol()
    for fmt in ("png", "svg"):
        render(mol, fmt)


if __name__ == "__main__":
    main()
