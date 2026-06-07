#!/usr/bin/env python
"""Summarize Tabula Sapiens donor IDs in the cached H5AD."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import h5py
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
OUT = LAB_DIR / "cache" / "tabula_sapiens_donor_summary.csv"

DONOR_COLUMNS = [
    "donor_id",
    "sex",
    "development_stage",
    "self_reported_ethnicity",
]


def decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_categorical(f: h5py.File, column: str) -> list[str]:
    node = f["obs"][column]
    categories = decode(node["categories"][:])
    codes = node["codes"][:]
    return ["missing" if code < 0 else categories[code] for code in codes]


def mode(values: list[str]) -> str:
    return Counter(values).most_common(1)[0][0]


def main() -> None:
    with h5py.File(H5AD, "r") as f:
        values = {column: read_categorical(f, column) for column in DONOR_COLUMNS}
        tissues = read_categorical(f, "tissue_in_publication")

    by_donor: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for i, donor in enumerate(values["donor_id"]):
        for column in DONOR_COLUMNS[1:]:
            by_donor[donor][column].append(values[column][i])
        by_donor[donor]["tissue_in_publication"].append(tissues[i])

    rows = []
    for donor, cols in by_donor.items():
        tissue_counts = Counter(cols["tissue_in_publication"])
        rows.append(
            {
                "donor_id": donor,
                "n_cells": sum(tissue_counts.values()),
                "sex": mode(cols["sex"]),
                "development_stage": mode(cols["development_stage"]),
                "self_reported_ethnicity": mode(cols["self_reported_ethnicity"]),
                "n_broad_tissues": len(tissue_counts),
                "broad_tissues": "; ".join(sorted(tissue_counts)),
            }
        )

    df = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"Wrote {OUT}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
