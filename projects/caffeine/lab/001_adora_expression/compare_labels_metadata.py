#!/usr/bin/env python
"""Compare label metadata between Tabula Sapiens and HBCA non-neuronal."""

from __future__ import annotations

import argparse
import json
import re
from difflib import SequenceMatcher
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_TABULA = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_HBCA = LAB_DIR / "cache" / "human_brain_cell_atlas" / "hbca_all_non_neuronal_b165f033.h5ad"
DEFAULT_OUT_DIR = LAB_DIR / "figures"

TABULA_LABEL_FIELDS = ("cell_type", "broad_cell_class", "tissue_in_publication", "tissue")
HBCA_LABEL_FIELDS = (
    "cell_type",
    "supercluster_term",
    "ROIGroup",
    "ROIGroupCoarse",
    "ROIGroupFine",
    "roi",
    "tissue",
)
STOPWORDS = {
    "cell",
    "cells",
    "human",
    "primary",
    "of",
    "the",
}


def decode_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_obs_values(h5ad_path: Path, obs_column: str) -> np.ndarray:
    with h5py.File(h5ad_path, "r") as f:
        node = f["obs"][obs_column]
        if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
            codes = node["codes"][:]
            categories = np.array(decode_values(node["categories"][:]), dtype=object)
            values = np.full(codes.shape, None, dtype=object)
            valid = codes >= 0
            values[valid] = categories[codes[valid]]
            return values
        if isinstance(node, h5py.Dataset):
            return np.array(decode_values(node[:]), dtype=object)
    raise ValueError(f"Unsupported obs column encoding: {h5ad_path}: obs/{obs_column}")


def obs_overview(h5ad_path: Path, dataset: str, likely_fields: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        n_obs = len(obs[obs.attrs["_index"]])
        for column in obs.keys():
            node = obs[column]
            encoding = node.attrs.get("encoding-type", "dataset") if hasattr(node, "attrs") else "dataset"
            n_categories = None
            if isinstance(node, h5py.Group) and "categories" in node:
                n_categories = int(len(node["categories"]))
            rows.append(
                {
                    "dataset": dataset,
                    "obs_column": column,
                    "is_requested_label_field": column in likely_fields,
                    "encoding_type": str(encoding),
                    "n_obs": n_obs,
                    "n_categories": n_categories,
                }
            )
    return pd.DataFrame(rows).sort_values(["is_requested_label_field", "obs_column"], ascending=[False, True])


def category_counts(h5ad_path: Path, dataset: str, fields: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for field in fields:
        values = read_obs_values(h5ad_path, field)
        counts = pd.Series(values, dtype="object").value_counts(dropna=False)
        total = int(counts.sum())
        for label, n_cells in counts.items():
            rows.append(
                {
                    "dataset": dataset,
                    "field": field,
                    "label": "" if label is None else str(label),
                    "n_cells": int(n_cells),
                    "pct_cells": float(n_cells / total * 100.0),
                }
            )
    return pd.DataFrame(rows).sort_values(["dataset", "field", "n_cells", "label"], ascending=[True, True, False, True])


def normalize_label(label: str) -> str:
    label = label.lower()
    label = label.replace("+", " positive ")
    label = re.sub(r"[^a-z0-9]+", " ", label)
    tokens = []
    for token in label.split():
        if token.endswith("s") and not token.endswith(("ss", "us", "ous")) and len(token) > 4:
            token = token[:-1]
        if token not in STOPWORDS:
            tokens.append(token)
    return " ".join(tokens)


def token_jaccard(left: str, right: str) -> float:
    left_tokens = set(left.split())
    right_tokens = set(right.split())
    if not left_tokens and not right_tokens:
        return 1.0
    if not left_tokens or not right_tokens:
        return 0.0
    return len(left_tokens & right_tokens) / len(left_tokens | right_tokens)


def compare_cell_type_overlap(counts: pd.DataFrame) -> pd.DataFrame:
    tabula = counts[(counts["dataset"] == "tabula_sapiens") & (counts["field"] == "cell_type")].copy()
    hbca = counts[(counts["dataset"] == "hbca_non_neuronal") & (counts["field"] == "cell_type")].copy()
    tabula["normalized_label"] = tabula["label"].map(normalize_label)
    hbca["normalized_label"] = hbca["label"].map(normalize_label)

    rows = []
    tabula_records = tabula[["label", "normalized_label", "n_cells", "pct_cells"]].to_dict("records")
    for hbca_row in hbca.to_dict("records"):
        best_score = None
        best_payload = None
        for tabula_row in tabula_records:
            exact = hbca_row["normalized_label"] == tabula_row["normalized_label"]
            seq = SequenceMatcher(None, hbca_row["normalized_label"], tabula_row["normalized_label"]).ratio()
            jac = token_jaccard(hbca_row["normalized_label"], tabula_row["normalized_label"])
            score = max(seq, jac, 1.0 if exact else 0.0)
            candidate_score = (score, exact, seq, jac, int(tabula_row["n_cells"]))
            if best_score is None or candidate_score > best_score:
                best_score = candidate_score
                best_payload = (score, exact, seq, jac, tabula_row)
        assert best_payload is not None
        score, exact, seq, jac, tabula_row = best_payload
        rows.append(
            {
                "hbca_cell_type": hbca_row["label"],
                "hbca_normalized": hbca_row["normalized_label"],
                "hbca_n_cells": hbca_row["n_cells"],
                "hbca_pct_cells": hbca_row["pct_cells"],
                "tabula_best_cell_type": tabula_row["label"],
                "tabula_best_normalized": tabula_row["normalized_label"],
                "tabula_best_n_cells": tabula_row["n_cells"],
                "tabula_best_pct_cells": tabula_row["pct_cells"],
                "normalized_exact_match": bool(exact),
                "match_quality": "exact" if exact else "similar" if score >= 0.75 or jac >= 0.5 else "weak_or_no_direct_match",
                "similarity_score": float(score),
                "sequence_similarity": float(seq),
                "token_jaccard": float(jac),
            }
        )

    overlap = pd.DataFrame(rows).sort_values(
        ["normalized_exact_match", "similarity_score", "hbca_n_cells"],
        ascending=[False, False, False],
    )
    return overlap


def write_summary(
    out_dir: Path,
    overview: pd.DataFrame,
    counts: pd.DataFrame,
    overlap: pd.DataFrame,
) -> None:
    tabula_cell_types = int(
        overview[(overview["dataset"] == "tabula_sapiens") & (overview["obs_column"] == "cell_type")]["n_categories"].iloc[0]
    )
    hbca_cell_types = int(
        overview[(overview["dataset"] == "hbca_non_neuronal") & (overview["obs_column"] == "cell_type")]["n_categories"].iloc[0]
    )
    exact_matches = int(overlap["normalized_exact_match"].sum())
    similar_matches = int((overlap["match_quality"] == "similar").sum())

    recommendations = [
        {
            "comparison": "cell identity",
            "tabula_label": "cell_type",
            "hbca_label": "cell_type or supercluster_term",
            "reason": (
                "Tabula has fine whole-body cell types, while HBCA non-neuronal has 11 broad cell_type labels "
                "and 10 supercluster_term labels. For ADORA summaries, use HBCA cell_type/supercluster_term "
                "inside the brain dataset and compare cautiously to Tabula cell_type."
            ),
        },
        {
            "comparison": "broad class",
            "tabula_label": "broad_cell_class",
            "hbca_label": "supercluster_term",
            "reason": (
                "These are the closest broad biological groupings across datasets, useful for comparing "
                "immune, endothelial, epithelial-like, and glial/non-neuronal compartments."
            ),
        },
        {
            "comparison": "anatomical context",
            "tabula_label": "tissue_in_publication or tissue",
            "hbca_label": "ROIGroupCoarse, ROIGroupFine, roi, tissue",
            "reason": (
                "Tabula tissue labels describe organs/tissues; HBCA labels describe brain regions. "
                "They should be presented side by side rather than normalized into one shared tissue axis."
            ),
        },
    ]
    summary = {
        "inputs": {
            "tabula_sapiens": str(DEFAULT_TABULA),
            "hbca_non_neuronal": str(DEFAULT_HBCA),
        },
        "n_label_fields": {
            "tabula_requested_fields": len(TABULA_LABEL_FIELDS),
            "hbca_requested_fields": len(HBCA_LABEL_FIELDS),
        },
        "cell_type_category_counts": {
            "tabula_sapiens": tabula_cell_types,
            "hbca_non_neuronal": hbca_cell_types,
            "normalized_exact_hbca_to_tabula_matches": exact_matches,
            "similar_non_exact_hbca_to_tabula_matches": similar_matches,
        },
        "recommendations": recommendations,
        "notes": [
            "Existing UMAP coordinates are dataset-specific and should not be merged as if they share axes.",
            "For ADORA pattern comparison, dotplots/pseudobulk tables by label are more defensible than comparing raw UMAP geometry.",
            "HBCA non-neuronal is brain-specific; Tabula Sapiens is whole-body and lacks brain tissue in this cache.",
        ],
    }

    json_path = out_dir / "compare_labels_summary.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    md_lines = [
        "# Label comparison: Tabula Sapiens vs HBCA non-neuronal",
        "",
        f"- Tabula `cell_type` categories: {tabula_cell_types}",
        f"- HBCA non-neuronal `cell_type` categories: {hbca_cell_types}",
        f"- Exact normalized HBCA-to-Tabula cell-type matches: {exact_matches} / {len(overlap)}",
        f"- Similar non-exact HBCA-to-Tabula cell-type matches: {similar_matches} / {len(overlap)}",
        "",
        "## Best labels for ADORA comparisons",
        "",
    ]
    for item in recommendations:
        md_lines.extend(
            [
                f"### {item['comparison'].title()}",
                "",
                f"- Tabula: `{item['tabula_label']}`",
                f"- HBCA: `{item['hbca_label']}`",
                f"- Why: {item['reason']}",
                "",
            ]
        )
    md_lines.extend(["## Notes", ""])
    md_lines.extend([f"- {note}" for note in summary["notes"]])
    md_path = out_dir / "compare_labels_summary.md"
    md_path.write_text("\n".join(md_lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tabula-h5ad", type=Path, default=DEFAULT_TABULA)
    parser.add_argument("--hbca-h5ad", type=Path, default=DEFAULT_HBCA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    overview = pd.concat(
        [
            obs_overview(args.tabula_h5ad, "tabula_sapiens", TABULA_LABEL_FIELDS),
            obs_overview(args.hbca_h5ad, "hbca_non_neuronal", HBCA_LABEL_FIELDS),
        ],
        ignore_index=True,
    )
    counts = pd.concat(
        [
            category_counts(args.tabula_h5ad, "tabula_sapiens", TABULA_LABEL_FIELDS),
            category_counts(args.hbca_h5ad, "hbca_non_neuronal", HBCA_LABEL_FIELDS),
        ],
        ignore_index=True,
    )
    overlap = compare_cell_type_overlap(counts)

    overview_path = args.out_dir / "compare_labels_obs_columns.csv"
    counts_path = args.out_dir / "compare_labels_category_counts.csv"
    overlap_path = args.out_dir / "compare_labels_cell_type_overlap.csv"

    overview.to_csv(overview_path, index=False)
    counts.to_csv(counts_path, index=False)
    overlap.to_csv(overlap_path, index=False)
    write_summary(args.out_dir, overview, counts, overlap)

    for path in [
        overview_path,
        counts_path,
        overlap_path,
        args.out_dir / "compare_labels_summary.json",
        args.out_dir / "compare_labels_summary.md",
    ]:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
