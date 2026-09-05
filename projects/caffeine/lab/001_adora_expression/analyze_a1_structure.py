#!/usr/bin/env python3
"""Find coherent A1-pathway structure hidden by atlas-wide raw-score UMAPs."""
from __future__ import annotations

import html
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.spatial import cKDTree

from visualize_a1_pathway_all_cells import ATLASES, OUT, categorical, stratified_sample


def enrichments(scores: np.ndarray, cell_types: np.ndarray, tissues: np.ndarray) -> pd.DataFrame:
    q95, q99 = np.quantile(scores, [.95, .99])
    rows: list[pd.DataFrame] = []
    frame = pd.DataFrame({"cell_type": cell_types, "region": tissues, "score": scores})
    for dimension in ("cell_type", "region"):
        grouped = frame.groupby(dimension, observed=True, sort=False)["score"]
        table = grouped.agg(cells="size", median="median", mean="mean", p95=lambda x: x.quantile(.95))
        table["top_5pct_fraction"] = grouped.apply(lambda x: float((x >= q95).mean()), include_groups=False)
        table["top_1pct_fraction"] = grouped.apply(lambda x: float((x >= q99).mean()), include_groups=False)
        table["top_5pct_enrichment"] = table.top_5pct_fraction / .05
        table["top_1pct_enrichment"] = table.top_1pct_fraction / .01
        table = table.reset_index(names="group")
        table.insert(0, "dimension", dimension)
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def smooth_map(coords: np.ndarray, scores: np.ndarray, cell_types: np.ndarray, tissues: np.ndarray,
               label: str, output: Path, sample_size: int = 80_000, k: int = 30) -> dict:
    idx = stratified_sample(cell_types, sample_size, seed=41)
    xy = coords[idx].astype(np.float64)
    raw = scores[idx].astype(np.float64)
    distances, neighbors = cKDTree(xy).query(xy, k=k)
    scale = np.maximum(distances[:, -1], 1e-8)
    weights = np.exp(-0.5 * (distances / scale[:, None]) ** 2)
    smooth = (scores[idx[neighbors]] * weights).sum(axis=1) / weights.sum(axis=1)
    local_support = (scores[idx[neighbors]] >= np.quantile(scores, .95)).mean(axis=1)

    custom = np.column_stack([cell_types[idx], tissues[idx], raw, smooth, local_support])
    common = dict(
        x=xy[:, 0], y=xy[:, 1], mode="markers", customdata=custom,
        hovertemplate=("<b>%{customdata[0]}</b><br>region: %{customdata[1]}"
                       "<br>raw score: %{customdata[2]:.4f}<br>kNN score: %{customdata[3]:.4f}"
                       "<br>top-5% neighbors: %{customdata[4]:.0%}<extra></extra>"),
    )
    traces = [
        go.Scattergl(**common, marker={"size": 3, "color": raw, "colorscale": "Magma", "opacity": .48,
                                           "colorbar": {"title": "raw"}}, name="Raw score", visible=False),
        go.Scattergl(**common, marker={"size": 3, "color": smooth, "colorscale": "Magma", "opacity": .72,
                                           "colorbar": {"title": f"kNN mean (k={k})"}}, name="Smoothed score"),
        go.Scattergl(**common, marker={"size": 3, "color": local_support, "colorscale": "Turbo", "opacity": .72,
                                           "cmin": 0, "cmax": 1, "colorbar": {"title": "local support"}},
                     name="Top-tail support", visible=False),
    ]
    fig = go.Figure(traces)
    fig.update_layout(
        title=f"{label}: coherent A1 neighborhoods ({len(idx):,}-cell stratified sample)",
        template="plotly_dark", paper_bgcolor="#08080c", plot_bgcolor="#08080c", dragmode="pan",
        margin={"l": 15, "r": 15, "t": 115, "b": 15},
        xaxis={"visible": False}, yaxis={"visible": False, "scaleanchor": "x", "scaleratio": 1},
        updatemenus=[{"type": "buttons", "direction": "right", "x": 0, "y": 1.09,
                      "buttons": [
                          {"label": "Raw", "method": "update", "args": [{"visible": [True, False, False]}]},
                          {"label": f"Smoothed k={k}", "method": "update", "args": [{"visible": [False, True, False]}]},
                          {"label": "Top-tail support", "method": "update", "args": [{"visible": [False, False, True]}]},
                      ]}],
    )
    fig.write_html(output, include_plotlyjs=True, full_html=True,
                   config={"responsive": True, "scrollZoom": True, "displaylogo": False})
    return {"sample_cells": len(idx), "k": k, "raw_p99": float(np.quantile(raw, .99)),
            "smooth_p99": float(np.quantile(smooth, .99))}


def enrichment_page(table: pd.DataFrame, label: str, output: Path) -> None:
    parts = []
    for dimension in ("cell_type", "region"):
        view = table[(table.dimension == dimension) & (table.cells >= 100)].nlargest(40, "top_5pct_enrichment")
        rows = "".join(
            f"<tr><td>{html.escape(str(r.group))}</td><td>{r.cells:,.0f}</td><td>{r.top_5pct_enrichment:.2f}×</td>"
            f"<td>{r.top_1pct_enrichment:.2f}×</td><td>{r.median:.4f}</td><td>{r.p95:.4f}</td></tr>"
            for r in view.itertuples()
        )
        parts.append(f"<section><h2>{dimension.replace('_', ' ').title()}</h2><div class='scroll'><table>"
                     "<thead><tr><th>Group</th><th>Cells</th><th>Top 5% enrichment</th><th>Top 1% enrichment</th>"
                     f"<th>Median</th><th>95th percentile</th></tr></thead><tbody>{rows}</tbody></table></div></section>")
    output.write_text(f"""<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>{html.escape(label)} A1 enrichment</title><style>
:root{{color-scheme:dark;font-family:system-ui,sans-serif}}body{{margin:0;background:#08080c;color:#f4f4f5;padding:22px}}
main{{max-width:1200px;margin:auto}}p{{color:#b8b8c6;line-height:1.5}}section{{margin:28px 0}}.scroll{{overflow:auto;border:1px solid #30303c;border-radius:14px}}
table{{border-collapse:collapse;width:100%;background:#111119}}th,td{{padding:11px 13px;text-align:right;border-bottom:1px solid #292933;white-space:nowrap}}
th:first-child,td:first-child{{text-align:left;white-space:normal}}th{{position:sticky;top:0;background:#1a1a25;color:#f0b4ff}}
</style></head><body><main><h1>{html.escape(label)}: A1 enrichment</h1>
<p>Groups ranked by overrepresentation among the atlas-wide highest-scoring cells. Values above 1× indicate enrichment. Groups with fewer than 100 cells are excluded here but remain in the downloadable CSV.</p>
{''.join(parts)}</main></body></html>""", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for slug, config in ATLASES.items():
        print(f"Loading {config['label']}", flush=True)
        scores = np.load(OUT / f"{slug}_a1_pathway_scores.npz")["score"]
        with h5py.File(config["path"], "r") as handle:
            coords = handle[f"obsm/{config['umap']}"][:].astype(np.float32)
            cell_types = categorical(handle["obs"], config["cell_type"])
            tissues = categorical(handle["obs"], config["tissue"])
        table = enrichments(scores, cell_types, tissues)
        table.to_csv(OUT / f"{slug}_a1_enrichment.csv", index=False)
        enrichment_page(table, config["label"], OUT / f"{slug}_enrichment.html")
        manifest[slug] = smooth_map(coords, scores, cell_types, tissues, config["label"],
                                    OUT / f"{slug}_neighborhoods.html")
    (OUT / "structure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
