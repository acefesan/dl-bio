#!/usr/bin/env python3
"""Interactive UMAP explorer for 8 mammalian root HOGs.

Generates an HTML file with plotly scatter plots — hover over any point
to see omaid, species, scientific name, HOG, root HOG, sub-HOG depth,
and sequence length. Tabs switch between ESM2 model sizes.

Usage:
    python hog_umap_interactive.py
    # Opens: assets/proteins/mammalia/hog_study/umap_interactive.html
"""

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import umap
from plotly.subplots import make_subplots

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "assets" / "proteins" / "mammalia" / "hog_study"

HOG_COLORS = {
    "Rhodopsin": "#E69F00",
    "Cone opsin LW": "#D55E00",
    "Crystallin gamma": "#56B4E9",
    "Sodium channel": "#009E73",
    "Synaptotagmin": "#CC79A7",
    "Keratin type I": "#F0E442",
    "Myosin heavy chain": "#0072B2",
    "Actin": "#999999",
}

SHAPE_TAXA = {
    "HUMAN": "circle",
    "MOUSE": "square",
    "BOVIN": "triangle-up",
    "CANLF": "diamond",
    "ORNAN": "cross",
}

ESM2_MODELS = [
    "esm2_t6_8M_UR50D",
    "esm2_t30_150M_UR50D",
    "esm2_t33_650M_UR50D",
    "esm2_t36_3B_UR50D",
]

MODEL_LABELS = {
    "esm2_t6_8M_UR50D": "ESM2 8M",
    "esm2_t30_150M_UR50D": "ESM2 150M",
    "esm2_t33_650M_UR50D": "ESM2 650M",
    "esm2_t36_3B_UR50D": "ESM2 3B",
}


def compute_umap(embeddings: np.ndarray) -> np.ndarray:
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine",
                        random_state=42)
    return reducer.fit_transform(embeddings)


def build_figure(proteins: pd.DataFrame) -> go.Figure:
    # Derive sub-HOG depth
    proteins = proteins.copy()
    proteins["sub_hog_depth"] = proteins["oma_hog_id"].str.count(r"\.")
    proteins["marker_symbol"] = proteins["species_code"].map(SHAPE_TAXA).fillna("circle-open")
    proteins["is_featured"] = proteins["species_code"].isin(SHAPE_TAXA)

    # Compute UMAP for each model
    umap_coords = {}
    for model in ESM2_MODELS:
        emb_path = OUTPUT_DIR / f"embeddings_{model}.npy"
        if not emb_path.exists():
            print(f"Skipping {model}: no embeddings")
            continue
        print(f"Computing UMAP for {model}...")
        embeddings = np.load(emb_path)
        umap_coords[model] = compute_umap(embeddings)

    if not umap_coords:
        raise RuntimeError("No embedding files found")

    # Build plotly figure with dropdown to switch models
    fig = go.Figure()

    hover_template = (
        "<b>%{customdata[0]}</b><br>"
        "Species: %{customdata[1]} (%{customdata[2]})<br>"
        "Root HOG: %{customdata[3]} (%{customdata[4]})<br>"
        "Full HOG: %{customdata[5]}<br>"
        "Sub-HOG depth: %{customdata[6]}<br>"
        "Seq length: %{customdata[7]} AA"
        "<extra></extra>"
    )

    customdata = proteins[["omaid", "species_code", "scientific_name",
                           "roothog_id", "hog_name", "oma_hog_id",
                           "sub_hog_depth", "sequence_length"]].values

    first_model = True
    for model_idx, (model, coords) in enumerate(umap_coords.items()):
        # One trace per HOG (for legend coloring)
        for hog_name, color in HOG_COLORS.items():
            hog_mask = (proteins["hog_name"] == hog_name).values

            # Featured taxa (solid markers)
            for species_code, symbol in SHAPE_TAXA.items():
                species_mask = (proteins["species_code"] == species_code).values
                mask = hog_mask & species_mask
                if not mask.any():
                    continue
                fig.add_trace(go.Scatter(
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode="markers",
                    marker=dict(
                        color=color,
                        symbol=symbol,
                        size=10,
                        line=dict(width=1, color="black"),
                    ),
                    customdata=customdata[mask],
                    hovertemplate=hover_template,
                    legendgroup=hog_name,
                    showlegend=False,
                    visible=first_model,
                    name=f"{hog_name} ({species_code})",
                ))

            # Other taxa (small open circles)
            other_mask = hog_mask & ~proteins["is_featured"].values
            if other_mask.any():
                fig.add_trace(go.Scatter(
                    x=coords[other_mask, 0],
                    y=coords[other_mask, 1],
                    mode="markers",
                    marker=dict(
                        color=color,
                        symbol="circle",
                        size=6,
                        opacity=0.5,
                        line=dict(width=0.5, color=color),
                    ),
                    customdata=customdata[other_mask],
                    hovertemplate=hover_template,
                    legendgroup=hog_name,
                    showlegend=(first_model and True),
                    visible=first_model,
                    name=hog_name,
                ))

        first_model = False

    # Count traces per model to build visibility toggles
    traces_per_model = len(fig.data) // len(umap_coords)
    buttons = []
    for i, model in enumerate(umap_coords):
        visibility = [False] * len(fig.data)
        for j in range(i * traces_per_model, (i + 1) * traces_per_model):
            visibility[j] = True
        buttons.append(dict(
            label=MODEL_LABELS.get(model, model),
            method="update",
            args=[{"visible": visibility}],
        ))

    # Shape legend annotation
    shape_legend_text = "  ".join(
        f"{'●' if s == 'circle' else '■' if s == 'square' else '▲' if s == 'triangle-up' else '◆' if s == 'diamond' else '✚'} {code}"
        for code, s in SHAPE_TAXA.items()
    )

    fig.update_layout(
        title=dict(
            text="Interactive UMAP — 8 Mammalian Root HOGs across ESM2 scales",
            font=dict(size=16),
        ),
        updatemenus=[dict(
            type="buttons",
            direction="right",
            x=0.5,
            xanchor="center",
            y=1.12,
            showactive=True,
            buttons=buttons,
            font=dict(size=12),
        )],
        annotations=[dict(
            text=f"Shapes: {shape_legend_text} &nbsp; (other species = small dots)",
            xref="paper", yref="paper",
            x=0.5, y=-0.06,
            showarrow=False,
            font=dict(size=11, color="gray"),
            xanchor="center",
        )],
        xaxis=dict(title="UMAP 1", showticklabels=False, zeroline=False),
        yaxis=dict(title="UMAP 2", showticklabels=False, zeroline=False),
        legend=dict(
            title="Root HOG",
            font=dict(size=11),
            itemsizing="constant",
        ),
        width=1100,
        height=750,
        template="plotly_white",
        hovermode="closest",
    )

    return fig


def main():
    proteins = pd.read_feather(OUTPUT_DIR / "hog_proteins.feather")
    print(f"Loaded {len(proteins):,} proteins")

    fig = build_figure(proteins)

    out_path = OUTPUT_DIR / "umap_interactive.html"
    fig.write_html(out_path, include_plotlyjs=True)
    print(f"\nSaved to {out_path}")
    print(f"File size: {out_path.stat().st_size / 1024**2:.1f} MB")


if __name__ == "__main__":
    main()
