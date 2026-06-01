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


COMPARE_JS = """
<style>
#compare-panel {
    position: fixed; bottom: 0; left: 0; right: 0;
    max-height: 45vh; overflow-y: auto;
    background: #fff; border-top: 2px solid #333;
    font-family: 'Segoe UI', system-ui, sans-serif; font-size: 13px;
    padding: 10px 16px; z-index: 9999;
    box-shadow: 0 -2px 12px rgba(0,0,0,0.15);
    display: none;
}
#compare-panel h3 { margin: 0 0 8px 0; font-size: 15px; }
#compare-panel .actions { margin-bottom: 8px; }
#compare-panel button {
    padding: 4px 12px; margin-right: 6px; cursor: pointer;
    border: 1px solid #888; border-radius: 4px; background: #f0f0f0;
    font-size: 12px;
}
#compare-panel button:hover { background: #ddd; }
#compare-table {
    width: 100%; border-collapse: collapse; table-layout: fixed;
}
#compare-table th, #compare-table td {
    border: 1px solid #ddd; padding: 4px 8px; text-align: left;
    vertical-align: top; word-break: break-all;
}
#compare-table th { background: #f5f5f5; font-weight: 600; width: 110px; }
.hog-match { color: #888; }
.hog-diverge { color: #d32f2f; font-weight: 700; }
.hog-level { font-family: monospace; }
.click-hint {
    position: fixed; bottom: 12px; left: 50%; transform: translateX(-50%);
    background: rgba(0,0,0,0.75); color: #fff; padding: 8px 18px;
    border-radius: 6px; font-family: system-ui; font-size: 13px; z-index: 9998;
}
</style>

<div class="click-hint" id="click-hint">Click points to compare HOG paths. Click same point again to remove.</div>
<div id="compare-panel">
    <h3>HOG Comparison <span id="compare-count"></span></h3>
    <div class="actions">
        <button onclick="clearCompare()">Clear All</button>
        <button onclick="closePanel()">Close</button>
    </div>
    <table id="compare-table"><tbody id="compare-body"></tbody></table>
</div>

<script>
var selectedPoints = [];

// Listen for click on any plotly trace
var plotDiv = document.getElementsByClassName('plotly-graph-div')[0];
plotDiv.on('plotly_click', function(data) {
    var pt = data.points[0];
    var cd = pt.customdata;
    if (!cd) return;
    var omaid = cd[0];

    // Toggle: remove if already selected
    var idx = selectedPoints.findIndex(p => p.omaid === omaid);
    if (idx >= 0) {
        selectedPoints.splice(idx, 1);
    } else {
        selectedPoints.push({
            omaid: cd[0],
            species: cd[1],
            sciname: cd[2],
            roothog: cd[3],
            hogname: cd[4],
            fullhog: cd[5],
            depth: cd[6],
            seqlen: cd[7],
        });
    }
    renderCompare();
});

function renderCompare() {
    var panel = document.getElementById('compare-panel');
    var hint = document.getElementById('click-hint');
    if (selectedPoints.length === 0) {
        panel.style.display = 'none';
        hint.style.display = 'block';
        return;
    }
    hint.style.display = 'none';
    panel.style.display = 'block';
    document.getElementById('compare-count').textContent =
        '(' + selectedPoints.length + ' proteins)';

    // Split HOG paths into levels
    var maxDepth = 0;
    var hogLevels = selectedPoints.map(function(p) {
        var parts = p.fullhog ? p.fullhog.split('.') : [''];
        if (parts.length > maxDepth) maxDepth = parts.length;
        return parts;
    });

    // Find first diverging level between all pairs
    var divergeLevel = maxDepth;
    if (selectedPoints.length >= 2) {
        for (var lv = 0; lv < maxDepth; lv++) {
            var vals = new Set(hogLevels.map(h => h[lv] || ''));
            if (vals.size > 1) { divergeLevel = lv; break; }
        }
    }

    // Build table
    var rows = [
        {label: 'OMA ID', key: 'omaid'},
        {label: 'Species', fn: function(p) { return p.species + ' (' + p.sciname + ')'; }},
        {label: 'Root HOG', fn: function(p) { return p.hogname + ' (' + p.roothog + ')'; }},
        {label: 'Full HOG', key: 'fullhog'},
        {label: 'Sub-HOG depth', key: 'depth'},
        {label: 'Seq length', fn: function(p) { return p.seqlen + ' AA'; }},
    ];

    // Add HOG level rows
    for (var lv = 0; lv < maxDepth; lv++) {
        (function(level) {
            rows.push({
                label: level === 0 ? 'Root' : 'Level ' + level,
                fn: function(p, i) {
                    var val = hogLevels[i][level] || '—';
                    var cls = level < divergeLevel ? 'hog-match' :
                              level === divergeLevel ? 'hog-diverge' : 'hog-level';
                    return '<span class="hog-level ' + cls + '">' + val + '</span>';
                },
                isHtml: true
            });
        })(lv);
    }

    var html = '';
    rows.forEach(function(row) {
        html += '<tr><th>' + row.label + '</th>';
        selectedPoints.forEach(function(p, i) {
            var val = row.fn ? row.fn(p, i) : p[row.key];
            if (row.isHtml) {
                html += '<td>' + val + '</td>';
            } else {
                html += '<td>' + String(val) + '</td>';
            }
        });
        html += '</tr>';
    });
    document.getElementById('compare-body').innerHTML = html;
}

function clearCompare() {
    selectedPoints = [];
    renderCompare();
}

function closePanel() {
    document.getElementById('compare-panel').style.display = 'none';
    document.getElementById('click-hint').style.display = 'block';
}
</script>
"""


def main():
    proteins = pd.read_feather(OUTPUT_DIR / "hog_proteins.feather")
    print(f"Loaded {len(proteins):,} proteins")

    fig = build_figure(proteins)

    out_path = OUTPUT_DIR / "umap_interactive.html"
    fig.write_html(out_path, include_plotlyjs=True, post_script="",
                   full_html=True)

    # Inject comparison JS before closing </body>
    html = out_path.read_text()
    html = html.replace("</body>", COMPARE_JS + "\n</body>")
    out_path.write_text(html)

    print(f"\nSaved to {out_path}")
    print(f"File size: {out_path.stat().st_size / 1024**2:.1f} MB")


if __name__ == "__main__":
    main()
