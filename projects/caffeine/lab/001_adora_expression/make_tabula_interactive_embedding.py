#!/usr/bin/env python
"""Build an interactive Tabula Sapiens embedding browser."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import plotly.graph_objects as go


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_ADORA_CACHE = LAB_DIR / "cache" / "tabula_sapiens_adora_expression.npz"
DEFAULT_OUT = LAB_DIR / "figures" / "tabula_sapiens_X_umap_interactive.html"
ADORA_GENES = ("ADORA1", "ADORA2A", "ADORA2B", "ADORA3")

METADATA_COLUMNS = [
    "donor_id",
    "tissue_in_publication",
    "tissue",
    "cell_type",
    "broad_cell_class",
    "assay",
    "sex",
    "development_stage",
    "self_reported_ethnicity",
    "disease",
    "compartment",
]


def decode_array(values: np.ndarray) -> list[str]:
    out: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8"))
        else:
            out.append(str(value))
    return out


def read_obs_column(f: h5py.File, column: str, idx: np.ndarray) -> np.ndarray:
    node = f["obs"][column]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        categories = np.array(decode_array(node["categories"][:]), dtype=object)
        codes = node["codes"][idx]
        values = np.full(codes.shape, "missing", dtype=object)
        ok = codes >= 0
        values[ok] = categories[codes[ok]]
        return values

    values = node[idx]
    if values.dtype.kind in {"S", "O"}:
        return np.array(decode_array(values), dtype=object)
    return values.astype(str).astype(object)


def categorical_codes(f: h5py.File, column: str, idx: np.ndarray) -> tuple[np.ndarray, list[str], np.ndarray]:
    node = f["obs"][column]
    if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
        raise ValueError(f"obs/{column!r} is not a categorical column")
    codes = node["codes"][idx]
    categories = decode_array(node["categories"][:])
    values = np.full(codes.shape, "missing", dtype=object)
    ok = codes >= 0
    values[ok] = np.array(categories, dtype=object)[codes[ok]]
    return codes, categories, values


def choose_indices(n_obs: int, max_points: int | None, seed: int) -> np.ndarray:
    if max_points is None or max_points >= n_obs:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_points, replace=False))


def make_palette(n: int) -> list[str]:
    # Plotly qualitative palettes top out below the 75-label fine tissue case.
    import plotly.express as px

    colors = []
    for palette in [
        px.colors.qualitative.Plotly,
        px.colors.qualitative.Dark24,
        px.colors.qualitative.Light24,
        px.colors.qualitative.Alphabet,
    ]:
        colors.extend(palette)
    if n > len(colors):
        colors.extend(px.colors.sample_colorscale("Turbo", np.linspace(0, 1, n - len(colors))))
    return colors[:n]


def build_html(
    h5ad_path: Path,
    output_path: Path,
    embedding: str,
    color_by: str,
    adora_cache: Path | None,
    max_points: int | None,
    seed: int,
    marker_size: float,
    max_marker_size: float,
    adora_floor: float,
    adora_low_opacity: float,
    adora_high_opacity: float,
) -> None:
    with h5py.File(h5ad_path, "r") as f:
        coords_node = f[f"obsm/{embedding}"]
        n_obs = coords_node.shape[0]
        idx = choose_indices(n_obs, max_points, seed)
        coords = coords_node[idx, :2]

        color_codes, categories, color_values = categorical_codes(f, color_by, idx)
        metadata = {column: read_obs_column(f, column, idx) for column in METADATA_COLUMNS if column in f["obs"]}
        metadata_names = list(metadata)
        adora_expression = None
        adora_genes: list[str] = []
        if adora_cache is not None and adora_cache.exists():
            cached = np.load(adora_cache, allow_pickle=False)
            adora_expression = cached["expression"][idx]
            adora_genes = [str(gene) for gene in cached["genes"]]

        customdata = np.column_stack(
            [idx.astype(object), *[metadata[column] for column in metadata_names]]
        )

    used_codes = [int(code) for code in sorted(set(color_codes.tolist())) if code >= 0]
    palette = make_palette(len(used_codes))
    color_for_code = {code: palette[i] for i, code in enumerate(used_codes)}

    fig = go.Figure()
    categorical_trace_indices = []
    for code in used_codes:
        mask = color_codes == code
        label = categories[code]
        categorical_trace_indices.append(len(fig.data))
        fig.add_trace(
            go.Scattergl(
                x=coords[mask, 0],
                y=coords[mask, 1],
                mode="markers",
                name=label,
                marker={
                    "size": marker_size,
                    "color": color_for_code[code],
                    "opacity": 0.72,
                    "line": {"width": 0},
                },
                customdata=customdata[mask],
                hovertemplate=(
                    "<b>%{customdata[4]}</b><br>"
                    "cell index: %{customdata[0]}<br>"
                    "donor: %{customdata[1]}<br>"
                    "broad tissue: %{customdata[2]}<br>"
                    "fine tissue: %{customdata[3]}<br>"
                    "broad class: %{customdata[5]}<br>"
                    "assay: %{customdata[6]}<br>"
                    "sex: %{customdata[7]}<br>"
                    "age/stage: %{customdata[8]}<br>"
                    "<extra></extra>"
                ),
            )
        )

    adora_trace_indices: dict[str, list[int]] = {}
    if adora_expression is not None:
        for i, gene in enumerate(adora_genes):
            expr = adora_expression[:, i]
            low_mask = expr < adora_floor
            high_mask = ~low_mask
            positive = expr[expr > 0]
            cmax = float(np.quantile(positive, 0.99)) if len(positive) else 1.0
            adora_trace_indices[gene] = []
            adora_trace_indices[gene].append(len(fig.data))
            fig.add_trace(
                go.Scattergl(
                    x=coords[low_mask, 0],
                    y=coords[low_mask, 1],
                    mode="markers",
                    name=f"{gene} < {adora_floor:g}",
                    visible=False,
                    hoverinfo="skip",
                    marker={
                        "size": marker_size,
                        "color": "#d2d2d2",
                        "opacity": adora_low_opacity,
                        "line": {"width": 0},
                    },
                )
            )
            adora_trace_indices[gene].append(len(fig.data))
            fig.add_trace(
                go.Scattergl(
                    x=coords[high_mask, 0],
                    y=coords[high_mask, 1],
                    mode="markers",
                    name=f"{gene} >= {adora_floor:g}",
                    visible=False,
                    marker={
                        "size": marker_size,
                        "color": expr[high_mask],
                        "colorscale": "Viridis",
                        "cmin": adora_floor,
                        "cmax": cmax,
                        "opacity": adora_high_opacity,
                        "line": {"width": 0},
                        "colorbar": {"title": gene},
                    },
                    customdata=np.column_stack([idx[high_mask].astype(object), expr[high_mask].astype(object)]),
                    hovertemplate=(
                        f"<b>{gene}</b>: %{{customdata[1]:.4g}}<br>"
                        "cell index: %{customdata[0]}<br>"
                        "<extra></extra>"
                    ),
                )
            )

    shown = len(idx)
    title = f"Tabula Sapiens {embedding}, colored by {color_by} ({shown:,}/{n_obs:,} cells)"
    buttons = []
    if adora_trace_indices:
        n_traces = len(fig.data)
        cell_type_visible = [i in categorical_trace_indices for i in range(n_traces)]
        buttons.append(
            {
                "label": color_by,
                "method": "update",
                "args": [
                    {"visible": cell_type_visible},
                    {
                        "title": title,
                        "showlegend": True,
                    },
                ],
            }
        )
        for gene, trace_indices in adora_trace_indices.items():
            visible = [False] * n_traces
            for trace_idx in trace_indices:
                visible[trace_idx] = True
            buttons.append(
                {
                    "label": gene,
                    "method": "update",
                    "args": [
                        {"visible": visible},
                        {
                            "title": f"Tabula Sapiens {embedding}, colored by {gene} expression ({shown:,}/{n_obs:,} cells)",
                            "showlegend": False,
                        },
                    ],
                }
            )

    fig.update_layout(
        title=title,
        template="plotly_white",
        width=1250,
        height=850,
        dragmode="pan",
        legend={"itemsizing": "constant", "font": {"size": 10}},
        margin={"l": 35, "r": 260, "t": 70, "b": 35},
        xaxis={"title": "dim 1", "zeroline": False},
        yaxis={"title": "dim 2", "zeroline": False, "scaleanchor": "x", "scaleratio": 1},
        annotations=[
            {
                "text": "Click a cell to pin metadata below. Use the buttons to toggle cell type vs ADORA expression.",
                "xref": "paper",
                "yref": "paper",
                "x": 0,
                "y": 1.04,
                "showarrow": False,
                "font": {"size": 12},
                "align": "left",
            }
        ],
        updatemenus=(
            [
                {
                    "type": "buttons",
                    "direction": "right",
                    "buttons": buttons,
                    "x": 0,
                    "y": 1.11,
                    "xanchor": "left",
                    "yanchor": "top",
                    "pad": {"r": 8, "t": 4},
                    "showactive": True,
                }
            ]
            if buttons
            else []
        ),
    )

    field_names = ["cell_index", *metadata_names]
    x_range = float(np.nanmax(coords[:, 0]) - np.nanmin(coords[:, 0]))
    y_range = float(np.nanmax(coords[:, 1]) - np.nanmin(coords[:, 1]))
    initial_span = max(x_range, y_range)
    post_script = f"""
const fieldNames = {json.dumps(field_names)};
const baseMarkerSize = {marker_size};
const maxMarkerSize = {max_marker_size};
const initialSpan = {initial_span};
const gd = document.getElementById('{{plot_id}}');
const panel = document.createElement('pre');
panel.id = 'cell-info-panel';
panel.style.cssText = [
  'font: 13px/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace',
  'white-space: pre-wrap',
  'background: #f7f7f7',
  'border: 1px solid #d0d0d0',
  'border-radius: 6px',
  'padding: 12px',
  'margin: 12px 260px 24px 35px',
  'min-height: 92px'
].join(';');
panel.textContent = 'Click a point to show cell metadata here.';
gd.parentNode.insertBefore(panel, gd.nextSibling);
gd.on('plotly_click', function(eventData) {{
  const point = eventData.points[0];
  const values = point.customdata;
  if (!Array.isArray(values) || values.length !== fieldNames.length) {{
    const traceName = point.data.name || 'expression';
    const expr = Array.isArray(values) && values.length > 1 ? values[1] : values;
    panel.textContent = [
      `view: ${{traceName}}`,
      `cell_index: ${{Array.isArray(values) ? values[0] : 'unknown'}}`,
      `expression: ${{expr}}`,
      `embedding_x: ${{point.x}}`,
      `embedding_y: ${{point.y}}`
    ].join('\\n');
    return;
  }}
  const lines = fieldNames.map((name, i) => `${{name}}: ${{values[i]}}`);
  lines.push(`embedding_x: ${{point.x}}`);
  lines.push(`embedding_y: ${{point.y}}`);
  panel.textContent = lines.join('\\n');
}});

function currentSpan(layout) {{
  const xr = layout.xaxis.range;
  const yr = layout.yaxis.range;
  if (!xr || !yr) return initialSpan;
  return Math.max(Math.abs(xr[1] - xr[0]), Math.abs(yr[1] - yr[0]));
}}

let pendingSizeUpdate = false;
function updateMarkerSize() {{
  pendingSizeUpdate = false;
  const span = currentSpan(gd._fullLayout);
  const zoom = Math.max(1, initialSpan / span);
  const size = Math.min(maxMarkerSize, baseMarkerSize * Math.sqrt(zoom));
  Plotly.restyle(gd, {{'marker.size': size}});
}}

gd.on('plotly_relayout', function() {{
  if (!pendingSizeUpdate) {{
    pendingSizeUpdate = true;
    window.requestAnimationFrame(updateMarkerSize);
  }}
}});
"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_path, include_plotlyjs="cdn", full_html=True, post_script=post_script)
    print(f"Wrote {output_path}")
    print(title)
    print("Metadata columns:", ", ".join(metadata_names))
    if adora_trace_indices:
        print("ADORA expression toggles:", ", ".join(adora_trace_indices))
        print(
            f"ADORA low-expression layer: expression < {adora_floor:g}, "
            f"opacity {adora_low_opacity:g}; high layer opacity {adora_high_opacity:g}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--embedding", default="X_umap")
    parser.add_argument("--color-by", default="tissue_in_publication")
    parser.add_argument(
        "--adora-cache",
        type=Path,
        default=DEFAULT_ADORA_CACHE,
        help="Optional .npz cache with ADORA expression. Use 'none' to disable expression toggles.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=250_000,
        help="Number of cells to include. Use 0 for all cells, but expect a very large HTML file.",
    )
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--marker-size", type=float, default=2.0)
    parser.add_argument(
        "--max-marker-size",
        type=float,
        default=8.0,
        help="Largest screen-pixel marker size after zoom-responsive scaling.",
    )
    parser.add_argument(
        "--adora-floor",
        type=float,
        default=0.04,
        help="ADORA expression below this value is drawn as a very faint context layer.",
    )
    parser.add_argument("--adora-low-opacity", type=float, default=0.025)
    parser.add_argument("--adora-high-opacity", type=float, default=0.92)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    max_points = None if args.max_points == 0 else args.max_points
    adora_cache = None if str(args.adora_cache).lower() == "none" else args.adora_cache
    build_html(
        h5ad_path=args.h5ad,
        output_path=args.out,
        embedding=args.embedding,
        color_by=args.color_by,
        adora_cache=adora_cache,
        max_points=max_points,
        seed=args.seed,
        marker_size=args.marker_size,
        max_marker_size=args.max_marker_size,
        adora_floor=args.adora_floor,
        adora_low_opacity=args.adora_low_opacity,
        adora_high_opacity=args.adora_high_opacity,
    )


if __name__ == "__main__":
    main()
