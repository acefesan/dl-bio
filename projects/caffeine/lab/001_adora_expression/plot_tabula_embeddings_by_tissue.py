#!/usr/bin/env python
"""Plot all Tabula Sapiens H5AD embeddings colored by tissue."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D


LAB_DIR = Path(__file__).resolve().parent
DEFAULT_H5AD = LAB_DIR / "cache" / "tabula_sapiens_all_cells.h5ad"
DEFAULT_OUT = LAB_DIR / "figures" / "tabula_sapiens_embeddings_by_tissue.png"


def decode_labels(values: np.ndarray) -> list[str]:
    labels: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            labels.append(value.decode("utf-8"))
        else:
            labels.append(str(value))
    return labels


def read_categorical(f: h5py.File, obs_column: str) -> tuple[np.ndarray, list[str]]:
    node = f["obs"][obs_column]
    if not isinstance(node, h5py.Group) or "codes" not in node or "categories" not in node:
        raise ValueError(f"obs/{obs_column!r} is not an AnnData categorical column")
    codes = node["codes"][:]
    categories = decode_labels(node["categories"][:])
    return codes, categories


def embedding_names(f: h5py.File) -> list[str]:
    names: list[str] = []
    for name, node in f["obsm"].items():
        if isinstance(node, h5py.Dataset) and node.ndim == 2 and node.shape[1] >= 2:
            names.append(name)
    return sorted(names)


def make_palette(n: int) -> ListedColormap:
    base_names = ["tab20", "tab20b", "tab20c"]
    colors = []
    for name in base_names:
        cmap = plt.get_cmap(name)
        colors.extend(cmap(i) for i in range(cmap.N))
    if n > len(colors):
        colors.extend(plt.get_cmap("hsv")(np.linspace(0, 1, n - len(colors), endpoint=False)))
    return ListedColormap(colors[:n])


def maybe_subsample(n_obs: int, max_points: int | None, seed: int) -> np.ndarray | slice:
    if max_points is None or max_points >= n_obs:
        return slice(None)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_points, replace=False))


def plot_embeddings(
    h5ad_path: Path,
    output_path: Path,
    color_by: str,
    max_points: int | None,
    seed: int,
    point_size: float,
    alpha: float,
    dpi: int,
) -> None:
    with h5py.File(h5ad_path, "r") as f:
        names = embedding_names(f)
        if not names:
            raise ValueError(f"No 2D-or-wider embeddings found in {h5ad_path}")

        codes, categories = read_categorical(f, color_by)
        n_obs = codes.shape[0]
        idx = maybe_subsample(n_obs, max_points, seed)
        plot_codes = codes[idx]

        valid_codes = np.array(sorted(set(int(c) for c in plot_codes if c >= 0)))
        used_categories = [categories[i] for i in valid_codes]
        code_to_dense = {code: i for i, code in enumerate(valid_codes)}
        dense_codes = np.array([code_to_dense.get(int(code), -1) for code in plot_codes], dtype=np.int16)
        cmap = make_palette(len(used_categories))

        ncols = 3 if len(names) > 2 else len(names)
        nrows = int(np.ceil(len(names) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.4 * ncols, 4.8 * nrows), squeeze=False)
        axes_flat = axes.ravel()

        for ax, name in zip(axes_flat, names, strict=False):
            coords = f["obsm"][name][idx, :2]
            draw_order = np.random.default_rng(seed).permutation(coords.shape[0])
            ax.scatter(
                coords[draw_order, 0],
                coords[draw_order, 1],
                c=dense_codes[draw_order],
                cmap=cmap,
                s=point_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )
            ax.set_title(name, fontsize=11)
            ax.set_xlabel("dim 1")
            ax.set_ylabel("dim 2")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_linewidth(0.6)
                spine.set_color("#b8b8b8")

        for ax in axes_flat[len(names) :]:
            ax.axis("off")

        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=5,
                markerfacecolor=cmap(i),
                markeredgecolor="none",
                label=label,
            )
            for i, label in enumerate(used_categories)
        ]
        fig.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(0.875, 0.5),
            frameon=False,
            title=color_by,
            fontsize=8,
            title_fontsize=9,
            markerscale=1.2,
        )
        shown = n_obs if isinstance(idx, slice) else len(idx)
        fig.suptitle(
            f"Tabula Sapiens embeddings colored by {color_by} ({shown:,}/{n_obs:,} cells)",
            fontsize=14,
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 0.86, 0.965))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {output_path}")
        print("Embeddings:", ", ".join(names))
        print(f"Color categories ({len(used_categories)}):", ", ".join(used_categories))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--color-by",
        default="tissue_in_publication",
        help="Categorical obs column. Use 'tissue' for the finer 75-category label.",
    )
    parser.add_argument("--max-points", type=int, default=None, help="Optional random subsample for faster drafts.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--point-size", type=float, default=0.05)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot_embeddings(
        h5ad_path=args.h5ad,
        output_path=args.out,
        color_by=args.color_by,
        max_points=args.max_points,
        seed=args.seed,
        point_size=args.point_size,
        alpha=args.alpha,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
