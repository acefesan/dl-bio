#!/usr/bin/env python3
"""Extract proteins for 8 root HOGs, compute ESM2 embeddings, and plot UMAPs.

Usage:
    python hog_umap.py                  # Run all steps
    python hog_umap.py extract          # Step 1: extract proteins + sequences
    python hog_umap.py embed            # Step 2: compute embeddings (all models)
    python hog_umap.py plot             # Step 3: UMAP plots
"""

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import umap
from Bio import SeqIO
from tqdm import tqdm
from transformers import AutoModel, AutoTokenizer

PROJECT_ROOT = Path(__file__).parent.parent.parent
MAMMALIA_DIR = PROJECT_ROOT / "assets" / "proteins" / "mammalia"
OUTPUT_DIR = MAMMALIA_DIR / "hog_study"

# 8 root HOGs: 2 vision, 1 lens, 2 neuron, 1 skin, 2 muscle
TARGET_HOGS = {
    "HOG:E0754125": "Rhodopsin",
    "HOG:E0747130": "Cone opsin LW",
    "HOG:E0736973": "Crystallin gamma",
    "HOG:E0801852": "Sodium channel",
    "HOG:E0781053": "Synaptotagmin",
    "HOG:E0738002": "Keratin type I",
    "HOG:E0793067": "Myosin heavy chain",
    "HOG:E1027835": "Actin",
}

# 8 colors — one per HOG (colorblind-friendly palette)
HOG_COLORS = {
    "Rhodopsin": "#E69F00",        # orange
    "Cone opsin LW": "#D55E00",    # vermillion
    "Crystallin gamma": "#56B4E9",  # sky blue
    "Sodium channel": "#009E73",    # bluish green
    "Synaptotagmin": "#CC79A7",     # reddish purple
    "Keratin type I": "#F0E442",    # yellow
    "Myosin heavy chain": "#0072B2",# blue
    "Actin": "#999999",             # grey
}

# 5 taxa with distinct shapes (pick representatives across mammalian diversity)
SHAPE_TAXA = {
    "HUMAN": "o",   # circle — Primates
    "MOUSE": "s",   # square — Rodentia
    "BOVIN": "^",   # triangle up — Artiodactyla
    "CANLF": "D",   # diamond — Carnivora
    "ORNAN": "P",   # plus — Monotremata (platypus)
}

ESM2_MODELS = [
    "facebook/esm2_t6_8M_UR50D",
    "facebook/esm2_t30_150M_UR50D",
    "facebook/esm2_t33_650M_UR50D",
    "facebook/esm2_t36_3B_UR50D",
]


# =============================================================================
# Step 1: Extract proteins and sequences
# =============================================================================

def step_extract():
    """Extract proteins matching the 8 HOGs and retrieve sequences from bulk FASTA."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lists_dir = MAMMALIA_DIR / "protein_lists"
    bulk_fasta = MAMMALIA_DIR / "oma-seqs.fa.gz"

    # Load species metadata
    with open(MAMMALIA_DIR / "species.json") as f:
        species_meta = json.load(f)

    # Collect matching proteins from all species
    records = []
    for feather_path in sorted(lists_dir.glob("*.feather")):
        code = feather_path.stem
        df = pd.read_feather(feather_path)
        for hog_id, hog_name in TARGET_HOGS.items():
            mask = df["oma_hog_id"].str.startswith(hog_id, na=False)
            matched = df[mask].copy()
            if len(matched) > 0:
                matched["species_code"] = code
                matched["hog_name"] = hog_name
                matched["roothog_id"] = hog_id
                meta = species_meta.get(code, {})
                matched["taxon_id"] = meta.get("taxon_id")
                matched["scientific_name"] = meta.get("species")
                records.append(matched)

    proteins = pd.concat(records, ignore_index=True)
    print(f"Found {len(proteins):,} proteins across {proteins['species_code'].nunique()} species")
    print(f"HOG breakdown:")
    for name, count in proteins["hog_name"].value_counts().items():
        print(f"  {name}: {count}")

    # Extract sequences from bulk FASTA
    target_ids = set(proteins["omaid"])
    sequences = {}
    print(f"\nScanning oma-seqs.fa.gz for {len(target_ids):,} IDs...")
    with gzip.open(bulk_fasta, "rt") as f:
        for record in tqdm(SeqIO.parse(f, "fasta"), desc="Scanning FASTA",
                           unit=" seqs"):
            if record.id in target_ids:
                sequences[record.id] = str(record.seq)
                if len(sequences) >= len(target_ids):
                    break

    print(f"Found {len(sequences):,} / {len(target_ids):,} sequences")

    proteins["sequence"] = proteins["omaid"].map(sequences)
    proteins = proteins[proteins["sequence"].notna()].reset_index(drop=True)
    proteins["sequence_length"] = proteins["sequence"].str.len()

    cols = ["omaid", "species_code", "taxon_id", "scientific_name",
            "oma_hog_id", "roothog_id", "hog_name", "sequence",
            "sequence_length"]
    proteins[cols].to_feather(OUTPUT_DIR / "hog_proteins.feather")

    print(f"\nSaved {len(proteins):,} proteins with sequences to hog_proteins.feather")
    return proteins


# =============================================================================
# Step 2: Compute embeddings
# =============================================================================

def compute_embeddings_for_model(model_name: str, sequences: list[str],
                                 batch_size: int = 8, max_len: int = 1024):
    """Compute mean-pooled embeddings for a list of sequences."""
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\n  Loading {model_name} on {device}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device).eval()

    # Truncate to max_len
    truncated = [s[:max_len] for s in sequences]

    all_embeddings = []
    for i in tqdm(range(0, len(truncated), batch_size),
                  desc=f"  {model_name.split('/')[-1]}"):
        batch = truncated[i:i + batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True,
                           truncation=True, max_length=max_len + 2).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        # Mean pool over sequence length (exclude special tokens)
        attention_mask = inputs["attention_mask"]
        hidden = outputs.last_hidden_state
        # Mask padding
        mask_expanded = attention_mask.unsqueeze(-1).float()
        summed = (hidden * mask_expanded).sum(dim=1)
        counts = mask_expanded.sum(dim=1)
        mean_pooled = summed / counts
        all_embeddings.append(mean_pooled.cpu().numpy())

    del model
    torch.cuda.empty_cache()

    return np.vstack(all_embeddings)


def step_embed():
    """Compute embeddings for all ESM2 models."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    proteins = pd.read_feather(OUTPUT_DIR / "hog_proteins.feather")
    sequences = proteins["sequence"].tolist()

    print(f"Computing embeddings for {len(sequences):,} proteins")

    for model_name in ESM2_MODELS:
        short_name = model_name.split("/")[-1]
        out_path = OUTPUT_DIR / f"embeddings_{short_name}.npy"
        if out_path.exists():
            print(f"\n  {short_name}: already exists, skipping")
            continue

        # Adjust batch size by model size
        if "8M" in short_name:
            batch_size = 64
        elif "150M" in short_name:
            batch_size = 32
        elif "650M" in short_name:
            batch_size = 16
        elif "3B" in short_name:
            batch_size = 4
        else:
            batch_size = 8

        embeddings = compute_embeddings_for_model(model_name, sequences,
                                                   batch_size=batch_size)
        np.save(out_path, embeddings)
        print(f"  Saved {out_path.name}: shape {embeddings.shape}")

    print("\nAll embeddings computed.")


# =============================================================================
# Step 3: UMAP plots
# =============================================================================

def step_plot():
    """Create UMAP plots for each ESM2 model."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    proteins = pd.read_feather(OUTPUT_DIR / "hog_proteins.feather")

    # Pre-compute which rows belong to shape taxa vs "other"
    shape_mask = proteins["species_code"].isin(SHAPE_TAXA)

    for model_name in ESM2_MODELS:
        short_name = model_name.split("/")[-1]
        emb_path = OUTPUT_DIR / f"embeddings_{short_name}.npy"
        if not emb_path.exists():
            print(f"Skipping {short_name}: no embeddings found")
            continue

        print(f"\nPlotting {short_name}...")
        embeddings = np.load(emb_path)

        # UMAP reduction
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine",
                            random_state=42)
        coords = reducer.fit_transform(embeddings)

        fig, ax = plt.subplots(figsize=(12, 9))

        # Plot "other" taxa first (small, faded dots)
        other_idx = ~shape_mask
        if other_idx.any():
            other_colors = [HOG_COLORS[h] for h in proteins.loc[other_idx, "hog_name"]]
            ax.scatter(coords[other_idx, 0], coords[other_idx, 1],
                       c=other_colors, marker=".", s=15, alpha=0.3,
                       linewidths=0, zorder=1)

        # Plot shape taxa on top (larger, distinct markers)
        for species_code, marker in SHAPE_TAXA.items():
            species_mask = proteins["species_code"] == species_code
            if not species_mask.any():
                continue
            idx = species_mask.values
            colors = [HOG_COLORS[h] for h in proteins.loc[idx, "hog_name"]]
            ax.scatter(coords[idx, 0], coords[idx, 1],
                       c=colors, marker=marker, s=60, alpha=0.85,
                       edgecolors="black", linewidths=0.5, zorder=2,
                       label=species_code)

        # HOG color legend
        from matplotlib.patches import Patch
        hog_handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.5, label=n)
                       for n, c in HOG_COLORS.items()]

        # Shape legend
        from matplotlib.lines import Line2D
        shape_handles = [
            Line2D([0], [0], marker=m, color="w", markerfacecolor="gray",
                   markeredgecolor="black", markersize=8, linewidth=0,
                   label=f"{code}")
            for code, m in SHAPE_TAXA.items()
        ]

        leg1 = ax.legend(handles=hog_handles, title="Root HOG",
                         loc="upper left", fontsize=8, title_fontsize=9,
                         framealpha=0.9)
        ax.add_artist(leg1)
        ax.legend(handles=shape_handles, title="Species",
                  loc="lower left", fontsize=8, title_fontsize=9,
                  framealpha=0.9)

        ax.set_title(f"UMAP of 8 Mammalian Root HOGs — {short_name}",
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])

        fig.tight_layout()
        out_path = OUTPUT_DIR / f"umap_{short_name}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {out_path.name}")

    print("\nAll plots saved.")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extract 8 HOGs, compute ESM2 embeddings, plot UMAPs"
    )
    parser.add_argument("step", nargs="?", default="all",
                        choices=["extract", "embed", "plot", "all"])
    args = parser.parse_args()

    steps = {
        "extract": step_extract,
        "embed": step_embed,
        "plot": step_plot,
    }

    if args.step == "all":
        for name, fn in steps.items():
            print(f"\n{'=' * 60}")
            print(f"Step: {name}")
            print(f"{'=' * 60}")
            fn()
    else:
        steps[args.step]()

    print("\nDone!")


if __name__ == "__main__":
    main()
