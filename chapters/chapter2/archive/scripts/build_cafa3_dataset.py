#!/usr/bin/env python3
"""Build merged CAFA3 dataset with pre-computed ESM2 embeddings and HOGs.

This script loads the CAFA3 protein dataset components, merges them with
pre-computed ESM2 embeddings and Hierarchical Orthologous Groups (HOGs) from OMA,
and saves the result with provenance metadata.

Usage:
    python build_cafa3_dataset.py                    # Full dataset
    python build_cafa3_dataset.py --skip-taxonomy    # Skip UniProt API calls
    python build_cafa3_dataset.py --fetch-hogs       # Fetch HOGs from OMA (slow)
    python build_cafa3_dataset.py --hog-cache hogs.csv  # Use cached HOG data
    python build_cafa3_dataset.py --output-dir /path # Custom output directory

Based on: chapter2.ipynb
"""

import argparse
import json
import time
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from tqdm import tqdm


# =============================================================================
# Configuration
# =============================================================================

BASE_PATH = Path(__file__).parent.parent
ASSETS_PATH = BASE_PATH / "assets/proteins/datasets"

# Input files
SEQUENCES_FILE = ASSETS_PATH / "train_sequences.fasta"
TAXONOMY_FILE = ASSETS_PATH / "train_taxonomy.tsv.zip"
LABELS_FILE = ASSETS_PATH / "train_terms.tsv.zip"
GO_DESCRIPTIONS_FILE = ASSETS_PATH / "go_term_descriptions.csv"
EMBEDDINGS_FILE = ASSETS_PATH / "all_species_embeddings.feather"

# Embedding model used
EMBEDDING_MODEL = "facebook/esm2_t33_650M_UR50D"


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_sequences(fasta_path: Path) -> pd.DataFrame:
    """Load protein sequences from FASTA file."""
    print(f"Loading sequences from {fasta_path.name}...")

    data = []
    with open(fasta_path) as f:
        for fasta in SeqIO.parse(f, "fasta"):
            data.append({
                "EntryID": fasta.id,
                "Sequence": str(fasta.seq),
                "Length": len(fasta.seq),
            })

    df = pd.DataFrame(data)
    print(f"  Loaded {len(df):,} sequences")
    return df


def load_taxonomy(taxonomy_path: Path) -> pd.DataFrame:
    """Load taxonomy data (EntryID -> taxonomyID mapping)."""
    print(f"Loading taxonomy from {taxonomy_path.name}...")
    df = pd.read_csv(taxonomy_path, sep="\t", compression="infer")
    print(f"  Loaded {len(df):,} taxonomy entries")
    return df


def fetch_taxonomy_names_batch(
    tax_ids: list[int],
    batch_size: int = 100
) -> pd.DataFrame:
    """Fetch species names from UniProt in batches."""
    url = "https://rest.uniprot.org/taxonomy/stream"
    all_results = []

    for i in tqdm(range(0, len(tax_ids), batch_size), desc="Fetching taxonomy names"):
        batch = tax_ids[i:i + batch_size]
        params = {
            "query": " OR ".join(f"id:{tid}" for tid in batch),
            "fields": "id,scientific_name,common_name",
            "format": "tsv"
        }
        try:
            response = requests.get(url, params=params, timeout=30)
            if response.ok and not response.text.startswith("<"):
                df = pd.read_csv(StringIO(response.text), sep="\t")
                all_results.append(df)
        except Exception as e:
            print(f"  Batch {i // batch_size} failed: {e}")
            continue

    if all_results:
        result = pd.concat(all_results, ignore_index=True)
        result.columns = ["taxonomyID", "scientific_name", "common_name"]
        result["taxonomyID"] = result["taxonomyID"].astype(int)
        return result

    return pd.DataFrame(columns=["taxonomyID", "scientific_name", "common_name"])


def load_labels(labels_path: Path) -> pd.DataFrame:
    """Load GO term labels."""
    print(f"Loading GO term labels from {labels_path.name}...")
    df = pd.read_csv(labels_path, sep="\t", compression="infer")
    print(f"  Loaded {len(df):,} label entries")
    print(f"  GO aspects: {df['aspect'].value_counts().to_dict()}")
    return df


def load_go_descriptions(descriptions_path: Path) -> pd.DataFrame:
    """Load GO term descriptions."""
    print(f"Loading GO descriptions from {descriptions_path.name}...")
    df = pd.read_csv(descriptions_path)
    print(f"  Loaded {len(df):,} GO term descriptions")
    return df


def load_embeddings(embeddings_path: Path) -> pd.DataFrame:
    """Load pre-computed ESM2 embeddings."""
    print(f"Loading embeddings from {embeddings_path.name}...")
    df = pd.read_feather(embeddings_path)

    # Keep only EntryID and embedding columns (ME:*)
    embedding_cols = [c for c in df.columns if c.startswith("ME:")]
    df = df[["EntryID"] + embedding_cols]

    print(f"  Loaded {len(df):,} embeddings with {len(embedding_cols)} dimensions")
    return df


# =============================================================================
# HOG (Hierarchical Orthologous Groups) Functions
# =============================================================================

def fetch_hogs_from_oma(
    entry_ids: list[str],
    rate_limit_delay: float = 0.2,
    cache_path: Path | None = None,
    save_interval: int = 500,
) -> pd.DataFrame:
    """Fetch HOG IDs for UniProt accessions using OMA API.

    Note: Uses single-protein lookups since the bulk API requires OMA IDs.
    This is slower but works directly with UniProt accession IDs.

    Args:
        entry_ids: List of UniProt accession IDs
        rate_limit_delay: Delay between requests (seconds)
        cache_path: Path to save incremental results (CSV)
        save_interval: Save to cache every N proteins

    Returns:
        DataFrame with columns: EntryID, oma_id, hog_id, roothog_id
    """
    base_url = "https://omabrowser.org/api"
    results = []
    not_found = 0
    errors = 0

    # Check for existing partial results
    already_fetched = set()
    if cache_path and cache_path.exists():
        existing_df = pd.read_csv(cache_path)
        already_fetched = set(existing_df["EntryID"].tolist())
        results = existing_df.to_dict("records")
        print(f"  Resuming from cache: {len(already_fetched):,} already fetched")

    # Filter out already-fetched entries
    remaining_ids = [eid for eid in entry_ids if eid not in already_fetched]

    print(f"Fetching HOGs for {len(remaining_ids):,} proteins from OMA...")
    if len(already_fetched) > 0:
        print(f"  (Skipping {len(already_fetched):,} already cached)")
    print(f"  Rate limit delay: {rate_limit_delay}s")
    print(f"  Estimated time: ~{len(remaining_ids) * rate_limit_delay / 60:.1f} minutes")
    if cache_path:
        print(f"  Saving to cache every {save_interval} proteins: {cache_path}")

    for i, entry_id in enumerate(tqdm(remaining_ids, desc="Fetching HOGs")):
        try:
            response = requests.get(
                f"{base_url}/protein/{entry_id}/",
                timeout=30
            )

            if response.ok:
                data = response.json()
                results.append({
                    "EntryID": entry_id,
                    "oma_id": data.get("omaid"),
                    "hog_id": data.get("oma_hog_id"),  # Full HOG ID with hierarchy
                    "roothog_id": data.get("roothog_id"),  # Root HOG (LUCA level)
                })
            elif response.status_code == 404:
                not_found += 1
                # Still record as not found so we don't retry
                results.append({
                    "EntryID": entry_id,
                    "oma_id": None,
                    "hog_id": None,
                    "roothog_id": None,
                })
            elif response.status_code == 429:  # Rate limited
                print(f"\n  Rate limited, waiting 10s...")
                time.sleep(10)
                # Retry
                response = requests.get(f"{base_url}/protein/{entry_id}/", timeout=30)
                if response.ok:
                    data = response.json()
                    results.append({
                        "EntryID": entry_id,
                        "oma_id": data.get("omaid"),
                        "hog_id": data.get("oma_hog_id"),
                        "roothog_id": data.get("roothog_id"),
                    })
            else:
                errors += 1

        except requests.exceptions.Timeout:
            errors += 1
        except Exception as e:
            errors += 1
            if errors <= 5:
                print(f"\n  Error for {entry_id}: {e}")

        time.sleep(rate_limit_delay)

        # Incremental save
        if cache_path and (i + 1) % save_interval == 0:
            pd.DataFrame(results).to_csv(cache_path, index=False)
            tqdm.write(f"  [Saved {len(results):,} entries to cache]")

    # Final save
    if cache_path and results:
        pd.DataFrame(results).to_csv(cache_path, index=False)

    if not results:
        print("  Warning: No HOG data retrieved from OMA")
        return pd.DataFrame(columns=["EntryID", "oma_id", "hog_id", "roothog_id"])

    df = pd.DataFrame(results)

    print(f"  Retrieved HOG info for {len(df):,} proteins")
    print(f"  Proteins with HOG assignment: {df['hog_id'].notna().sum():,}")
    if not_found > 0:
        print(f"  Proteins not found in OMA: {not_found:,}")
    if errors > 0:
        print(f"  Errors: {errors:,}")

    return df


def load_hog_cache(cache_path: Path) -> pd.DataFrame | None:
    """Load cached HOG data if available."""
    if cache_path.exists():
        print(f"Loading HOG cache from {cache_path.name}...")
        df = pd.read_csv(cache_path)
        print(f"  Loaded {len(df):,} cached HOG entries")
        return df
    return None


def save_hog_cache(hog_df: pd.DataFrame, cache_path: Path) -> None:
    """Save HOG data to cache file."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    hog_df.to_csv(cache_path, index=False)
    print(f"  Saved HOG cache to {cache_path}")


# =============================================================================
# Dataset Building
# =============================================================================

def build_dataset(
    fetch_taxonomy_names: bool = True,
    fetch_hogs: bool = False,
    hog_cache_path: Path | None = None,
    output_dir: Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Build the merged CAFA3 dataset with embeddings and optional HOGs.

    Args:
        fetch_taxonomy_names: Whether to fetch species names from UniProt
        fetch_hogs: Whether to fetch HOGs from OMA (slow)
        hog_cache_path: Path to cached HOG data (load if exists, save if fetched)
        output_dir: Output directory for saving HOG cache

    Returns:
        Tuple of (merged_dataframe, statistics_dict)
    """
    stats = {}

    # Load all components
    sequence_df = load_sequences(SEQUENCES_FILE)
    stats["n_sequences"] = len(sequence_df)

    taxonomy_df = load_taxonomy(TAXONOMY_FILE)
    stats["n_taxonomy_entries"] = len(taxonomy_df)

    # Optionally fetch taxonomy names from UniProt
    if fetch_taxonomy_names:
        unique_tax_ids = taxonomy_df["taxonomyID"].unique().tolist()
        print(f"\nFetching species names for {len(unique_tax_ids):,} taxonomy IDs...")
        taxonomy_names = fetch_taxonomy_names_batch(unique_tax_ids, batch_size=100)
        print(f"  Retrieved {len(taxonomy_names):,} taxonomy names")

        taxonomy_df = taxonomy_df.merge(taxonomy_names, on="taxonomyID", how="left")
        stats["n_taxonomy_names_fetched"] = len(taxonomy_names)
    else:
        print("\nSkipping taxonomy name fetch (--skip-taxonomy)")
        taxonomy_df["scientific_name"] = None
        taxonomy_df["common_name"] = None

    labels_df = load_labels(LABELS_FILE)
    stats["n_label_entries"] = len(labels_df)

    go_descriptions_df = load_go_descriptions(GO_DESCRIPTIONS_FILE)
    stats["n_go_terms"] = len(go_descriptions_df)

    embeddings_df = load_embeddings(EMBEDDINGS_FILE)
    stats["n_embeddings"] = len(embeddings_df)
    stats["embedding_dim"] = len([c for c in embeddings_df.columns if c.startswith("ME:")])

    # Load or fetch HOGs
    hog_df = None
    if hog_cache_path and hog_cache_path.exists():
        hog_df = load_hog_cache(hog_cache_path)
    elif fetch_hogs:
        print("\nFetching HOGs from OMA browser...")
        unique_entries = sequence_df["EntryID"].unique().tolist()
        # Use cache path for incremental saving
        incremental_cache = hog_cache_path or (output_dir / "hog_cache.csv" if output_dir else None)
        hog_df = fetch_hogs_from_oma(unique_entries, cache_path=incremental_cache)

        # Save to cache
        if hog_cache_path:
            save_hog_cache(hog_df, hog_cache_path)
        elif output_dir:
            default_cache = output_dir / "hog_cache.csv"
            save_hog_cache(hog_df, default_cache)
    else:
        print("\nSkipping HOG fetch (use --fetch-hogs or --hog-cache)")

    if hog_df is not None:
        stats["n_hog_entries"] = len(hog_df)
        stats["n_proteins_with_hog"] = int(hog_df["hog_id"].notna().sum())
        stats["n_unique_hogs"] = int(hog_df["hog_id"].nunique())
        stats["n_unique_roothogs"] = int(hog_df["roothog_id"].nunique())

    # Merge datasets
    print("\nMerging datasets...")

    # Step 1: Sequences + Taxonomy
    merged_df = sequence_df.merge(taxonomy_df, on="EntryID", how="left")
    print(f"  After taxonomy merge: {len(merged_df):,} rows")

    # Step 2: Merge with HOGs (if available)
    if hog_df is not None and len(hog_df) > 0:
        hog_cols = ["EntryID", "oma_id", "hog_id", "roothog_id"]
        hog_df_dedup = hog_df[hog_cols].drop_duplicates(subset=["EntryID"])
        merged_df = merged_df.merge(hog_df_dedup, on="EntryID", how="left")
        n_with_hog = merged_df["hog_id"].notna().sum()
        print(f"  After HOG merge: {len(merged_df):,} rows ({n_with_hog:,} with HOG)")
    else:
        merged_df["oma_id"] = None
        merged_df["hog_id"] = None
        merged_df["roothog_id"] = None

    # Step 3: Add GO term descriptions to labels
    labels_df = labels_df.merge(go_descriptions_df, on="term", how="left")

    # Step 4: Merge with labels (this expands rows - one per protein-GO term pair)
    merged_df = merged_df.merge(labels_df, on="EntryID", how="inner")
    print(f"  After labels merge: {len(merged_df):,} rows")

    # Step 5: Merge with embeddings
    merged_df = merged_df.merge(embeddings_df, on="EntryID", how="left")
    n_with_embeddings = merged_df[merged_df.columns[merged_df.columns.str.startswith("ME:")]].notna().any(axis=1).sum()
    print(f"  After embeddings merge: {len(merged_df):,} rows ({n_with_embeddings:,} with embeddings)")

    stats["n_final_rows"] = len(merged_df)
    stats["n_unique_proteins"] = merged_df["EntryID"].nunique()
    stats["n_unique_go_terms"] = merged_df["term"].nunique()
    stats["n_rows_with_embeddings"] = int(n_with_embeddings)
    stats["n_rows_with_hog"] = int(merged_df["hog_id"].notna().sum())

    return merged_df, stats


def create_metadata(stats: dict, output_path: Path, hog_included: bool = False) -> dict:
    """Create provenance metadata for the dataset."""
    metadata = {
        "name": "CAFA3 Dataset with ESM2 Embeddings" + (" and HOGs" if hog_included else ""),
        "description": "Merged CAFA3 protein function prediction dataset with pre-computed ESM2 embeddings"
                       + (" and Hierarchical Orthologous Groups (HOGs) from OMA" if hog_included else ""),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "generator_script": "chapters/build_cafa3_dataset.py",
        "output_file": str(output_path.name),
        "source_datasets": {
            "sequences": {
                "file": str(SEQUENCES_FILE.relative_to(BASE_PATH)),
                "description": "CAFA3 training protein sequences (FASTA format)",
            },
            "taxonomy": {
                "file": str(TAXONOMY_FILE.relative_to(BASE_PATH)),
                "description": "Protein to taxonomy ID mapping",
            },
            "labels": {
                "file": str(LABELS_FILE.relative_to(BASE_PATH)),
                "description": "Protein GO term annotations (BPO, MFO, CCO)",
            },
            "go_descriptions": {
                "file": str(GO_DESCRIPTIONS_FILE.relative_to(BASE_PATH)),
                "description": "GO term ID to description mapping",
            },
            "embeddings": {
                "file": str(EMBEDDINGS_FILE.relative_to(BASE_PATH)),
                "description": "Pre-computed ESM2 mean-pooled embeddings",
            },
        },
        "embedding_model": {
            "name": EMBEDDING_MODEL,
            "type": "ESM2 Protein Language Model",
            "source": "https://github.com/facebookresearch/esm",
            "pooling": "mean",
            "dimensions": stats.get("embedding_dim", 1280),
        },
        "statistics": stats,
        "columns": {
            "EntryID": "UniProt accession ID",
            "Sequence": "Amino acid sequence",
            "Length": "Sequence length in amino acids",
            "taxonomyID": "NCBI taxonomy ID",
            "scientific_name": "Species scientific name (from UniProt)",
            "common_name": "Species common name (from UniProt)",
            "oma_id": "OMA database identifier",
            "hog_id": "Hierarchical Orthologous Group ID (species-specific level)",
            "roothog_id": "Root HOG ID (LUCA level, most ancient)",
            "term": "GO term ID (e.g., GO:0003674)",
            "aspect": "GO aspect: BPO (Biological Process), MFO (Molecular Function), CCO (Cellular Component)",
            "description": "GO term description",
            "ME:1-ME:1280": "ESM2 embedding dimensions",
        },
    }

    if hog_included:
        metadata["source_datasets"]["hogs"] = {
            "source": "OMA Browser API (https://omabrowser.org/api)",
            "description": "Hierarchical Orthologous Groups mapping proteins to evolutionary lineages",
        }

    return metadata


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build merged CAFA3 dataset with ESM2 embeddings and HOGs"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: assets/proteins/datasets/cafa3_merged)"
    )
    parser.add_argument(
        "--skip-taxonomy",
        action="store_true",
        help="Skip fetching taxonomy names from UniProt API"
    )
    parser.add_argument(
        "--fetch-hogs",
        action="store_true",
        help="Fetch HOGs from OMA browser API (slow, ~24h for full dataset)"
    )
    parser.add_argument(
        "--hog-cache",
        type=str,
        default=None,
        help="Path to HOG cache CSV file (load if exists, save if fetching)"
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["feather", "parquet", "csv"],
        default="feather",
        help="Output format (default: feather)"
    )
    args = parser.parse_args()

    # Set output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = ASSETS_PATH / "cafa3_merged"

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Set HOG cache path
    hog_cache_path = Path(args.hog_cache) if args.hog_cache else None

    # Build dataset
    print("\n" + "=" * 60)
    print("Building CAFA3 dataset with embeddings" + (" and HOGs" if args.fetch_hogs or hog_cache_path else ""))
    print("=" * 60 + "\n")

    merged_df, stats = build_dataset(
        fetch_taxonomy_names=not args.skip_taxonomy,
        fetch_hogs=args.fetch_hogs,
        hog_cache_path=hog_cache_path,
        output_dir=output_dir,
    )

    # Save dataset
    print("\n" + "=" * 60)
    print("Saving dataset")
    print("=" * 60 + "\n")

    if args.format == "feather":
        output_file = output_dir / "cafa3_with_embeddings.feather"
        merged_df.to_feather(output_file)
    elif args.format == "parquet":
        output_file = output_dir / "cafa3_with_embeddings.parquet"
        merged_df.to_parquet(output_file)
    else:
        output_file = output_dir / "cafa3_with_embeddings.csv"
        merged_df.to_csv(output_file, index=False)

    print(f"Saved dataset to: {output_file}")
    print(f"  File size: {output_file.stat().st_size / 1024**2:.1f} MB")

    # Save metadata
    hog_included = args.fetch_hogs or (hog_cache_path is not None and hog_cache_path.exists())
    metadata = create_metadata(stats, output_file, hog_included=bool(hog_included))
    metadata_file = output_dir / "metadata.json"
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved metadata to: {metadata_file}")

    # Print summary
    print("\n" + "=" * 60)
    print("Dataset Summary")
    print("=" * 60)
    print(f"  Total rows: {stats['n_final_rows']:,}")
    print(f"  Unique proteins: {stats['n_unique_proteins']:,}")
    print(f"  Unique GO terms: {stats['n_unique_go_terms']:,}")
    print(f"  Rows with embeddings: {stats['n_rows_with_embeddings']:,}")
    print(f"  Rows with HOG: {stats.get('n_rows_with_hog', 0):,}")
    if "n_unique_hogs" in stats:
        print(f"  Unique HOGs: {stats['n_unique_hogs']:,}")
        print(f"  Unique root HOGs: {stats['n_unique_roothogs']:,}")
    print(f"  Embedding dimensions: {stats['embedding_dim']}")
    print("\n  Columns:")
    for col in merged_df.columns[:13]:
        print(f"    - {col}")
    if len(merged_df.columns) > 13:
        print(f"    ... and {len(merged_df.columns) - 13} more (including ME:1 to ME:{stats['embedding_dim']})")

    print("\nDone!")


if __name__ == "__main__":
    main()
