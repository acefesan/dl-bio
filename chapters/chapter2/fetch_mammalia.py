#!/usr/bin/env python3
"""Download and prepare a HOG-balanced mammalian protein dataset from OMA.

Fetches protein lists, performs HOG-balanced sampling, retrieves sequences,
and merges into a single feather dataset for 30 mammalian species.

Usage:
    python fetch_mammalia.py                          # Run all phases
    python fetch_mammalia.py fetch-lists               # Phase 1 only
    python fetch_mammalia.py sample --target-n 100     # Phase 2 only
    python fetch_mammalia.py fetch-sequences            # Phase 3 only
    python fetch_mammalia.py merge                      # Phase 4 only
    python fetch_mammalia.py all --species HUMAN,MOUSE  # Test with 2 species
    python fetch_mammalia.py all --bulk-fasta oma-seqs.fa.gz  # Use local FASTA
"""

import argparse
import gzip
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUTPUT_BASE = PROJECT_ROOT / "assets" / "proteins" / "mammalia"

OMA_API = "https://omabrowser.org/api"

# 30 species across 14 mammalian orders
ALL_SPECIES_CODES = [
    "ORNAN", "TACAU",                          # Monotremata
    "MONDO", "SARHA",                          # Metatheria
    "LOXAF", "ECHTE",                          # Afrotheria
    "DASNO",                                   # Xenarthra
    "ERIEU",                                   # Eulipotyphla
    "RHIFE", "MYOLU",                          # Chiroptera
    "FELCA", "CANLF", "AILME",                 # Carnivora
    "HORSE",                                   # Perissodactyla
    "BOVIN", "PIGXX", "TURTR",                 # Artiodactyla
    "MANJA",                                   # Pholidota
    "HUMAN", "MACMU", "CALJA", "NOMLE", "MICMU",  # Primates
    "MOUSE", "RATNO", "CAVPO", "HETGA",        # Rodentia
    "RABIT",                                   # Lagomorpha
    "TUPBE",                                   # Scandentia
    "BALMU",                                   # Cetacea
]


# =============================================================================
# OMA API helper
# =============================================================================

def oma_get(url: str, params: dict | None = None, rate_limit: float = 1.0,
            max_retries: int = 5) -> requests.Response:
    """GET with retry on 429/5xx and exponential backoff."""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.ok:
                time.sleep(rate_limit)
                return r
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2 ** attempt * 5, 120)
                tqdm.write(f"  HTTP {r.status_code}, retrying in {wait}s (attempt {attempt + 1})")
                time.sleep(wait)
                continue
            r.raise_for_status()
        except requests.exceptions.Timeout:
            wait = min(2 ** attempt * 5, 120)
            tqdm.write(f"  Timeout, retrying in {wait}s (attempt {attempt + 1})")
            time.sleep(wait)
        except requests.exceptions.ConnectionError:
            wait = min(2 ** attempt * 5, 120)
            tqdm.write(f"  Connection error, retrying in {wait}s (attempt {attempt + 1})")
            time.sleep(wait)
    raise RuntimeError(f"Failed after {max_retries} retries: {url}")


# =============================================================================
# Phase 1: Fetch protein lists
# =============================================================================

def fetch_species_metadata(species_codes: list[str], rate_limit: float) -> dict:
    """Fetch metadata for each species from /api/genome/{code}/."""
    species_meta = {}
    for code in tqdm(species_codes, desc="Fetching species metadata"):
        r = oma_get(f"{OMA_API}/genome/{code}/", rate_limit=rate_limit)
        data = r.json()
        species_meta[code] = {
            "code": code,
            "species": data["species"],
            "taxon_id": data["taxon_id"],
            "nr_entries": data["nr_entries"],
            "lineage": data.get("lineage", [])[:10],
        }
    return species_meta


def fetch_protein_list(species_code: str, rate_limit: float) -> pd.DataFrame:
    """Paginate through /api/genome/{code}/proteins/ and return all proteins."""
    records = []
    page = 1
    per_page = 500
    total = None

    with tqdm(desc=f"  {species_code}", unit=" proteins") as pbar:
        while True:
            r = oma_get(
                f"{OMA_API}/genome/{species_code}/proteins/",
                params={"per_page": per_page, "page": page},
                rate_limit=rate_limit,
            )
            data = r.json()

            if total is None:
                total = int(r.headers.get("x-total-count", 0))
                pbar.total = total

            for entry in data:
                records.append({
                    "omaid": entry.get("omaid"),
                    "canonicalid": entry.get("canonicalid"),
                    "oma_hog_id": entry.get("oma_hog_id"),
                    "sequence_length": entry.get("sequence_length"),
                })

            pbar.update(len(data))

            if len(data) < per_page:
                break
            page += 1

    return pd.DataFrame(records)


def phase_fetch_lists(species_codes: list[str], output_dir: Path,
                      rate_limit: float) -> None:
    """Phase 1: fetch species metadata and protein lists."""
    lists_dir = output_dir / "protein_lists"
    lists_dir.mkdir(parents=True, exist_ok=True)

    # Species metadata
    species_json = output_dir / "species.json"
    if species_json.exists():
        print(f"Loading existing species metadata from {species_json.name}")
        with open(species_json) as f:
            species_meta = json.load(f)
        # Fetch any missing species
        missing = [c for c in species_codes if c not in species_meta]
        if missing:
            print(f"Fetching metadata for {len(missing)} new species...")
            new_meta = fetch_species_metadata(missing, rate_limit)
            species_meta.update(new_meta)
            with open(species_json, "w") as f:
                json.dump(species_meta, f, indent=2)
    else:
        print("Fetching species metadata...")
        species_meta = fetch_species_metadata(species_codes, rate_limit)
        with open(species_json, "w") as f:
            json.dump(species_meta, f, indent=2)

    # Protein lists (resume: skip species whose feather exists)
    print(f"\nFetching protein lists for {len(species_codes)} species...")
    for code in species_codes:
        feather_path = lists_dir / f"{code}.feather"
        if feather_path.exists():
            n = len(pd.read_feather(feather_path))
            print(f"  {code}: already have {n:,} proteins, skipping")
            continue

        df = fetch_protein_list(code, rate_limit)
        df.to_feather(feather_path)
        print(f"  {code}: saved {len(df):,} proteins")

    print("\nPhase 1 complete.")


# =============================================================================
# Phase 2: HOG-balanced sampling
# =============================================================================

def hog_balanced_sample(protein_df: pd.DataFrame, target_n: int,
                        seed: int = 42) -> pd.DataFrame:
    """Sample proteins with maximum HOG diversity.

    Round 1 (diversity): 1 protein per root HOG.
    Round 2 (fill): proportional sampling from largest HOGs.
    """
    rng = np.random.default_rng(seed)

    has_hog = protein_df[
        protein_df["oma_hog_id"].notna() & (protein_df["oma_hog_id"] != "")
    ].copy()
    has_hog["roothog_id"] = has_hog["oma_hog_id"].str.split(".").str[0]

    # Round 1: 1 protein per root HOG
    round1 = (
        has_hog.groupby("roothog_id")
        .sample(1, random_state=int(rng.integers(1e9)))
        .reset_index(drop=True)
    )

    if len(round1) >= target_n:
        return round1.sample(target_n, random_state=int(rng.integers(1e9)))

    # Round 2: fill from largest HOGs
    remaining = target_n - len(round1)
    already = set(round1["omaid"])
    pool = has_hog[~has_hog["omaid"].isin(already)]
    hog_sizes = pool.groupby("roothog_id").size().sort_values(ascending=False)

    fill_parts = []
    for hog_id in hog_sizes.index:
        if remaining <= 0:
            break
        hog_prots = pool[pool["roothog_id"] == hog_id]
        n_take = min(len(hog_prots),
                     max(1, remaining * len(hog_prots) // len(pool)),
                     remaining)
        if n_take > 0:
            fill_parts.append(
                hog_prots.sample(n_take, random_state=int(rng.integers(1e9)))
            )
            remaining -= n_take

    if fill_parts:
        sampled = pd.concat([round1, pd.concat(fill_parts)])
    else:
        sampled = round1

    return sampled


def phase_sample(species_codes: list[str], output_dir: Path,
                 target_n: int, seed: int) -> None:
    """Phase 2: HOG-balanced sampling across all species."""
    lists_dir = output_dir / "protein_lists"
    sampled_dir = output_dir / "sampled"
    sampled_dir.mkdir(parents=True, exist_ok=True)

    # Load all protein lists
    all_dfs = []
    for code in species_codes:
        feather_path = lists_dir / f"{code}.feather"
        if not feather_path.exists():
            print(f"ERROR: {feather_path} not found. Run fetch-lists first.")
            sys.exit(1)
        df = pd.read_feather(feather_path)
        df["species_code"] = code
        all_dfs.append(df)

    # Find effective target (min of target and smallest species HOG-bearing count)
    min_hog_count = float("inf")
    min_species = None
    for df in all_dfs:
        n_hog = df["oma_hog_id"].notna().sum()
        if n_hog < min_hog_count:
            min_hog_count = n_hog
            min_species = df["species_code"].iloc[0]

    effective_target = min(target_n, min_hog_count)
    print(f"Target per species: {target_n}")
    print(f"Smallest HOG-bearing count: {min_hog_count:,} ({min_species})")
    print(f"Effective target: {effective_target:,}")

    # Sample each species
    sampled_parts = []
    stats_per_species = {}
    for df in tqdm(all_dfs, desc="Sampling species"):
        code = df["species_code"].iloc[0]
        sampled = hog_balanced_sample(df, effective_target, seed=seed)
        sampled["species_code"] = code
        sampled_parts.append(sampled)
        stats_per_species[code] = {
            "total_proteins": len(df),
            "hog_bearing": int(df["oma_hog_id"].notna().sum()),
            "sampled": len(sampled),
            "unique_roothogs": int(sampled["roothog_id"].nunique()),
        }

    all_sampled = pd.concat(sampled_parts, ignore_index=True)
    all_sampled.to_feather(sampled_dir / "sampled_proteins.feather")

    stats = {
        "effective_target": effective_target,
        "requested_target": target_n,
        "seed": seed,
        "n_species": len(species_codes),
        "total_sampled": len(all_sampled),
        "unique_roothogs": int(all_sampled["roothog_id"].nunique()),
        "per_species": stats_per_species,
    }
    with open(sampled_dir / "sampling_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    print(f"\nSampled {len(all_sampled):,} proteins across {len(species_codes)} species")
    print(f"Unique root HOGs: {all_sampled['roothog_id'].nunique():,}")
    print(f"Phase 2 complete.")


# =============================================================================
# Phase 3: Fetch sequences
# =============================================================================

def fetch_sequences_bulk_api(omaids: list[str], rate_limit: float,
                             batch_size: int = 1000) -> dict[str, str]:
    """Try POST /api/protein/bulk_retrieve/ for sequences."""
    sequences = {}
    for i in tqdm(range(0, len(omaids), batch_size), desc="  Bulk API"):
        batch = omaids[i:i + batch_size]
        try:
            r = requests.post(
                f"{OMA_API}/protein/bulk_retrieve/",
                json={"ids": batch},
                timeout=120,
            )
            if r.ok:
                for entry in r.json():
                    omaid = entry.get("omaid")
                    seq = entry.get("sequence")
                    if omaid and seq:
                        sequences[omaid] = seq
            else:
                tqdm.write(f"  Bulk API returned {r.status_code}, falling back")
                return sequences  # partial results; caller handles fallback
        except Exception as e:
            tqdm.write(f"  Bulk API error: {e}")
            return sequences
        time.sleep(rate_limit)
    return sequences


def fetch_sequences_from_bulk_fasta(omaids: set[str],
                                    fasta_path: Path) -> dict[str, str]:
    """Stream-extract matching IDs from a pre-downloaded oma-seqs.fa.gz."""
    sequences = {}
    open_fn = gzip.open if fasta_path.suffix == ".gz" else open
    mode = "rt" if fasta_path.suffix == ".gz" else "r"

    print(f"  Scanning {fasta_path.name} for {len(omaids):,} IDs...")
    found = 0
    with open_fn(fasta_path, mode) as f:
        for record in tqdm(SeqIO.parse(f, "fasta"), desc="  Scanning FASTA",
                           unit=" seqs"):
            if record.id in omaids:
                sequences[record.id] = str(record.seq)
                found += 1
                if found >= len(omaids):
                    break
    print(f"  Found {found:,} / {len(omaids):,} sequences in bulk FASTA")
    return sequences


def fetch_sequences_individual(omaids: list[str],
                               rate_limit: float) -> dict[str, str]:
    """Last resort: individual GET /api/protein/{omaid}/ calls."""
    sequences = {}
    for omaid in tqdm(omaids, desc="  Individual API"):
        try:
            r = oma_get(f"{OMA_API}/protein/{omaid}/", rate_limit=rate_limit)
            data = r.json()
            seq = data.get("sequence")
            if seq:
                sequences[omaid] = seq
        except Exception as e:
            tqdm.write(f"  Error fetching {omaid}: {e}")
    return sequences


def write_fasta(sequences: dict[str, str], output_path: Path) -> None:
    """Write sequences dict to FASTA file."""
    records = [
        SeqRecord(Seq(seq), id=omaid, description="")
        for omaid, seq in sequences.items()
    ]
    with open(output_path, "w") as f:
        SeqIO.write(records, f, "fasta")


def phase_fetch_sequences(species_codes: list[str], output_dir: Path,
                          rate_limit: float,
                          bulk_fasta: Path | None = None) -> None:
    """Phase 3: fetch sequences for sampled proteins."""
    sampled_path = output_dir / "sampled" / "sampled_proteins.feather"
    if not sampled_path.exists():
        print("ERROR: sampled_proteins.feather not found. Run sample first.")
        sys.exit(1)

    seq_dir = output_dir / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)

    sampled = pd.read_feather(sampled_path)

    for code in species_codes:
        fasta_path = seq_dir / f"{code}.fasta"
        if fasta_path.exists():
            n = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
            expected = len(sampled[sampled["species_code"] == code])
            if n >= expected:
                print(f"  {code}: already have {n:,} sequences, skipping")
                continue

        species_ids = sampled[sampled["species_code"] == code]["omaid"].tolist()
        if not species_ids:
            continue

        print(f"\n{code}: fetching {len(species_ids):,} sequences...")

        # Strategy 1: bulk API
        sequences = fetch_sequences_bulk_api(species_ids, rate_limit)

        # Strategy 2: bulk FASTA file
        missing = [oid for oid in species_ids if oid not in sequences]
        if missing and bulk_fasta and bulk_fasta.exists():
            extra = fetch_sequences_from_bulk_fasta(set(missing), bulk_fasta)
            sequences.update(extra)

        # Strategy 3: individual calls for remaining
        missing = [oid for oid in species_ids if oid not in sequences]
        if missing:
            print(f"  {len(missing):,} still missing, fetching individually...")
            extra = fetch_sequences_individual(missing, rate_limit)
            sequences.update(extra)

        write_fasta(sequences, fasta_path)
        print(f"  {code}: saved {len(sequences):,} / {len(species_ids):,} sequences")

    print("\nPhase 3 complete.")


# =============================================================================
# Phase 4: Merge
# =============================================================================

def phase_merge(species_codes: list[str], output_dir: Path) -> None:
    """Phase 4: merge sampled proteins with sequences into final dataset."""
    sampled_path = output_dir / "sampled" / "sampled_proteins.feather"
    if not sampled_path.exists():
        print("ERROR: sampled_proteins.feather not found. Run sample first.")
        sys.exit(1)

    sampled = pd.read_feather(sampled_path)
    seq_dir = output_dir / "sequences"

    # Load species metadata for taxon_id and scientific_name
    species_json = output_dir / "species.json"
    species_meta = {}
    if species_json.exists():
        with open(species_json) as f:
            species_meta = json.load(f)

    # Load all sequences
    all_seqs = {}
    for code in species_codes:
        fasta_path = seq_dir / f"{code}.fasta"
        if not fasta_path.exists():
            print(f"WARNING: {fasta_path} not found, sequences will be missing for {code}")
            continue
        for record in SeqIO.parse(fasta_path, "fasta"):
            all_seqs[record.id] = str(record.seq)

    print(f"Loaded {len(all_seqs):,} sequences from FASTA files")

    # Add sequence column
    sampled["sequence"] = sampled["omaid"].map(all_seqs)
    n_with_seq = sampled["sequence"].notna().sum()
    print(f"Proteins with sequence: {n_with_seq:,} / {len(sampled):,}")

    # Add species metadata
    sampled["taxon_id"] = sampled["species_code"].map(
        {code: meta.get("taxon_id") for code, meta in species_meta.items()}
    )
    sampled["scientific_name"] = sampled["species_code"].map(
        {code: meta.get("species") for code, meta in species_meta.items()}
    )

    # Select and order final columns
    cols = ["omaid", "species_code", "taxon_id", "scientific_name",
            "oma_hog_id", "roothog_id", "sequence", "sequence_length"]
    for col in cols:
        if col not in sampled.columns:
            sampled[col] = None
    final = sampled[cols].copy()

    output_path = output_dir / "mammalia_dataset.feather"
    final.to_feather(output_path)

    print(f"\nSaved {len(final):,} proteins to {output_path}")
    print(f"  Species: {final['species_code'].nunique()}")
    print(f"  Unique root HOGs: {final['roothog_id'].nunique():,}")
    print(f"  With sequence: {final['sequence'].notna().sum():,}")
    print(f"  File size: {output_path.stat().st_size / 1024**2:.1f} MB")
    print("\nPhase 4 complete.")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare a HOG-balanced mammalian dataset from OMA"
    )
    parser.add_argument("phase", nargs="?", default="all",
                        choices=["fetch-lists", "sample", "fetch-sequences",
                                 "merge", "all"],
                        help="Pipeline phase to run (default: all)")
    parser.add_argument("--target-n", type=int, default=12000,
                        help="Target proteins per species (default: 12000)")
    parser.add_argument("--rate-limit", type=float, default=1.0,
                        help="Seconds between API calls (default: 1.0)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for sampling (default: 42)")
    parser.add_argument("--species", type=str, default=None,
                        help="Comma-separated species codes for testing (e.g. HUMAN,MOUSE)")
    parser.add_argument("--bulk-fasta", type=str, default=None,
                        help="Path to pre-downloaded oma-seqs.fa.gz")
    parser.add_argument("--output-dir", type=str, default=None,
                        help=f"Output directory (default: {OUTPUT_BASE})")
    args = parser.parse_args()

    species_codes = (
        args.species.split(",") if args.species
        else ALL_SPECIES_CODES
    )
    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_BASE
    output_dir.mkdir(parents=True, exist_ok=True)
    bulk_fasta = Path(args.bulk_fasta) if args.bulk_fasta else None

    print(f"Mammalia dataset pipeline")
    print(f"  Phase: {args.phase}")
    print(f"  Species: {len(species_codes)} ({', '.join(species_codes[:5])}{'...' if len(species_codes) > 5 else ''})")
    print(f"  Target per species: {args.target_n:,}")
    print(f"  Output: {output_dir}")
    print()

    phases = {
        "fetch-lists": lambda: phase_fetch_lists(species_codes, output_dir, args.rate_limit),
        "sample": lambda: phase_sample(species_codes, output_dir, args.target_n, args.seed),
        "fetch-sequences": lambda: phase_fetch_sequences(species_codes, output_dir, args.rate_limit, bulk_fasta),
        "merge": lambda: phase_merge(species_codes, output_dir),
    }

    if args.phase == "all":
        for name, fn in phases.items():
            print(f"\n{'=' * 60}")
            print(f"Phase: {name}")
            print(f"{'=' * 60}\n")
            fn()
    else:
        phases[args.phase]()

    print("\nDone!")


if __name__ == "__main__":
    main()
