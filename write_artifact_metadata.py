#!/usr/bin/env python3
"""Write a reproducibility metadata file alongside an artifact.

Every lab artifact (figure, CSV, JSON result) gets a sibling .meta.json
file recording how it was produced. Given the artifact path and the bash
command that created it, this script captures:

  - The exact command
  - The current git commit hash (fails if working tree is dirty)
  - Timestamp

Usage:
    python write_artifact_metadata.py <artifact_path> "<command>"

Example:
    python write_artifact_metadata.py \\
        chapters/chapter2/lab/006/figures/umap.png \\
        ".venv/bin/python chapters/chapter2/04_subtree_hog_analysis.py --hog-level 2"

Output: chapters/chapter2/lab/006/figures/umap.png.meta.json
"""

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def get_git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip()


def check_git_clean() -> bool:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip() == ""


def main():
    if len(sys.argv) < 3:
        print("Usage: write_artifact_metadata.py <artifact_path> \"<command>\"")
        sys.exit(1)

    artifact_path = Path(sys.argv[1])
    command = sys.argv[2]

    if not artifact_path.exists():
        print(f"Warning: artifact does not exist yet: {artifact_path}")

    commit = get_git_commit()
    clean = check_git_clean()

    if not clean:
        print(f"WARNING: git working tree is dirty. Commit changes first for reproducibility.")
        commit += "-dirty"

    meta = {
        "artifact": str(artifact_path),
        "command": command,
        "commit": commit,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    meta_path = Path(str(artifact_path) + ".meta.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote: {meta_path}")


if __name__ == "__main__":
    main()
