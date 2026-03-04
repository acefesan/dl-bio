#!/usr/bin/env python3
"""Config-driven pipeline runner for chapter 2.

Reads a JSON config file and executes the selected pipeline steps
(embeddings, dataset building, clustering) as subprocesses.

Usage:
    python chapters/chapter2/run.py                           # runs config.json (all steps)
    python chapters/chapter2/run.py --config my_run.json      # custom config
    python chapters/chapter2/run.py --steps 3                 # only step 3
    python chapters/chapter2/run.py --steps 2,3               # steps 2 and 3
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHAPTER_DIR = Path(__file__).parent
VENV_PYTHON = PROJECT_ROOT / ".venv" / "bin" / "python"

SCRIPTS = {
    1: CHAPTER_DIR / "01_compute_embeddings.py",
    2: CHAPTER_DIR / "02_build_dataset.py",
    3: CHAPTER_DIR / "03_clustering_analysis.py",
}


def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return json.load(f)


def build_step1_args(cfg: dict, run_dir: Path) -> list[str]:
    """Build CLI args for 01_compute_embeddings.py."""
    emb = cfg["embeddings"]
    args = ["--batch-size", str(emb["batch_size"])]
    if emb.get("model"):
        args += ["--model", emb["model"]]
    if emb.get("output"):
        output_path = PROJECT_ROOT / emb["output"]
        args += ["--output", str(output_path)]
    if emb.get("checkpoint_dir"):
        args += ["--checkpoint-dir", str(PROJECT_ROOT / emb["checkpoint_dir"])]
    if emb.get("max_seq_length") is not None:
        args += ["--max-seq-length", str(emb["max_seq_length"])]
    if emb.get("limit") is not None:
        args += ["--limit", str(emb["limit"])]
    return args


def build_step2_args(cfg: dict, run_dir: Path) -> list[str]:
    """Build CLI args for 02_build_dataset.py."""
    ds = cfg["dataset"]
    emb = cfg.get("embeddings", {})
    output_dir = run_dir / "dataset"
    args = ["--output-dir", str(output_dir), "--format", ds.get("format", "feather")]
    if ds.get("skip_taxonomy"):
        args.append("--skip-taxonomy")
    if ds.get("fetch_hogs"):
        args.append("--fetch-hogs")
    if ds.get("hog_cache"):
        hog_path = PROJECT_ROOT / ds["hog_cache"]
        args += ["--hog-cache", str(hog_path)]
    return args


def build_step3_args(cfg: dict, run_dir: Path) -> list[str]:
    """Build CLI args for 03_clustering_analysis.py."""
    cl = cfg["clustering"]
    emb = cfg.get("embeddings", {})
    output_dir = run_dir / "clustering"
    args = [
        "--output-dir", str(output_dir),
        "--data", str(run_dir / "dataset" / "cafa3_annotations.feather"),
        "--sample-size", str(cl["sample_size"]),
        "--min-taxa", str(cl["min_taxa"]),
        "--min-hog", str(cl["min_hog"]),
        "--n-clusters", str(cl["n_clusters"]),
    ]
    if emb.get("output"):
        args += ["--embeddings", str(PROJECT_ROOT / emb["output"])]
    return args


def run_step(step: int, cfg: dict, run_dir: Path) -> bool:
    """Run a single pipeline step. Returns True on success."""
    script = SCRIPTS[step]
    if not script.exists():
        print(f"ERROR: Script not found: {script}")
        return False

    if step == 1:
        args = build_step1_args(cfg, run_dir)
    elif step == 2:
        args = build_step2_args(cfg, run_dir)
    elif step == 3:
        args = build_step3_args(cfg, run_dir)
    else:
        print(f"ERROR: Unknown step {step}")
        return False

    python = str(VENV_PYTHON) if VENV_PYTHON.exists() else sys.executable
    cmd = [python, str(script)] + args
    print(f"\n{'=' * 60}")
    print(f"STEP {step}: {script.name}")
    print(f"{'=' * 60}")
    print(f"Command: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Config-driven pipeline runner for chapter 2")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to JSON config file (default: config.json next to this script)")
    parser.add_argument("--steps", type=str, default=None,
                        help="Comma-separated step numbers to run (overrides config)")
    args = parser.parse_args()

    # Resolve config path
    if args.config:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = Path.cwd() / config_path
    else:
        config_path = CHAPTER_DIR / "config.json"

    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    cfg = load_config(config_path)
    name = cfg.get("name", "run")

    # Determine which steps to run
    if args.steps:
        steps = [int(s.strip()) for s in args.steps.split(",")]
    else:
        steps = cfg.get("steps", [1, 2, 3])

    for s in steps:
        if s not in SCRIPTS:
            print(f"ERROR: Invalid step {s}. Valid steps: {list(SCRIPTS.keys())}")
            sys.exit(1)

    # Create timestamped run directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = CHAPTER_DIR / "runs" / f"{name}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Copy config for provenance
    shutil.copy2(config_path, run_dir / "config.json")

    print(f"Pipeline: {name}")
    print(f"Config:   {config_path}")
    print(f"Steps:    {steps}")
    print(f"Run dir:  {run_dir}")

    # Execute steps
    results = {}
    for step in steps:
        ok = run_step(step, cfg, run_dir)
        results[step] = ok
        if not ok:
            print(f"\nSTEP {step} FAILED — aborting remaining steps.")
            break

    # Summary
    print(f"\n{'=' * 60}")
    print("PIPELINE SUMMARY")
    print(f"{'=' * 60}")
    print(f"Run directory: {run_dir}")
    for step, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  Step {step} ({SCRIPTS[step].name}): {status}")

    # Show output locations
    if 2 in results and results[2]:
        print(f"\n  Dataset output:    {run_dir / 'dataset'}")
    if 3 in results and results[3]:
        print(f"  Clustering output: {run_dir / 'clustering'}")
    if 1 in results and results[1]:
        print(f"  Embeddings:        (default shared location)")

    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
