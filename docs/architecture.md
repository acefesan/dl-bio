# Architecture

## Repository layout

```
dl_bio/
├── README.md              # entry point
├── WORKSTREAMS.md         # threads of investigation + status
├── AGENTS.md              # pointers for agents
├── docs/                  # architecture, running, data
├── chapters/chapter2/     # all current work
│   ├── 0N_*.py           # numbered pipeline steps
│   ├── run.py            # config-driven runner for 01→03
│   ├── config.json       # default pipeline config
│   ├── fetch_mammalia.py # OMA data ingest (T5)
│   ├── hog_umap*.py      # focused UMAP + interactive viz (T5)
│   ├── lab/              # per-experiment notebooks (001..014)
│   ├── runs/             # timestamped runner outputs
│   ├── results/          # analysis outputs (whitening, etc.)
│   └── figures/          # publishable figures
├── dlfb/                  # submodule: dlfb library
├── notebooks/             # submodule: official book notebooks
└── assets/                # data (git-ignored)
```

## Chapter 2 pipelines

### Broad CAFA3 pipeline (T1, T2, T6)
`01_compute_embeddings.py` → `02_build_dataset.py` → `03_clustering_analysis.py`, orchestrated by `run.py` with `config.json`.

Inputs: `assets/proteins/datasets/train_sequences.fasta` (~141k CAFA3 proteins).
Output: `cafa3_with_embeddings.feather` (~1.8 GB) + clustering figures in `results/`.

Additional analyses consume the merged feather:
- `04_subtree_hog_analysis.py` — sub-HOG coloring (T4)
- `05_go_prediction.py` / `06_go_grid_search.py` — supervised GO (T2)
- `07_whitening_analysis.py` — scaling rescue (T6)

### Pretraining pipeline (T3)
`04_pretrain_esm2.py` / `pretrain_esm2.py` — train ESM2 from scratch at multiple scales. Drivers: `run_1m_sweep*.sh`, `run_200epoch_sweep.sh`, `run_60epoch_scaled_lr.sh`.

### Mammalian focused pipeline (T5)
`fetch_mammalia.py` — 4-phase OMA ingest (fetch-lists → sample → fetch-sequences → merge).
`hog_umap.py` — extract 8 root HOGs × 30 species, embed with 4 ESM2 scales, plot static UMAPs.
`hog_umap_interactive.py` — plotly HTML with click-to-compare HOG paths.

Inputs: OMA API + local `oma-seqs.fa.gz` (5.0 GB). Output: `mammalia_dataset.feather`.

## Data flow

```
raw sources ──▶ fetch scripts ──▶ feather caches ──▶ analysis scripts ──▶ figures/results
                                       │
                                       └──▶ lab entries (entry.md + metadata.json)
```

Lab entries (`chapters/chapter2/lab/NNN_name/`) record the narrative; `runs/` records the raw config+output of each `run.py` invocation. `lab/summary.md` is the cumulative results table.

## Model scales in use

| Model | Params | Dim | Purpose |
|---|---|---|---|
| esm2_t6_8M_UR50D | 8M | 320 | minimal baseline (T5) |
| esm2_t30_150M_UR50D | 150M | 640 | best broad clustering (T1) |
| esm2_t33_650M_UR50D | 650M | 1280 | best supervised GO (T2) |
| esm2_t36_3B_UR50D | 3B | 2560 | worst on most tasks (T1, T2) |
| ~1.47B from-scratch | 1.47B | 1920 | pretraining upper end (T3) |
