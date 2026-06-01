# Chapter 2 — Protein Embeddings & Clustering

How well does ESM2 capture evolutionary (HOG) and functional (GO) structure? Core finding so far: larger models often do *worse* on broad clustering ("inverse scaling"). See [../../WORKSTREAMS.md](../../WORKSTREAMS.md) for the current threads.

## Layout

```
chapter2/
├── 01_compute_embeddings.py    # step 1: ESM2 embeddings
├── 02_build_dataset.py         # step 2: merge sequences + taxa + GO + HOG + emb
├── 03_clustering_analysis.py   # step 3: KMeans + silhouettes + UMAP
├── 04_pretrain_esm2.py         # train ESM2 from scratch (T3)
├── 04_subtree_hog_analysis.py  # sub-HOG coloring (T4)
├── 05_go_prediction.py         # supervised GO (T2)
├── 05_model_size_clustering.py # scale-vs-clustering sweep (T1)
├── 06_go_grid_search.py        # GO hyperparam sweep (T2)
├── 07_whitening_analysis.py    # whitening rescue (T6)
├── fetch_mammalia.py           # OMA ingest for focused dataset (T5)
├── hog_umap.py                 # extract + embed + plot 8 HOGs × 30 species (T5)
├── hog_umap_interactive.py     # plotly click-to-compare (T5)
├── run.py                      # config-driven runner for 01→03
├── config.json                 # default pipeline config
├── run_*.sh                    # pretraining sweep drivers
├── lab/                        # per-experiment narratives (NNN_name/entry.md)
├── runs/                       # timestamped runner outputs
├── results/                    # analysis outputs (whitening, etc.)
├── figures/                    # publishable figures
└── archive/                    # old scripts/analyses preserved for reference
```

## Running

See [../../docs/running.md](../../docs/running.md) for all commands. Quick start:

```bash
python run.py                     # full broad CAFA3 pipeline
python hog_umap.py plot           # regenerate focused UMAPs
python hog_umap_interactive.py    # regenerate interactive HTML
```

## Results

Per-experiment table: [lab/summary.md](lab/summary.md) — covers all 14+ entries with key findings.

Key figures live in `figures/` (broad CAFA3) and `lab/NNN_*/figures/` (per-experiment).
