# Running the pipelines

Activate the venv first: `source .venv/bin/activate`. All paths are relative to the repo root.

## Broad CAFA3 pipeline

### Full run via the config-driven runner

```bash
python chapters/chapter2/run.py                      # all steps from config.json
python chapters/chapter2/run.py --config my_run.json
python chapters/chapter2/run.py --steps 3            # clustering only
python chapters/chapter2/run.py --steps 2,3          # dataset + clustering
```

Each invocation writes a timestamped dir under `chapters/chapter2/runs/` containing the config used and per-step outputs.

### Individual steps

```bash
# Step 1 — compute ESM2 embeddings (supports --resume)
python chapters/chapter2/01_compute_embeddings.py
python chapters/chapter2/01_compute_embeddings.py --limit 100     # smoke test

# Step 2 — merge sequences + taxonomy + GO + HOG + embeddings
python chapters/chapter2/02_build_dataset.py \
  --hog-cache assets/proteins/datasets/cafa3_merged/hog_cache.csv

# Step 3 — KMeans + silhouettes + UMAP
python chapters/chapter2/03_clustering_analysis.py
```

Flags: `--skip-taxonomy` skips UniProt calls; `--fetch-hogs` re-fetches HOGs from OMA.

### Other analyses on the merged feather

```bash
python chapters/chapter2/04_subtree_hog_analysis.py    # sub-HOG coloring (T4)
python chapters/chapter2/05_go_prediction.py           # supervised GO (T2)
python chapters/chapter2/06_go_grid_search.py          # GO hyperparam sweep (T2)
python chapters/chapter2/07_whitening_analysis.py      # whitening rescue (T6)
```

## Pretraining pipeline (T3)

```bash
python chapters/chapter2/04_pretrain_esm2.py   # single run
bash   chapters/chapter2/run_1m_sweep.sh       # multi-scale sweep
bash   chapters/chapter2/run_200epoch_sweep.sh # long-training rerun
```

## Mammalian focused pipeline (T5)

```bash
# 1. Download from OMA (resumable; cached under assets/)
python chapters/chapter2/fetch_mammalia.py fetch-lists
python chapters/chapter2/fetch_mammalia.py sample --target-n 12000
python chapters/chapter2/fetch_mammalia.py fetch-sequences
python chapters/chapter2/fetch_mammalia.py merge

# 2. Extract 8 HOGs, embed with 4 models, plot
python chapters/chapter2/hog_umap.py extract
python chapters/chapter2/hog_umap.py embed
python chapters/chapter2/hog_umap.py plot

# 3. Interactive HTML with click-to-compare
python chapters/chapter2/hog_umap_interactive.py
```

## Writing up a new experiment

1. Create `chapters/chapter2/lab/NNN_short_name/` with `entry.md` + `metadata.json`.
2. Record: hypothesis, setup (script + config), results table, interpretation, next steps.
3. Put figures under `lab/NNN_short_name/figures/`.
4. Add a row to `chapters/chapter2/lab/summary.md`.
5. If the experiment opens or closes a thread, update `WORKSTREAMS.md`.
