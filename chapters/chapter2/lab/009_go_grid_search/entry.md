# 009 — GO Prediction Classifier Grid Search

**Date:** 2026-03-08
**Models:** ESM2 150M, 650M, 3B
**Status:** complete

## Hypothesis
Entry 007 used a fixed hidden_dim=256, 2-layer MLP for all embeddings. The 3B input
(2560D) is 4× wider than 150M (640D) — maybe a 256-dim bottleneck is too severe for
the 3B embeddings. If classifier capacity is the bottleneck, scaling up hidden_dim
and/or depth should disproportionately help 3B.

## Setup
- **Script:** `chapters/chapter2/06_go_grid_search.py`
- **Grid:** hidden_dim ∈ {128, 256, 512, 1024} × num_layers ∈ {1, 2, 3} × 3 embeddings = 36 runs
- **Training:** 300 steps, batch_size=64, lr=1e-3, Adam, seed=42
- **Eval:** every 50 steps, final test metrics on held-out set

## Results

### Best per embedding model

| Model | Best config | Params | Test auPRC | Test auROC |
|-------|-------------|--------|------------|------------|
| 150M  | 2L-1024H    | 2.0M   | 0.2971     | 0.8694     |
| 650M  | 2L-1024H    | 2.7M   | **0.3027** | **0.8794** |
| 3B    | 1L-1024H    | 2.9M   | 0.1139     | 0.6503     |

### Key patterns

| Observation | Detail |
|-------------|--------|
| 3B ceiling  | Even best 3B (0.114 auPRC) is **worse than worst 150M** (0.120, 1L-128H) |
| Capacity helps 150M/650M | 150M: 0.120→0.297 (2.5×), 650M: 0.160→0.303 (1.9×) |
| Capacity barely helps 3B | 3B: 0.076→0.114 (1.5×), and most gain from hidden_dim not depth |
| Depth hurts 3B | For 3B, 1-layer consistently beats 2 and 3 layers at every hidden_dim |
| 2 layers optimal for 150M/650M | Diminishing returns at 3 layers |

### 3B: auPRC by config

| Hidden | 1L | 2L | 3L |
|--------|-----|-----|-----|
| 128    | 0.083 | 0.076 | 0.076 |
| 256    | 0.093 | 0.087 | 0.080 |
| 512    | 0.106 | 0.090 | 0.082 |
| 1024   | **0.114** | 0.101 | 0.092 |

Depth actively hurts 3B — deeper networks perform worse at every hidden_dim.

## Figures
- `figures/go_grid_search_heatmap.png` — auPRC heatmaps (hidden × layers per model)
- `figures/go_grid_search_scaling.png` — auPRC vs classifier params, colored by embedding

## Interpretation
**The 3B problem is NOT classifier capacity.** Even with a 2.9M-param classifier (5× the
original), 3B auPRC (0.114) remains below the weakest 150M config (0.120). The information
simply isn't there — or isn't accessible to a linear probe.

The fact that depth hurts 3B specifically suggests the 3B embeddings have a different
geometry: deeper networks may amplify noise in the high-dimensional but uninformative
3B space. This supports the **anisotropy hypothesis** (open question #1) — the 3B
embeddings may be concentrated in a narrow cone where additional nonlinear layers
can't extract useful directions.

650M remains the sweet spot — same architecture wins (2L-1024H) and achieves the
highest metrics across the entire grid.

## Next steps
- **PCA/whitening on 3B** before training — if anisotropy is the issue, this should help
- **Linear probe** (0 hidden layers) to establish a floor — is the 3B linear signal weaker?
- **Dropout/regularization sweep** — 3B may benefit from regularization given its high input dim
