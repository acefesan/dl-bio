# Local CELLxGENE Census mirror — progress notes

Status as of 2026-09-24: **mirror complete, benchmark done.** Written up so the
blog posts (caffeine/ADORA series, see [`../post.md`](../post.md) and
[`../../caffeine-adenosine-receptors/`](../../caffeine-adenosine-receptors/)) can
re-run their single-cell queries locally instead of over S3.

## What exists

Census release **2025-11-08** (current `stable`), SOMA store, mirrored to the
shared bulk drive:

    /mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma/

| Path under `soma/` | Size |
|---|---:|
| `census_data/homo_sapiens` | 1224 GB |
| `census_data/mus_musculus` | 281 GB |
| `census_data/callithrix_jacchus` (marmoset) | 18.5 GB |
| `census_data/macaca_mulatta` | 12.1 GB |
| `census_data/pan_troglodytes` | 1.3 GB |
| `census_spatial_sequencing/*` | 9.2 GB |
| **Total** | **~1.55 TB** |

The primates and spatial data were not requested; the first sync's exclude only
matched `census_data/mus_musculus/*`. Harmless (`/mnt/bulk` had ~14 TB free),
delete if the space is wanted.

Why the harmonized SOMA store and not per-dataset H5ADs: it has one gene axis,
one global `soma_joinid`, and `is_primary_data` dedup flags, so cross-dataset
queries need no manual joins. See
`projects/caffeine/wiki/concepts/census-source-h5ads.md` for the trade-offs.

## Benchmark (the question Lab 001 got stuck on)

Query: all primary human brain cells x {ADORA1, ADORA2A, ADORA2B, ADORA3}
(`tissue_general == 'brain' and is_primary_data == True`), via
`cellxgene_census.get_anndata` on the local store —
[`bench_local_brain_adora.py`](bench_local_brain_adora.py).

| | Result |
|---|---|
| Local mirror | **1380 s (~23 min)**, 28,967,109 cells x 4 genes |
| Remote S3 (Lab 001 v2 probe, 2026-06-01) | 10 h 39 min for **11 cells** |

Caveats: it is not instant. The query still walks a lot of fragments (about
37 GB logically read in the first ~7 min, mostly page cache) — local disk
removes the network cost, not the decompression/filter cost. One run, cold
cache, no repeat timing. The result was not saved to disk; the script only
prints the AnnData summary.

## Reproducing

    # env: system Python here is 3.14 and tiledbsoma has no wheel for it,
    # so use the repo's uv-managed env (uv installed to ~/.local/bin)
    cd ~/src/dl_bio && uv sync
    ./.venv/bin/python blogs/a1-receptor-pathway/census-mirror/bench_local_brain_adora.py

    # re-mirror (needs the AWS CLI: `pip install awscli` in any venv; anonymous)
    ./census-mirror/sync_human.sh     # ~1.2 TB, ~4 h at ~50 MiB/s
    ./census-mirror/sync_mouse.sh     # ~281 GB

Open directly:

    import tiledbsoma, cellxgene_census
    census = tiledbsoma.Collection.open("/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")

## Gotchas hit

- `aws s3 sync` exits **2** on this box: `/mnt/bulk` is NTFS and every file
  triggers an "unable to update the last modified time" warning (~325k of them).
  Data is fine; check the log for real errors (`grep -iE "error|failed"` minus the
  utime lines — none found).
- The first mouse sync degraded to ~0.5 MiB/s on an 85 GB object (network and
  disk were each fine at ~30 MB/s when tested alone). Killing and restarting it
  fixed it. Partial objects are not resumed by `s3 sync`; delete the
  `*.tdb.<random>` temp file before restarting.
- Foreground `du`/`find` over the whole tree takes minutes on NTFS; run in the
  background.

## Not done / next

- Re-run the ADORA analyses from
  [`projects/caffeine/lab/001_adora_expression`](../../../projects/caffeine/lab/001_adora_expression/)
  against the full local data instead of the 1000-cells-per-type stratified
  sample (`fetch_adora_cache.py`, v3), and save the fetched AnnData to the
  lab `cache/` so it is not another 23 min each time.
- Mouse data was requested as potentially useful for cross-species comparison;
  no analysis has used it yet.
- Human Cell Atlas / Tabula Sapiens H5ADs already in that `cache/` (about
  83 GB) were not touched.
