# Raw: Lab 001 Fetch Stall Post-Mortem

## Provenance

This page summarizes the *observed* behavior of a single run of `lab/001_adora_expression/fetch_adora_cache.py` on the local workstation between **2026-05-31 00:06 PDT** and **2026-05-31 21:48 PDT**. The evidence is the log file the script wrote to disk plus live `ps`, `ss`, `df`, and `/proc/net/dev` snapshots taken before the process was terminated.

This is a source page (raw observations). Interpretation lives in [001 fetch stall post-mortem](../labs/001-fetch-stall-postmortem.md). Conceptual background lives in [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) and [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md).

## Source Files

- `lab/001_adora_expression/cache/fetch.log` (only file left in `cache/` after the run)
- `lab/001_adora_expression/fetch_adora_cache.py` (the script that produced the log)

The script was invoked as:

```bash
nice -n 10 python fetch_adora_cache.py --tissues brain --retries 5 --backoff 60
```

Run as foreground bash background-task ID `bwr9icldg`. Python PID `697993`.

## Log Excerpt

The full `fetch.log` for this run, verbatim:

```text
2026-05-31 00:06:45,124 INFO Cache dir: /home/acefsan/src/dl_bio/projects/caffeine/lab/001_adora_expression/cache
2026-05-31 00:06:45,124 INFO Tissues:   ['brain']
2026-05-31 00:06:45,520 INFO The "stable" release is currently 2025-11-08. Specify 'census_version="2025-11-08"' in future calls to open_soma() to ensure data consistency.
2026-05-31 07:08:10,159 INFO [brain] 186 datasets, 18539831 total cells
2026-05-31 07:08:10,180 INFO   [1/186] 00476f9f-ebc1-4b72-b541-32f912ce36ea
2026-05-31 18:22:26,957 INFO   [2/186] 0087cde2-967d-4f7c-8e6e-40e4c9ad1891
2026-05-31 21:29:10,526 INFO   [3/186] 03d38670-1444-4001-bc53-9936e61d9b20
```

## Timing Breakdown

| Phase | Wall-clock | Notes |
|---|---:|---|
| `open_soma()` to obs enumeration done | **7h 02m** | a single `obs.read(value_filter=...)` that returned 18,539,831 brain cells across 186 dataset_ids |
| dataset 1 fetched (`get_anndata`) | **11h 14m** | no retry warnings logged |
| dataset 2 fetched (`get_anndata`) | **3h 07m** | no retry warnings logged |
| dataset 3 started, killed before finishing | — | first checkpoint would have fired on completion of dataset 3 (i % 3 == 0) |

The 32-minute total of the script's own retry/backoff schedule (60, 120, 240, 480, 960 seconds) is far smaller than any of the per-dataset wall-clock figures above, and zero `attempt %d/%d failed` warnings appeared in the log, so the slow path is inside the `get_anndata` call itself, not in the script's exception handler.

## Process Snapshot Before Termination

Taken at 2026-05-31 21:48:27 PDT, ~21 hours into the run:

```text
PID     ELAPSED   %CPU %MEM   RSS    STAT  CMD
697993  21:03:35  39.2 1.1    955896 SNl   python fetch_adora_cache.py --tissues brain ...
```

| Field | Value | Reading |
|---|---|---|
| `STAT=SNl` | sleeping, niced, multithreaded | normal for an I/O-bound multithreaded Python process |
| `%CPU=39.2` | 21h average | sustained compute, not idle |
| `RSS=956 MB` | up from ~109 MB shortly after start | dominated by the two completed AnnData pieces held in the `pieces[]` list awaiting checkpoint |

## Disk Snapshot

```text
Filesystem      Size  Used Avail Use%
/dev/sdd       1007G  758G  199G  80%
```

`cache/` directory contained only `fetch.log` (636 bytes). No `.h5ad` file existed because the first checkpoint condition `i % 3 == 0 and pieces` had not yet fired.

## Network Snapshot

`ss -tnp` taken minutes after termination showed no remaining `*.amazonaws.com` connections from any user process. `/proc/net/dev` showed `eth0` lifetime `RX bytes = 50,265,341,506` (50.3 GB), which is **cumulative since boot**, not attributable to this run alone.

No other process on the box was visibly consuming sustained bandwidth in the snapshots taken during the run. The system has VS Code Remote, Codex extension, immich_ml, llama-swap, openclaw-gateway, and several other long-running services, but none appeared as a heavy network user when checked.

## Exit Disposition

The process was terminated by `kill -TERM 697993` and exited with code `143` (`128 + SIGTERM`). No core dump, no kernel OOM, no segfault, no on-disk checkpoint. All in-memory progress (two completed datasets) was lost.

## What This Page Does Not Establish

- A reproducible per-dataset throughput rate (we have only two completed datasets).
- The actual cause-of-slowness in the storage layer (we have circumstantial evidence consistent with cell-major fragment layout, not proof; see [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md) for the hypothesis).
- A historical bandwidth trace of the run. We never instrumented the process at the network or syscall level. The next run must, per [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md).

Related pages: [001 fetch stall post-mortem](../labs/001-fetch-stall-postmortem.md), [001 data flow](../labs/001-data-flow.md), [TileDB-SOMA storage](../concepts/tiledb-soma-storage.md), [network and I/O instrumentation](../concepts/network-and-io-instrumentation.md)
