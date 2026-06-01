# Network and I/O Instrumentation

## Summary

When a long-running fetch script is suspiciously slow, the only honest answer is *measure what it is doing*. "It feels slow" and "elapsed time is high" are not measurements. The categories of evidence that actually distinguish causes are:

| Question | Measurement |
|---|---|
| Is the process making forward progress? | per-process CPU time, RSS, application-level log lines |
| Is the process receiving bytes? | per-process or per-interface bandwidth |
| Is the process making many small requests, or few big ones? | request count, average request size |
| Is the process blocked on the network, on disk, or on CPU? | process state code, syscall trace, stack sample |
| Is something else on the box competing for the network? | per-process bandwidth ranked across all processes |

This page is a tool-by-tool reference for getting those measurements without rerunning the workload from scratch. For the incident that motivated it, see [001 fetch stall post-mortem](../labs/001-fetch-stall-postmortem.md) and the underlying [raw evidence](../raw/lab-001-stall-postmortem.md). For the storage-layer reason a slow fetch can be CPU-busy while sending almost no bytes per second, see [TileDB-SOMA storage](tiledb-soma-storage.md).

## The Bandwidth Misreading That Almost Fooled Us

The clue that broke open Lab 001's first fetch was that a system monitor reported **very low downlink bandwidth** while the Python process was still using ~40% CPU and steadily growing in RSS. The intuitive reading is "nothing is downloading, the process is stuck". That reading is wrong in three common situations and a measurement plan has to distinguish them:

1. **Bursty I/O.** A SOMA fetch may pull a fragment, then spend many seconds decompressing and filtering it in CPU before the next fetch. Average bandwidth over the burst is high; instantaneous bandwidth between bursts is near zero.
2. **Latency-bound transfer.** Many small S3 requests with ~100 ms round-trip time saturate request-rate, not bandwidth. A monitor that samples bandwidth once per second sees a low number; the actual constraint is requests per second.
3. **Stalled syscall.** The process is genuinely blocked. CPU is idle, RSS does not grow, no syscalls fire. Distinguishable from cases 1 and 2 by *RSS staying flat* and CPU sitting near 0%.

The right tools answer each of these in turn.

## Per-Process Bandwidth

### `nethogs`

`nethogs` ranks running processes by current network usage, like `top` for the network. The headline columns are `SENT KB/s` and `RECEIVED KB/s` per PID.

```bash
sudo nethogs -d 2 eth0       # refresh every 2 seconds, watch eth0
sudo nethogs -t              # trace mode: human-readable lines for logging
```

For Lab 001 use: open in a second terminal while the fetch runs, sort by RECEIVED to confirm whether the Python process is the dominant consumer and at what level.

Limitations:
- needs root,
- WSL2 interfaces (`eth0`, `veth*`, `tailscale0`) may all be present; pick the one your default route uses,
- the per-process attribution can lag a few seconds after a process exits.

### `iftop`

`iftop` shows traffic by **connection**, not by process: source IP, destination IP, and current bandwidth.

```bash
sudo iftop -i eth0 -N        # don't resolve service names
sudo iftop -i eth0 -P        # show ports
```

Useful for confirming the destination of the traffic. For Lab 001 the brain fetch should produce many connections to `*.amazonaws.com` (specifically `s3-us-west-2`, or whichever endpoint serves `cellxgene-census-public-us-west-2`). If `iftop` shows only background hosts (system updates, telemetry), the fetch is not downloading.

### `bmon` and `nload`

Per-interface bandwidth time series. Less useful for attribution, more useful for confirming the link itself is not saturated by something invisible to `nethogs`.

```bash
bmon -p eth0
nload eth0
```

## Per-Process Disk I/O

### `iotop`

Ranks processes by disk read/write throughput. Even a network-fetch workload writes to disk when it serializes checkpoints; a checkpoint that should be happening but is not shows up as zero disk write.

```bash
sudo iotop -o                # only show processes doing I/O
sudo iotop -ao               # accumulated bytes since launch
```

For Lab 001: the fetch script checkpoints every 3 datasets. `iotop -ao` between checkpoints should show a clear write burst from the Python process. If those bursts are missing, the partial cache is not being written.

### `pidstat`

Periodic per-process I/O snapshots, scriptable.

```bash
pidstat -d -p <PID> 2        # disk I/O every 2 seconds
pidstat -u -p <PID> 2        # CPU
pidstat -r -p <PID> 2        # memory
```

## Live Socket Inventory

### `ss`

Lists current TCP/UDP sockets. The flag set that matters here:

```bash
ss -tnp                      # tcp, numeric, with process
ss -tnp '( dst :443 )'       # only HTTPS
ss -tnpie                    # add socket buffer / congestion info
```

A SOMA fetch should show one or more `ESTAB` connections to S3 endpoints attributed to the Python process. The `Recv-Q` column on those sockets tells you how many bytes are sitting in the kernel receive buffer waiting for the application to read. Persistently high `Recv-Q` means the kernel has data but the Python side is too busy (likely decompressing the previous chunk) to consume it. Persistently zero `Recv-Q` plus high CPU means the application is CPU-bound between fetches.

## Process Activity

### `ps` state code

The `STAT` column in `ps -p <PID>` is a quick triage.

| Code | Meaning | Implication |
|---|---|---|
| `R` | running on CPU | the process is doing work right now |
| `S` | interruptible sleep | usually waiting on I/O or a timer |
| `D` | uninterruptible sleep | waiting on disk I/O |
| `Z` | zombie | dead, parent has not reaped |
| `T` / `t` | stopped / traced | paused |
| trailing `N` | niced | running at lower priority (we use `nice -n 10`) |
| trailing `l` | multithreaded | normal for `cellxgene_census` |
| trailing `+` | foreground | normal for bash background tasks |

A long-running fetch will spend most of its life in `S` (waiting on the next S3 chunk) with occasional `R` bursts (decompression). Sustained `D` over many minutes suggests local disk is the bottleneck. Sustained `R` at 100% CPU suggests pure compute.

### `py-spy`

Sampling profiler for a *running* Python process. Does not require restarting the script.

```bash
py-spy dump --pid <PID>      # one-shot stack dump per thread
py-spy top --pid <PID>       # live "top" view of Python functions
py-spy record -o flame.svg --pid <PID> --duration 60   # flamegraph
```

For Lab 001: a `py-spy dump` of the stuck process would have told us within seconds whether time was being spent in `cellxgene_census.get_anndata`, in `tiledbsoma` C++ bindings, in zstd decompression, in `numpy` sparse construction, or in `urllib3` reading S3 bytes. Each of those points to a different fix.

### `strace`

Last-resort syscall trace. Heavy: it slows the target process.

```bash
sudo strace -p <PID> -c -f      # summary of syscalls over the run, with threads
sudo strace -p <PID> -e network # only network syscalls, live
```

If `py-spy` is blocked (rare), `strace` confirms whether `recvfrom`, `read`, or no syscalls are happening.

## Cumulative Counters

`/proc/net/dev` is the cheapest possible measurement: lifetime byte and packet counters per interface, no permissions needed.

```bash
cat /proc/net/dev
```

Sampling this twice with a known interval gives an average bandwidth over that interval, without any tool installed:

```bash
awk '/eth0/ {print $2}' /proc/net/dev; sleep 10; awk '/eth0/ {print $2}' /proc/net/dev
```

Useful when nothing else is installed (containers, fresh boxes) and as a sanity check on `nethogs`.

`/proc/<PID>/io` is the per-process equivalent for disk:

```bash
cat /proc/<PID>/io
```

The fields `read_bytes` and `write_bytes` are cumulative for that process.

## Putting It Together For The Next Fetch

A minimal instrumentation plan for the next ADORA fetch:

| Layer | Tool | What it confirms |
|---|---|---|
| application | log every chunk start, end, and byte count from inside the script | progress is real and at the expected rate |
| process | `pidstat -u -d -r -p <PID> 30` redirected to a file | CPU, disk, memory time series alongside the application log |
| network | `nethogs -t -d 10` redirected to a file | per-process bandwidth ranked across the whole box |
| sockets | one-shot `ss -tnp '( dst :443 )' > sockets.txt` during a known stall | what S3 endpoints we are still talking to |
| stack | one-shot `py-spy dump --pid <PID> > stack.txt` during a known stall | where time is actually being spent in Python and C++ |

Without these the next slow run is no more diagnosable than the last one. With them, "is it the network, the storage layout, or our code?" becomes an answerable question.

## What This Page Does Not Cover

- Application-level metrics from `tiledb` itself (stats counters from `tiledb.stats_enable() / tiledb.stats_dump()`). These exist for the raw TileDB library and are partly exposed through `tiledbsoma`; a future page should document them.
- Continuous monitoring stacks (Prometheus, Grafana) for long-running services. Out of scope for a single-script overnight run.

Related pages: [TileDB-SOMA storage](tiledb-soma-storage.md), [001 fetch stall post-mortem](../labs/001-fetch-stall-postmortem.md), [001 data flow](../labs/001-data-flow.md), [001 H5AD and AnnData cache](../labs/001-h5ad-anndata-cache.md)
