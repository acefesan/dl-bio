# dl-bio

My work following [Deep Learning for Biology](https://github.com/deep-learning-for-biology). Focused on ESM2 protein language model embeddings — how well they capture evolutionary (HOG) and functional (GO) structure, and why larger models often do worse on broad clustering.

## Setup

```bash
git clone --recurse-submodules <repo-url>
uv sync
source .venv/bin/activate
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init
```

Submodules: `dlfb/` ([dlfb library](https://github.com/deep-learning-for-biology/dlfb)), `notebooks/` ([official notebooks](https://github.com/deep-learning-for-biology/notebooks)).

## Where to look

| You want to... | Read |
|---|---|
| see what's being investigated + status | [WORKSTREAMS.md](WORKSTREAMS.md) |
| see all experiments and results | [chapters/chapter2/lab/summary.md](chapters/chapter2/lab/summary.md) |
| understand the pipelines + repo layout | [docs/architecture.md](docs/architecture.md) |
| run a pipeline | [docs/running.md](docs/running.md) |
| find datasets and caches | [docs/data.md](docs/data.md) |
| work on chapter 2 specifically | [chapters/chapter2/README.md](chapters/chapter2/README.md) |

Agent context (for Claude Code / other coding agents): [AGENTS.md](AGENTS.md).

## Current focus

- **Area A — `chapters/chapter2/`:** ESM2 protein embeddings & clustering, following the DL for Biology book (threads T1–T6, active)
- **Area B — `projects/caffeine/`:** computational cross-tissue caffeine epigenome mapping, standalone (T7, proposal only; see [projects/caffeine/PROPOSAL.md](projects/caffeine/PROPOSAL.md))

See [WORKSTREAMS.md](WORKSTREAMS.md) for the live status of every thread.
