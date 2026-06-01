# Agent context — dl_bio

Pointers for coding agents. Keep this file short; put content in the documents it points to.

## Orientation

- **Status across threads:** [WORKSTREAMS.md](WORKSTREAMS.md) — start here for "what are we doing"
- **Per-experiment results:** [chapters/chapter2/lab/summary.md](chapters/chapter2/lab/summary.md)
- **Pipelines + layout:** [docs/architecture.md](docs/architecture.md)
- **How to run things:** [docs/running.md](docs/running.md)
- **Data locations:** [docs/data.md](docs/data.md)

## Environment

- Package manager: `uv` (this is a `pyproject.toml` project — use `uv add <pkg>`, not `pip install`)
- Activate venv: `source .venv/bin/activate`
- All current work is under `chapters/chapter2/`

## Conventions

- New experiments get a lab entry under `chapters/chapter2/lab/NNN_short_name/` with `entry.md` + `metadata.json`. Add a row to `lab/summary.md`. If a thread opens or closes, update `WORKSTREAMS.md`.
- Scripts are numbered by pipeline step (`01_*.py`, `02_*.py`, ...). The config-driven runner is `chapters/chapter2/run.py`.
- Timestamped runs go under `chapters/chapter2/runs/`. Narrative goes in `lab/`. Don't mix them.
- Figures for a lab entry live under `lab/NNN_.../figures/`. Publishable figures promoted to `chapters/chapter2/figures/`.
