---
name: orchestrate
description: >-
  Director/sub-agent workflow for non-trivial work in this repo. Use when a task
  is large enough to delegate (multi-file change, data regeneration, an analysis
  pass) and should be done in isolation before landing on main. The main session
  acts as DIRECTOR: it specs the work, spins up a git worktree, delegates to a
  sub-agent that works only in that worktree, then reviews the diff and merges to
  main. Triggers: "have a sub-agent do X", "delegate", "orchestrate", "director",
  "work on a worktree and merge".
---

# Orchestrate — Director & Sub-Agents on Worktrees

The main session is the **director**. The director does NOT do the bulk work
itself; it decomposes, delegates to sub-agents, reviews, and merges. Sub-agents
**always** work inside a dedicated git worktree on a branch, never on `main`.
Only the director merges to `main`.

## When to use
- Multi-file changes, data/figure regeneration, analysis passes, anything you'd
  want reviewed before it lands.
- NOT for trivial one-liners, a single doc edit, or a quick read — the director
  just does those directly.

## The loop

### 1. Spec the task (director)
Write a crisp brief the sub-agent can execute without guessing:
- **Inputs** (exact file paths), **Outputs** (exact artifacts), **Acceptance
  criteria** (how the director will know it's correct — numbers, files exist,
  script exits 0), and **Out of scope**.

### 2. Create the worktree + branch (director)
```bash
REPO=/home/acefsan/src/dl_bio
WT=$REPO/.worktrees/<task-slug>           # .worktrees/ is gitignored
git -C "$REPO" worktree add "$WT" -b agent/<task-slug>
```

### 3. Wire in gitignored large data (director) — IMPORTANT
A fresh worktree does **not** contain untracked/ignored files. This repo
gitignores big data: `projects/caffeine/lab/*/cache/` (H5ADs, npz, feather),
`checkpoints/`, model caches. If the task needs them, **symlink** the real ones
into the worktree before delegating:
```bash
ln -s "$REPO/projects/caffeine/lab/001_adora_expression/cache" \
      "$WT/projects/caffeine/lab/001_adora_expression/cache"
```
The script then finds inputs and writes ignored outputs (npz/feather) straight to
the real cache; only tracked changes (code, `figures/*.png`) live in the worktree
and get merged.

### 4. Delegate (director → sub-agent)
Spawn a `general-purpose` sub-agent. Tell it explicitly:
- `cd "$WT"` and do all work there; **never touch `main` or the primary checkout's
  tracked files**.
- Follow the spec; report back key numbers + the commit SHA.
- `git add` only the intended files and `git commit` on the current branch.
- Use the repo venv `/home/acefsan/src/dl_bio/.venv/bin/python` for Python.

### 5. Review (director)
```bash
git -C "$WT" log --oneline main..HEAD
git -C "$WT" diff --stat main..HEAD
git -C "$WT" diff main..HEAD -- '*.py' '*.md'    # eyeball code/doc changes
```
Check the acceptance criteria. Confirm no stray files, no large data committed
(`git -C "$WT" diff --stat main..HEAD` should not list H5AD/npz/feather), diff
matches the spec. If wrong, send the sub-agent back with specifics (reuse the
same agent via SendMessage) — do not silently fix large work yourself.

### 6. Merge (director)
```bash
git -C "$REPO" merge --no-ff agent/<task-slug> -m "Merge agent/<task-slug>: <summary>"
```
Resolve any conflicts on `main`. Push only if the user asked.

### 7. Clean up (director)
```bash
git -C "$REPO" worktree remove "$WT"
git -C "$REPO" branch -d agent/<task-slug>
```

## Conventions
- Branches: `agent/<task-slug>`. Worktrees: `.worktrees/<task-slug>` (gitignored).
- One task per worktree. Keep briefs and commits scoped to that task.
- The director owns `main` and every merge. Sub-agents own only their branch.
- A sub-agent that hits a session/tool limit may return empty — the director
  checks the worktree state and either retries, re-delegates, or finishes it.
