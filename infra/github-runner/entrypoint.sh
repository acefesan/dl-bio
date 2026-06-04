#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${RUNNER_REPO_URL:-}" ]]; then
  echo "RUNNER_REPO_URL is required, e.g. https://github.com/acefesan/dl-bio" >&2
  exit 1
fi

RUNNER_NAME="${RUNNER_NAME:-dl-bio-container-$(hostname)}"
RUNNER_LABELS="${RUNNER_LABELS:-dl-bio-container,linux}"
RUNNER_WORKDIR="${RUNNER_WORKDIR:-/runner/_work}"

cleanup() {
  if [[ -f .runner && -n "${RUNNER_TOKEN:-}" ]]; then
    ./config.sh remove --unattended --token "${RUNNER_TOKEN}" || true
  fi
}
trap cleanup EXIT INT TERM

if [[ ! -f .runner ]]; then
  if [[ -z "${RUNNER_TOKEN:-}" ]]; then
    echo "RUNNER_TOKEN is required for first-time registration." >&2
    echo "After registration, normal container restarts do not need RUNNER_TOKEN." >&2
    exit 1
  fi

  ./config.sh \
    --unattended \
    --url "${RUNNER_REPO_URL}" \
    --token "${RUNNER_TOKEN}" \
    --name "${RUNNER_NAME}" \
    --labels "${RUNNER_LABELS}" \
    --work "${RUNNER_WORKDIR}" \
    --replace
fi

exec ./run.sh
