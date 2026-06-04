# Containerized GitHub Actions Runner

This setup runs the self-hosted GitHub Actions runner inside one Docker container.

Network shape:

```text
GitHub Actions control plane
        ^
        |
outbound HTTPS only
        |
Squid egress proxy allowlist
        ^
        |
runner container on an internal Docker network
```

There are no published runner ports. GitHub does not need inbound access to the runner. The runner makes outbound HTTPS connections to GitHub Actions and polls/receives jobs over that channel.

## Why the Proxy Exists

Docker can put the runner on an `internal` network, which blocks direct internet access. The runner then gets only one route out: the Squid proxy. Squid allows only GitHub/GitHub-Actions domains needed for:

- runner registration and job polling,
- checking out the repo,
- downloading official actions,
- sending logs/results back to GitHub,
- runner updates.

If a tool ignores `HTTPS_PROXY`, it should fail closed because the runner container has no direct external network route.

## Start

Create a one-hour runner registration token:

```bash
gh api \
  --method POST \
  repos/acefesan/dl-bio/actions/runners/registration-token \
  --jq .token
```

Then start the runner:

```bash
cd infra/github-runner
RUNNER_TOKEN="<token from command above>" docker compose up -d --build
```

Check logs:

```bash
docker logs -f dl-bio-github-runner
```

## Start on Reboot with systemd

The Compose file already uses `restart: unless-stopped`, so Docker can restart the existing containers after reboot. A small systemd unit makes that explicit and gives you a normal service name to start/stop.

Install the unit:

```bash
sudo cp /home/acefsan/src/dl_bio/infra/github-runner/dl-bio-github-runner.service \
  /etc/systemd/system/dl-bio-github-runner.service
sudo systemctl daemon-reload
sudo systemctl enable dl-bio-github-runner.service
```

Start or stop it:

```bash
sudo systemctl start dl-bio-github-runner.service
sudo systemctl stop dl-bio-github-runner.service
sudo systemctl status dl-bio-github-runner.service
```

First-time registration still needs a fresh one-hour `RUNNER_TOKEN`. After the runner has registered, normal container restarts do not need that token because the runner config remains in the existing container.

If you delete/recreate the runner container from scratch, register it again:

```bash
cd /home/acefsan/src/dl_bio/infra/github-runner
RUNNER_TOKEN="$(gh api \
  --method POST \
  repos/acefesan/dl-bio/actions/runners/registration-token \
  --jq .token)" \
docker compose up -d --build
```

The runner should appear in:

```text
https://github.com/acefesan/dl-bio/settings/actions/runners
```

## Labels

The runner registers with:

```text
self-hosted
linux
dl-bio-container
```

Use it from a workflow with:

```yaml
runs-on: [self-hosted, linux, dl-bio-container]
```

## Notes

- No Docker socket is mounted into the runner.
- No WSL home/source directory is mounted into the runner.
- The only persistent runner storage is the named Docker volume `runner_work`.
- The image preinstalls Python, Node, `ruff`, `nbconvert`, `markdownlint-cli2`, and `lychee` so CI does not need PyPI/npm egress at job time.
- If a workflow needs a new tool, bake it into the image rather than downloading it during the job.

GitHub's self-hosted runner docs say runners require outbound HTTPS on port 443 and list the domains needed by function:

https://docs.github.com/en/actions/reference/self-hosted-runners-reference#communication
