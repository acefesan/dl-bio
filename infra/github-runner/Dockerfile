FROM ubuntu:24.04

ARG RUNNER_VERSION=2.328.0
ARG LYCHEE_VERSION=0.20.1

ENV DEBIAN_FRONTEND=noninteractive
ENV RUNNER_ALLOW_RUNASROOT=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        jq \
        nodejs \
        npm \
        python3 \
        python3-pip \
        python3-venv \
        tar \
        unzip \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --break-system-packages --no-cache-dir \
        nbconvert \
        ruff

RUN npm install --global \
        markdownlint-cli2

RUN curl -fsSL "https://github.com/lycheeverse/lychee/releases/download/lychee-v${LYCHEE_VERSION}/lychee-x86_64-unknown-linux-gnu.tar.gz" \
        -o lychee.tar.gz \
    && tar -xzf lychee.tar.gz \
    && mv lychee /usr/local/bin/lychee \
    && chmod +x /usr/local/bin/lychee \
    && rm lychee.tar.gz

WORKDIR /runner

RUN curl -fsSL "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz" \
        -o actions-runner.tar.gz \
    && tar -xzf actions-runner.tar.gz \
    && rm actions-runner.tar.gz \
    && ./bin/installdependencies.sh

COPY entrypoint.sh /runner/entrypoint.sh
RUN chmod +x /runner/entrypoint.sh

ENTRYPOINT ["/runner/entrypoint.sh"]
