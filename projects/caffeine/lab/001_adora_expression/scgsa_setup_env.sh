#!/usr/bin/env bash
# Build the R environment for the Wang & Thakar 2024 scGSA method comparison.
#
# There is no R in this repo's normal toolchain, and the seven methods are R
# packages, so this stands up a self-contained conda env. It is ~4 GB and
# entirely disposable: delete ~/.local/share/mamba to remove it.
#
# Usage:  bash scgsa_setup_env.sh
# Then:   SCGSA_WORK=<workdir> $R_BIN/Rscript scgsa_score_all_methods.R hbca
set -euo pipefail

MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-$HOME/.local/share/mamba}"
export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
ENV_NAME=scgsa
WORK="${SCGSA_WORK:-scgsa_work}"

mkdir -p "$MAMBA_ROOT/bin" "$WORK"

if [ ! -x "$MAMBA_ROOT/bin/micromamba" ]; then
  echo "== installing micromamba =="
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$MAMBA_ROOT" bin/micromamba
fi
MM="$MAMBA_ROOT/bin/micromamba"

echo "== creating env '$ENV_NAME' =="
"$MM" create -y -r "$MAMBA_ROOT" -n "$ENV_NAME" \
  -c conda-forge -c bioconda \
  r-base r-matrix r-seurat r-remotes r-data.table r-ggplot2 r-hdf5r \
  bioconductor-aucell bioconductor-ucell bioconductor-gsva \
  bioconductor-complexheatmap bioconductor-decoupler bioconductor-ggtree \
  bioconductor-nebulosa bioconductor-singscore r-gghalves

R_BIN="$MAMBA_ROOT/envs/$ENV_NAME/bin"

echo "== fetching the authors' own source =="
# scPS (the paper's method) and its comparison scripts
curl -sL -o "$WORK/scPS.tar.gz" \
  https://github.com/Thakar-Lab/scPS/archive/refs/heads/main.tar.gz
tar xzf "$WORK/scPS.tar.gz" -C "$WORK"

# JASMINE, from its own authors. The published file executes an example call
# at the top level, so strip it and keep only the function definitions.
curl -sL -o "$WORK/jasmine.tar.gz" \
  https://github.com/NNoureen/JASMINE/archive/refs/heads/main.tar.gz
tar xzf "$WORK/jasmine.tar.gz" -C "$WORK"
grep -vE "^Result[[:space:]]*=|^library\(GSA\)|^library\(stringr\)" \
  "$WORK/JASMINE-main/JASMINE_V1_11October2021.r" > "$WORK/jasmine_functions.R"

echo "== installing irGSEA (GitHub; the wrapper the paper used) =="
"$R_BIN/Rscript" -e '
  lib <- .libPaths()[1]; options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
  remotes::install_github("chuiqin/irGSEA", lib = lib, upgrade = "never")
  cat("irGSEA available:", requireNamespace("irGSEA", quietly = TRUE), "\n")'

echo
echo "done. R at: $R_BIN/Rscript"
echo "next:  SCGSA_WORK=$WORK python scgsa_extract_subsample.py"
echo "       SCGSA_WORK=$WORK $R_BIN/Rscript scgsa_score_all_methods.R hbca"
