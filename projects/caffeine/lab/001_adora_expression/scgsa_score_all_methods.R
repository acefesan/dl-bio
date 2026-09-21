#!/usr/bin/env Rscript
# Run all seven scGSA methods from Wang & Thakar 2024 (NARGAB, lqae124) on the
# A1 gene panels, using the authors' own invocations wherever possible.
#
# Structure mirrors the paper's own comparison scripts (scPS-main/):
#   - irGSEA::irGSEA.score()  -> AUCell, JASMINE, ssGSEA, UCell, scSE
#   - Seurat::AddModuleScore() with defaults
#   - scPS() sourced from the authors' scPS_2024.R
#
# Usage: Rscript score_all_methods.R <atlas>   (atlas = hbca | tabula)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(GSEABase)
})

args <- commandArgs(trailingOnly = TRUE)
atlas <- if (length(args)) args[1] else "hbca"

TMP <- Sys.getenv("SCGSA_WORK", unset = "scgsa_work")
DATA <- file.path(TMP, "data")
OUT <- file.path(TMP, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ---------------------------------------------------------------- gene panels
# Minimal panel: receptor + neuronal transducer + the two effector arms that
# have hippocampus-specific literature, plus one canonically Gi-inhibited AC.
PANEL_MIN <- c("ADORA1", "GNAO1", "GNAI1", "KCNJ3", "KCNJ6",
               "CACNA1A", "CACNA1B", "ADCY5")

# Full 61-gene exhaustive panel, carried over unchanged for comparison.
PANEL_FULL <- c(
  "ADORA1", "GNAI1", "GNAI2", "GNAI3", "GNAO1",
  "ADCY1", "ADCY5", "ADCY6", "ADCY8",
  "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKACA", "PRKACB",
  "CREB1",
  "KCNJ3", "KCNJ6", "KCNJ9", "KCNJ5",
  "CACNA1B", "CACNA1A",
  "PLCB1", "PLCB2", "PLCB3", "PLCB4", "PRKCA", "PRKCB", "PRKCD", "PRKCE",
  "CALM1", "CALM2", "CALM3",
  "NFKB1", "RELA",
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PREX1", "RAC1",
  "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
  "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK14",
  "AKT1", "AKT2", "AKT3", "VEGFA",
  "NOS1", "NOS2", "NOS3", "GUCY1A1", "GUCY1B1", "PRKG1", "PRKG2"
)

PANELS <- list(minimal8 = PANEL_MIN, full61 = PANEL_FULL)

# ------------------------------------------------------------------ load data
cat("== loading", atlas, "==\n")
counts <- Read10X_h5(file.path(DATA, paste0(atlas, "_counts.h5")))
meta <- read.csv(file.path(DATA, paste0(atlas, "_cells.csv")),
                 stringsAsFactors = FALSE)
stopifnot(ncol(counts) == nrow(meta))
colnames(counts) <- meta$barcode
rownames(meta) <- meta$barcode

obj <- CreateSeuratObject(counts = counts, meta.data = meta[, "group",
                                                            drop = FALSE])
obj <- NormalizeData(obj, verbose = FALSE)
cat("  object:", nrow(obj), "genes x", ncol(obj), "cells,",
    length(unique(meta$group)), "groups\n")

present <- rownames(obj)
for (nm in names(PANELS)) {
  miss <- setdiff(PANELS[[nm]], present)
  PANELS[[nm]] <- intersect(PANELS[[nm]], present)
  cat("  panel", nm, ":", length(PANELS[[nm]]), "genes present",
      if (length(miss)) paste0("(missing: ", paste(miss, collapse = ","), ")")
      else "", "\n")
}

scores <- list()   # method -> panel -> named numeric vector over cells

# ------------------------------------------------------- 1. AddModuleScore
# Paper's script: AddModuleScore(obj, predBind) with all defaults
# (nbin = 24, ctrl = 100).
cat("== AddModuleScore ==\n")
for (nm in names(PANELS)) {
  o <- AddModuleScore(obj, features = list(PANELS[[nm]]), name = "AMS",
                      seed = 42)
  scores[["AddModuleScore"]][[nm]] <- setNames(o$AMS1, colnames(o))
}

# ------------------------------------------------------------------ 2. scPS
# Sourced from the authors' scPS_2024.R. Their RunPCA call hardcodes npcs = 10,
# which cannot run on a panel smaller than 11 genes, so npcs is clamped to
# length(genes) - 1 here. That is a deviation from the published code and is
# recorded in the results as scps_npcs.
cat("== scPS ==\n")
scps_npcs <- list()
run_scPS <- function(obj, genes) {
  npcs <- min(10, length(genes) - 1)
  o <- ScaleData(obj, features = genes, verbose = FALSE)
  o <- RunPCA(o, features = genes, npcs = npcs, weight.by.var = FALSE,
              verbose = FALSE)
  emb <- Embeddings(o, "pca")
  sdev <- Stdev(o, "pca")
  varExp <- (sdev^2) / sum(sdev^2)
  cumVar <- cumsum(sdev^2) / sum(sdev^2)
  idx <- which(cumVar > 0.5)
  maxCombPC <- if (length(idx) == 0) 1 else max(min(idx), 1)
  PC.X <- emb - min(emb)                       # authors shift by global min
  w <- sapply(seq_len(maxCombPC), function(i) PC.X[, i] * varExp[i])
  agg <- if (is.matrix(w) && ncol(w) > 1) sqrt(rowSums(w)) else sqrt(as.vector(w))
  dat <- GetAssayData(o, layer = "data")
  meanExpr <- colMeans(dat[rownames(dat) %in% genes, , drop = FALSE])
  list(score = agg * meanExpr[names(agg)], npcs = npcs, maxCombPC = maxCombPC)
}
for (nm in names(PANELS)) {
  r <- run_scPS(obj, PANELS[[nm]])
  scores[["scPS"]][[nm]] <- r$score
  scps_npcs[[nm]] <- sprintf("npcs=%d,maxCombPC=%d", r$npcs, r$maxCombPC)
  cat("   ", nm, scps_npcs[[nm]], "\n")
}

# ------------------------------------------- 3-7. irGSEA: AUCell/JASMINE/etc.
have_irgsea <- requireNamespace("irGSEA", quietly = TRUE)
cat("== irGSEA available:", have_irgsea, "==\n")

if (have_irgsea) {
  suppressPackageStartupMessages(library(irGSEA))
  for (nm in names(PANELS)) {
    gs <- list(PANELS[[nm]])
    names(gs) <- nm
    o <- irGSEA::irGSEA.score(
      object = obj, assay = "RNA", slot = "data",
      seeds = 42, custom = TRUE, geneset = gs, ncores = 8,
      method = c("AUCell", "JASMINE", "ssgsea", "UCell", "scSE"),
      kcdf = "Gaussian",
      ucell.MaxRank = NULL, aucell.MaxRank = NULL,
      JASMINE.method = "oddsratio")
    for (a in c("AUCell", "JASMINE", "ssgsea", "UCell", "scSE")) {
      if (a %in% Seurat::Assays(o)) {
        m <- GetAssayData(o, assay = a, layer = "data")
        if (all(m == 0)) m <- GetAssayData(o, assay = a, layer = "counts")
        lbl <- c(AUCell = "AUCell", JASMINE = "JASMINE", ssgsea = "ssGSEA",
                 UCell = "UCell", scSE = "SCSE")[[a]]
        scores[[lbl]][[nm]] <- setNames(as.numeric(m[1, ]), colnames(m))
      }
    }
  }
} else {
  cat("!! irGSEA missing - see score_fallback.R\n")
}

# --------------------------------------------------------------- aggregate
cat("== aggregating ==\n")
rows <- list()
for (meth in names(scores)) {
  for (nm in names(scores[[meth]])) {
    v <- scores[[meth]][[nm]]
    v <- v[meta$barcode]
    agg <- tapply(v, meta$group, mean)
    rows[[length(rows) + 1]] <- data.frame(
      atlas = atlas, method = meth, panel = nm,
      group = names(agg), score = as.numeric(agg),
      n_cells = as.numeric(table(meta$group)[names(agg)]),
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, paste0(atlas, "_group_scores.csv")),
          row.names = FALSE)

percell <- do.call(cbind, lapply(names(scores), function(meth) {
  do.call(cbind, lapply(names(scores[[meth]]), function(nm) {
    m <- matrix(scores[[meth]][[nm]][meta$barcode], ncol = 1,
                dimnames = list(meta$barcode, paste(meth, nm, sep = "__")))
    m
  }))
}))
write.csv(data.frame(barcode = meta$barcode, group = meta$group, percell),
          file.path(OUT, paste0(atlas, "_percell_scores.csv")), row.names = FALSE)

cat("methods run:", paste(names(scores), collapse = ", "), "\n")
cat("wrote", file.path(OUT, paste0(atlas, "_group_scores.csv")), "\n")
