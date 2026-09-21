#!/usr/bin/env Rscript
# Independent implementation of the five irGSEA-wrapped methods, calling the
# authoritative packages directly with the parameters irGSEA resolves to when
# the paper's arguments (ucell.MaxRank = NULL, aucell.MaxRank = NULL,
# kcdf = "Gaussian", JASMINE.method = "oddsratio", seeds = 42) are used:
#
#   AUCell  aucMaxRank = ceiling(0.05 * n_genes)   (AUCell's own default)
#   UCell   maxRank    = 1500                      (UCell's own default)
#   ssGSEA  GSVA ssgsea
#   JASMINE authors' JASMINE_V1_11October2021.r, oddsratio variant
#   SCSE    colSums(geneset) / colSums(all)        (Pont et al. 2019 formula)
#
# Run alongside score_all_methods.R as a cross-check on the wrapper.

suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(AUCell); library(UCell)
  library(GSVA)
})

args <- commandArgs(trailingOnly = TRUE)
atlas <- if (length(args)) args[1] else "hbca"

TMP <- Sys.getenv("SCGSA_WORK", unset = "scgsa_work")
DATA <- file.path(TMP, "data"); OUT <- file.path(TMP, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# jasmine_functions.R is the authors' JASMINE_V1_11October2021.r with the
# top-level example call stripped so sourcing only defines the functions.
source(file.path(TMP, "jasmine_functions.R"))

PANEL_MIN <- c("ADORA1", "GNAO1", "GNAI1", "KCNJ3", "KCNJ6",
               "CACNA1A", "CACNA1B", "ADCY5")
PANEL_FULL <- c(
  "ADORA1","GNAI1","GNAI2","GNAI3","GNAO1","ADCY1","ADCY5","ADCY6","ADCY8",
  "PRKAR1A","PRKAR1B","PRKAR2A","PRKAR2B","PRKACA","PRKACB","CREB1",
  "KCNJ3","KCNJ6","KCNJ9","KCNJ5","CACNA1B","CACNA1A",
  "PLCB1","PLCB2","PLCB3","PLCB4","PRKCA","PRKCB","PRKCD","PRKCE",
  "CALM1","CALM2","CALM3","NFKB1","RELA",
  "PIK3CA","PIK3CB","PIK3CD","PIK3R1","PREX1","RAC1",
  "MAPK1","MAPK3","MAP2K1","MAP2K2","MAPK8","MAPK9","MAPK10","MAPK11","MAPK14",
  "AKT1","AKT2","AKT3","VEGFA","NOS1","NOS2","NOS3","GUCY1A1","GUCY1B1",
  "PRKG1","PRKG2")

cat("== loading", atlas, "==\n")
counts <- Read10X_h5(file.path(DATA, paste0(atlas, "_counts.h5")))
meta <- read.csv(file.path(DATA, paste0(atlas, "_cells.csv")),
                 stringsAsFactors = FALSE)
colnames(counts) <- meta$barcode
obj <- CreateSeuratObject(counts = counts)
obj <- NormalizeData(obj, verbose = FALSE)
expr <- GetAssayData(obj, layer = "data")   # log-normalised, genes x cells
cat("  expr:", nrow(expr), "x", ncol(expr), "\n")

PANELS <- list(minimal8 = intersect(PANEL_MIN, rownames(expr)),
               full61   = intersect(PANEL_FULL, rownames(expr)))

cat("== AUCell (buildRankings once) ==\n")
rankings <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
aucMaxRank <- ceiling(0.05 * nrow(rankings))
cat("   aucMaxRank =", aucMaxRank, "of", nrow(rankings), "genes\n")

res <- list()
for (nm in names(PANELS)) {
  genes <- PANELS[[nm]]
  cat("== panel", nm, "(", length(genes), "genes ) ==\n")
  gs <- list(genes); names(gs) <- nm

  auc <- AUCell_calcAUC(gs, rankings, aucMaxRank = aucMaxRank, verbose = FALSE)
  res[[paste0("AUCell__", nm)]] <- getAUC(auc)[1, ]
  cat("   AUCell ok\n")

  uc <- ScoreSignatures_UCell(expr, features = gs, maxRank = 1500,
                              ncores = 8, name = "")
  res[[paste0("UCell__", nm)]] <- setNames(uc[, 1], rownames(uc))
  cat("   UCell ok\n")

  ss <- tryCatch({
    p <- GSVA::ssgseaParam(as.matrix(expr), gs)
    GSVA::gsva(p, verbose = FALSE)[1, ]
  }, error = function(e) { cat("   ssGSEA FAILED:", conditionMessage(e), "\n")
                           rep(NA_real_, ncol(expr)) })
  res[[paste0("ssGSEA__", nm)]] <- ss
  cat("   ssGSEA ok\n")

  # SCSE: Pont et al. 2019 - sum over gene set / sum over all genes, per cell
  scse <- Matrix::colSums(expr[genes, , drop = FALSE]) / Matrix::colSums(expr)
  res[[paste0("SCSE__", nm)]] <- scse
  cat("   SCSE ok\n")

  jas <- JASMINE(as.matrix(expr), genes, method = "oddsratio")
  res[[paste0("JASMINE__", nm)]] <- setNames(jas$JAS_Scores, jas$SampleID)
  cat("   JASMINE ok\n")
}

cat("== aggregating ==\n")
rows <- list()
for (k in names(res)) {
  parts <- strsplit(k, "__")[[1]]
  v <- res[[k]][meta$barcode]
  agg <- tapply(v, meta$group, mean)
  rows[[length(rows) + 1]] <- data.frame(
    atlas = atlas, method = parts[1], panel = parts[2],
    group = names(agg), score = as.numeric(agg),
    n_cells = as.numeric(table(meta$group)[names(agg)]),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
write.csv(out, file.path(OUT, paste0(atlas, "_direct_group_scores.csv")),
          row.names = FALSE)
cat("wrote", file.path(OUT, paste0(atlas, "_direct_group_scores.csv")), "\n")
