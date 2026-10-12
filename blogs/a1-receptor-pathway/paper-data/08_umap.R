# UMAP of the healthy-donor cells from 01_naive_b.R, for the playground's
# "Build the dataset" node: shows the clusters, which were kept as naive B,
# which donor each cell came from, and the marker genes behind the choice.
# Same PCs as the clustering (1-15). Exports a 5,000-cell random sample.
# Writes scgsa-site/data/provenance_umap.json. Run from paper-data/.
suppressPackageStartupMessages({ library(Seurat); library(jsonlite) })
dir <- "/mnt/bulk/dl_bio/gse164381"
o <- readRDS(file.path(dir, "healthy_clustered.rds"))
o <- RunUMAP(o, dims = 1:15, seed.use = 42, verbose = FALSE)
keep <- c("0", "1", "3", "4", "7")
genes <- c("MS4A1", "IGHD", "TCL1A", "FCER2", "CD27", "IGHG1", "CD3E", "NKG7", "LYZ", "JCHAIN")
set.seed(42)
idx <- sort(sample(ncol(o), min(5000, ncol(o))))
um <- Embeddings(o, "umap")[idx, ]
ex <- as.matrix(LayerData(o, layer = "data")[genes, idx])
donors <- sort(unique(o$donor))
cl <- as.character(Idents(o))
tab <- table(cl, o$donor)
pct <- round(100 * prop.table(tab, 1))
out <- list(
  n_total = ncol(o), n_shown = length(idx),
  donors = paste("healthy donor", seq_along(donors)), genes = genes, kept = keep,
  x = round(um[, 1], 2), y = round(um[, 2], 2),
  cluster = as.integer(cl[idx]), donor = match(o$donor[idx], donors) - 1L,
  expr = setNames(lapply(genes, function(g) round(ex[g, ], 2)), genes),
  cluster_donor_pct = lapply(rownames(pct), function(k) list(cluster = k, cells = sum(tab[k, ]), pct = as.integer(pct[k, ])))
)
write_json(out, file.path(getwd(), "..", "scgsa-site", "data", "provenance_umap.json"), auto_unbox = TRUE, digits = NA)
cat("wrote provenance_umap.json:", length(idx), "cells\n")
