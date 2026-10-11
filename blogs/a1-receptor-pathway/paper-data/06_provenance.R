# Record how the playground's datasets were built, step by step, for the
# "Build the dataset" node: GEO samples -> healthy donors -> QC -> clusters
# (with the marker table used to pick naive B) -> kept cells; and the Splatter
# parameters fitted for SSD. Writes scgsa-site/data/provenance.json.
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(jsonlite); library(splatter) })
dir <- "/mnt/bulk/dl_bio/gse164381"
out <- file.path(getwd(), "..", "scgsa-site", "data", "provenance.json")  # run from paper-data/
samples <- sort(Sys.glob(file.path(dir, "*_G.filtered_feature_bc_matrix.h5")))
geo <- lapply(samples, function(f) {
  m <- Read10X_h5(f)
  id <- sub("_.*", "", basename(f))
  list(gsm = id, cells = ncol(m), healthy = id %in% c("GSM5008746", "GSM5008747", "GSM5008748", "GSM5008749"))
})
o <- readRDS(file.path(dir, "healthy_clustered.rds"))
nb <- readRDS(file.path(dir, "naive_b.rds"))
keep <- c("0", "1", "3", "4", "7")
naive <- c("IGHD", "TCL1A", "FCER2", "IL4R", "BACH2", "IGHM")
memory <- c("CD27", "IGHG1", "IGHA1", "AIM2", "TNFRSF13B", "CD80")
show <- c("MS4A1", "IGHD", "TCL1A", "CD27", "IGHG1", "CD3E", "NKG7", "JCHAIN")
avg <- as.matrix(AverageExpression(o, features = unique(c(naive, memory, show)), layer = "data")$RNA)
cl <- sub("^g", "", colnames(avg))
n <- as.integer(table(Idents(o))[cl])
score <- colMeans(avg[naive, , drop = FALSE]) - colMeans(avg[memory, , drop = FALSE])
clusters <- lapply(seq_along(cl), function(i) c(list(cluster = cl[i], cells = n[i],
  naive_minus_memory = round(unname(score[i]), 2), kept = cl[i] %in% keep),
  as.list(round(avg[show, i], 2))))
clusters <- clusters[order(-sapply(clusters, `[[`, "naive_minus_memory"))]
p <- readRDS(file.path(dir, "naive_b_splat_params.rds"))
prov <- list(
  geo = list(accession = "GSE164381", samples = geo, cells = sum(sapply(geo, `[[`, "cells"))),
  healthy_cells = sum(sapply(Filter(function(x) x$healthy, geo), `[[`, "cells")),
  qc = list(rule = "min.cells = 3, min.features = 200; percent.mt < 15; nFeature_RNA < 4000",
            cells_after = ncol(o)),
  clustering = list(steps = "LogNormalize; 2,000 variable genes; ScaleData; PCA (20); neighbours on PCs 1-15; Louvain resolution 0.5",
                    n_clusters = length(cl), markers_shown = show,
                    naive_markers = naive, memory_markers = memory, clusters = clusters),
  naive_b = list(cells = ncol(nb), genes = nrow(nb), donors = as.list(table(nb$donor))),
  splatter = list(nGenes = getParam(p, "nGenes"), mean_shape = round(getParam(p, "mean.shape"), 3),
                  mean_rate = round(getParam(p, "mean.rate"), 3), lib_loc = round(getParam(p, "lib.loc"), 3),
                  lib_scale = round(getParam(p, "lib.scale"), 3), bcv_common = round(getParam(p, "bcv.common"), 3),
                  bcv_df = round(getParam(p, "bcv.df"), 2), out_prob = round(getParam(p, "out.prob"), 4),
                  pool_cells = 4000)
)
write_json(prov, out, auto_unbox = TRUE, pretty = TRUE)
cat("wrote", normalizePath(out), "\n")
