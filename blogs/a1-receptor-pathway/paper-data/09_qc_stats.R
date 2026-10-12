# Numbers behind each QC rule in 01_naive_b.R, for the playground's
# step-by-step explanations: per-cell distributions of genes detected and %
# mitochondrial reads (healthy donors, before QC), how many cells/genes each
# rule removed. Writes scgsa-site/data/provenance_qc.json. Run from paper-data/.
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(jsonlite) })
dir <- "/mnt/bulk/dl_bio/gse164381"
healthy <- c("GSM5008746_2010113866", "GSM5008747_2010113873", "GSM5008748_2010113875", "GSM5008749_2010113904")
raw <- do.call(cbind, lapply(healthy, function(s) Read10X_h5(file.path(dir, paste0(s, "_G.filtered_feature_bc_matrix.h5")))))
genes_total <- nrow(raw); cells_total <- ncol(raw)
nfeat <- colSums(raw > 0); ncount <- colSums(raw)
mt <- grepl("^MT-", rownames(raw))
pmt <- 100 * colSums(raw[mt, ]) / ncount
gene_cells <- rowSums(raw > 0)
drop_minfeat <- nfeat < 200
k <- !drop_minfeat
drop_mt <- k & pmt >= 15
drop_hi <- k & nfeat >= 4000
hist_of <- function(v, br) { h <- hist(pmin(v, max(br) - 1e-9), breaks = br, plot = FALSE); list(breaks = br, counts = h$counts) }
out <- list(
  genes_total = genes_total, cells_total = cells_total,
  genes_in_lt3_cells = sum(gene_cells < 3), genes_in_0_cells = sum(gene_cells == 0),
  cells_lt200_genes = sum(drop_minfeat),
  cells_mt_ge15 = sum(drop_mt), cells_genes_ge4000 = sum(drop_hi), cells_both = sum(drop_mt & drop_hi),
  cells_after = sum(k & !drop_mt & !drop_hi),
  median_genes = median(nfeat[k]), median_umi = median(ncount[k]), median_pmt = round(median(pmt[k]), 2),
  hist_genes = hist_of(nfeat[k], seq(0, 6000, by = 100)),
  hist_pmt = hist_of(pmt[k], seq(0, 60, by = 1))
)
write_json(out, file.path(getwd(), "..", "scgsa-site", "data", "provenance_qc.json"), auto_unbox = TRUE)
str(out[setdiff(names(out), c("hist_genes", "hist_pmt"))])
