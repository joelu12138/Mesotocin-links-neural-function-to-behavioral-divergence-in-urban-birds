# Regional differential expression restricted to annotated neurons (cellbin).
# Methods: cellbins with segmented nuclei and their transcripts were annotated
# with RCTD (spatial_rctd.R, full mode, major cell-class level). Here neurons
# (main_celltype == "Neuron") are extracted within each brain region, and
# Urban vs Rural expression is tested within that same region with Seurat's
# Wilcoxon rank-sum test. Thresholds are identical to spatial_diff_genes.R:
# min.pct = 0.10, |log2FC| > 0.59, BH-adjusted P < 0.05.
#
# Cellbin QC (>50 detected genes, >500 UMIs, mitochondrial fraction < 20%) is
# applied when the cellbin matrix is produced upstream; the objects read here
# are the RCTD-annotated outputs. Only one individual per habitat was
# profiled, so these contrasts describe differences between the sampled brains
# (cellbins are NOT biological replicates of habitat).
#
# Input: RCTD outputs from spatial_rctd.R
#   {region}_{group}.annotated.Rdata -> main_celltype = major cell class
# region in {HPR, LSR, NCL, PVN}; group A1 = Rural, A2 = Urban.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

regions <- c("HPR", "LSR", "NCL", "PVN")
groups  <- c(A1 = "Rural", A2 = "Urban")

annot_file <- function(region, group) sprintf("%s_%s.annotated.Rdata", region, group)

load_obj <- function(f) {
  e <- new.env()
  load(f, envir = e)
  get(ls(e)[1], envir = e)
}

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

# annotated neurons of one region x habitat, tagged with habitat
neurons_of <- function(region, group_code) {
  obj <- load_obj(annot_file(region, group_code))
  obj <- subset(obj, main_celltype == "Neuron")
  obj$group <- groups[[group_code]]
  obj
}

deg_all <- list()
for (rg in regions) {
  a1 <- neurons_of(rg, "A1")   # Rural
  a2 <- neurons_of(rg, "A2")   # Urban
  if (ncol(a1) < 2 || ncol(a2) < 2) {
    cat(sprintf("%-4s skipped (Rural=%d, Urban=%d neurons)\n", rg, ncol(a1), ncol(a2)))
    next
  }

  st <- merge(a1, a2)
  st <- JoinLayers(st)
  DefaultAssay(st) <- "Spatial"
  st <- NormalizeData(st, verbose = FALSE)

  mk <- FindMarkers(
    st, ident.1 = "Urban", ident.2 = "Rural", group.by = "group",
    assay = "Spatial", test.use = "wilcox",
    min.pct = 0.10, logfc.threshold = 0.59, verbose = FALSE)
  mk$gene <- rownames(mk)
  mk <- mk[mk$p_val_adj < 0.05 & abs(mk$avg_log2FC) > 0.59, ]
  mk$region <- rg

  write_csv_bom(mk, sprintf("region_neuron_DEG_%s.csv", rg))
  deg_all[[rg]] <- mk
  cat(sprintf("%-4s Rural=%d Urban=%d neurons -> DEG=%d\n",
              rg, ncol(a1), ncol(a2), nrow(mk)))
}

deg_all <- do.call(rbind, deg_all)
row.names(deg_all) <- NULL
write_csv_bom(deg_all, "region_neuron_DEG_all.csv")
cat("\nWrote region_neuron_DEG_<region>.csv and region_neuron_DEG_all.csv\n")
