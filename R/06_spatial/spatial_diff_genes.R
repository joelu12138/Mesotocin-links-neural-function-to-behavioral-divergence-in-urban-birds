# Regional differential expression on Stereo-seq bin100 (Methods: Regional DE).
# bin100 = aggregation of 100 x 100 DNA-nanoball positions. Bins with fewer than
# 200 detected genes or fewer than 1,000 UMIs are excluded. Within each brain
# region, Urban vs Rural expression is tested with Seurat's Wilcoxon rank-sum
# test, min.pct = 0.10; genes with |log2FC| > 0.59 and BH-adjusted P < 0.05 are
# retained. Only one individual per habitat was profiled, so these contrasts
# describe differences between the sampled brains (bins are NOT biological
# replicates of habitat).
#
# Input contract (confirmed):
#   - bin100 Seurat objects, one file per region x habitat
#   - region in {HPR, LSR, NCL, PVN}; group A1 = Rural, A2 = Urban
#   - QC columns: nFeature_Spatial (detected genes), dnbCount (UMIs)
# EDIT bin100_file() below to match your actual file names / extension.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

regions   <- c("HPR", "LSR", "NCL", "PVN")
NFEAT_MIN <- 200     # bins with <200 genes excluded
DNB_MIN   <- 1000    # bins with <1000 UMIs excluded

# A1 = Rural, A2 = Urban
bin100_file <- function(region, group) {
  sprintf("%s_%s.bin100.Rdata", region, group)   # <-- EDIT to your file naming
}

load_obj <- function(f) {
  e <- new.env()
  load(f, envir = e)
  get(ls(e)[1], envir = e)   # robust to the internal object name
}

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

deg_all <- list()
for (rg in regions) {
  a1 <- load_obj(bin100_file(rg, "A1"))   # Rural
  a2 <- load_obj(bin100_file(rg, "A2"))   # Urban

  # bin-level QC (Methods: >=200 genes and >=1000 UMIs)
  a1 <- subset(a1, subset = nFeature_Spatial >= NFEAT_MIN & dnbCount >= DNB_MIN)
  a2 <- subset(a2, subset = nFeature_Spatial >= NFEAT_MIN & dnbCount >= DNB_MIN)
  a1$group <- "Rural"
  a2$group <- "Urban"

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

  write_csv_bom(mk, sprintf("region_DEG_%s.csv", rg))
  deg_all[[rg]] <- mk
  cat(sprintf("%-4s Rural=%d Urban=%d bins -> DEG=%d\n",
              rg, ncol(a1), ncol(a2), nrow(mk)))
}

deg_all <- do.call(rbind, deg_all)
row.names(deg_all) <- NULL
write_csv_bom(deg_all, "region_DEG_all.csv")
cat("\nWrote region_DEG_<region>.csv and region_DEG_all.csv\n")