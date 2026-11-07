suppressMessages(library(Seurat))

# ---- QC thresholds (Methods: cellbin >50 genes, >500 UMIs, mt fraction <20%) ----
NFEAT_MIN <- 50
DNB_MIN   <- 500

# ---- area lookup: same bbox area from bin100 CSV ----
csv <- list.files("data", pattern = "_bin100\\.csv$", full.names = TRUE)
stopifnot(length(csv) == 1)
area <- read.csv(csv, stringsAsFactors = FALSE, skip = 1, check.names = FALSE)
area_map <- setNames(area$bbox_mm2, area$dataset)

files <- list.files("data", pattern = "\\.sub\\.Rdata$", full.names = TRUE)
subtype_all <- list()
qc_summary  <- list()

for (f in files) {
  e <- new.env()
  load(f, envir = e)
  o <- get(ls(e)[1], envir = e)
  md <- o@meta.data
  ds <- sub(".sub.Rdata", "", basename(f), fixed = TRUE)

  n_total <- nrow(md)
  # Methods: cellbin retained if >50 genes and >500 UMIs (mt<20% applied upstream)
  pass <- !is.na(md$main_celltype) & md$nFeature_Spatial > NFEAT_MIN & md$dnbCount > DNB_MIN

  qc_summary[[ds]] <- data.frame(
    dataset     = ds,
    n_total     = n_total,
    n_NA_subtype= sum(is.na(md$main_celltype)),
    n_pass      = sum(pass),
    n_dropped   = n_total - sum(pass),
    area_mm2_bbox = area_map[[ds]],
    stringsAsFactors = FALSE
  )

  tb <- as.data.frame(table(md$main_celltype[pass]), stringsAsFactors = FALSE)
  names(tb) <- c("subtype", "n_cells")
  tb$dataset <- ds
  tb$area_mm2 <- area_map[[ds]]
  tb$density_per_mm2 <- tb$n_cells / tb$area_mm2
  subtype_all[[ds]] <- tb

  cat(sprintf("DONE %-8s total=%d pass=%d dropped=%d area=%.4f mm2\n",
      ds, n_total, sum(pass), n_total - sum(pass), area_map[[ds]]))
}

qs <- do.call(rbind, qc_summary)
st <- do.call(rbind, subtype_all)
st <- st[order(st$dataset, -st$n_cells), ]
row.names(st) <- NULL

write_utf8 <- function(x, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE)
  write.csv(x, con, row.names = FALSE)
  close(con)
}
write_utf8(qs, "data/neuron_subtype_QC_summary.csv")
write_utf8(st, "data/neuron_subtype_density.csv")

cat("\n===== QC SUMMARY =====\n"); print(qs, row.names = FALSE)
cat("\n===== SUBTYPE DENSITY =====\n"); print(st, row.names = FALSE, digits = 3)
