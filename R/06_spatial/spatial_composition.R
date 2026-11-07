# Spatial composition differences (Methods: Spatial annotation & composition).
# Within each region:
#   - cell-class proportions relative to all annotated cellbins;
#   - neuronal-subtype proportions relative to all annotated neuronal cellbins.
# Rural (A1) vs Urban (A2) composition is compared with chi-square tests, and
# P values are FDR-corrected (Benjamini-Hochberg) across all region x level
# tests. These tests describe composition within the sampled sections; they do
# NOT establish population-level replication (one individual per habitat).
#
# Input: RCTD outputs from spatial_rctd.R
#   {region}_{group}.annotated.Rdata       -> main_celltype = major cell class
#   {region}_{group}.neuron_subtype.Rdata  -> main_celltype = Ex#/In# subtype
# region in {HPR, LSR, NCL, PVN}; group A1 = Rural, A2 = Urban.

suppressPackageStartupMessages(library(Seurat))

regions <- c("HPR", "LSR", "NCL", "PVN")
groups  <- c(A1 = "Rural", A2 = "Urban")

annot_file   <- function(region, group) sprintf("%s_%s.annotated.Rdata", region, group)
subtype_file <- function(region, group) sprintf("%s_%s.neuron_subtype.Rdata", region, group)

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

# annotated identities (drop NA / Unknown) from one object
identities <- function(obj) {
  ct <- as.character(obj@meta.data$main_celltype)
  ct[!is.na(ct) & ct != "Unknown"]
}

# Build a group x category count matrix, proportions, and a chi-square test.
composition_test <- function(counts_by_group) {
  cats <- sort(unique(unlist(lapply(counts_by_group, names))))
  mat <- sapply(counts_by_group, function(tb) {
    v <- as.numeric(tb[cats])
    v[is.na(v)] <- 0
    v
  })                                   # categories x groups
  rownames(mat) <- cats
  mat <- t(mat)                        # groups x categories
  props <- mat / rowSums(mat)
  test <- suppressWarnings(chisq.test(mat))
  list(mat = mat, props = props, cats = cats,
       chisq = unname(test$statistic), df = unname(test$parameter),
       p = test$p.value)
}

prop_rows <- list()
test_rows <- list()

run_level <- function(rg, level, file_fun) {
  counts_by_group <- lapply(names(groups), function(g) {
    obj <- load_obj(file_fun(rg, g))
    table(identities(obj))
  })
  names(counts_by_group) <- unname(groups)   # Rural / Urban
  res <- composition_test(counts_by_group)

  for (i in seq_along(res$cats)) {
    prop_rows[[length(prop_rows) + 1]] <<- data.frame(
      region = rg, level = level, category = res$cats[i],
      Rural_n = res$mat["Rural", i], Urban_n = res$mat["Urban", i],
      Rural_prop = round(res$props["Rural", i], 5),
      Urban_prop = round(res$props["Urban", i], 5),
      stringsAsFactors = FALSE)
  }
  test_rows[[length(test_rows) + 1]] <<- data.frame(
    region = rg, level = level,
    chisq = round(res$chisq, 4), df = res$df, p_value = res$p,
    stringsAsFactors = FALSE)
}

for (rg in regions) {
  run_level(rg, "major_class", annot_file)
  run_level(rg, "neuronal_subtype", subtype_file)
}

props <- do.call(rbind, prop_rows)
tests <- do.call(rbind, test_rows)
tests$p_FDR <- p.adjust(tests$p_value, method = "BH")   # FDR across all tests
row.names(props) <- NULL
row.names(tests) <- NULL

write_csv_bom(props[props$level == "major_class", ], "composition_major_class.csv")
write_csv_bom(props[props$level == "neuronal_subtype", ], "composition_neuronal_subtype.csv")
write_csv_bom(tests, "composition_chisq_FDR.csv")

cat("===== chi-square + FDR =====\n")
print(tests, row.names = FALSE, digits = 4)
cat("\nWrote composition_major_class.csv / composition_neuronal_subtype.csv / composition_chisq_FDR.csv\n")