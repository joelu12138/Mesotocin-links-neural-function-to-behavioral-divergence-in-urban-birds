## ============================================================
## 03_state_vs_composition_decomposition.R  (core)
## Question: is the brain region functional difference (DeltaF = F_HA - F_LA) caused by
##       1) intrinsic functional change of cell types ("state"), or
##       2) change in cell composition/proportion ("composition")?
## Method: Oaxaca-Blinder mean decomposition (exact identity, no residual).
##
##   Total region function  F^g = Sum_c p_c^g * f_c^g
##   p_c^g : proportion of cell type c in group g
##   f_c^g : mean AUCell of this pathway for cell type c in group g
##
##   DeltaF = Sum_c p_bar_c (f_HA - f_LA)   <- state contribution (the cells themselves changed function)
##      + Sum_c f_bar_c (p_HA - p_LA)   <- composition contribution (proportions changed)
##
## Interpretation: state proportion ~ 1 -> intrinsic cell function is the main driver;
##                 composition proportion ~ 1 -> cell proportions are the main driver;
##                 intermediate value -> mixed, report both proportions.
##
## Input: *.sub.Rdata (Seurat object A2_sub, Spatial assay)
##        marker.csv (genesets), KEGG_pathways_categories_WGCNA.csv (focused pathways)
## Output: functional_difference_state_vs_composition_decomposition.csv
##         functional_difference_state_vs_composition_decomposition_by_subtype.csv
## ============================================================

suppressMessages(library(irGSEA))
suppressMessages(library(UCell))

regions  <- c("CMM","HVC","LAMN","NCM","RA")
groups   <- c("HA","LA")
subtypes <- paste0("Ex_", 1:6)

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

# ---- 1. Read focused pathways and genesets ----
focus <- read.csv("KEGG_pathways_categories_WGCNA.csv", stringsAsFactors = FALSE)
pw    <- focus$Pathway                                    # focused pathways (63 entries)
df <- read.csv("marker.csv", header = TRUE, stringsAsFactors = FALSE)
gene_list <- split(df$SYMBOL, df$names)
miss <- setdiff(pw, names(gene_list))                     # pathways without a matching geneset (e.g. RNA transport)
if (length(miss)) cat("Skipping pathways without genesets:", paste(miss, collapse = ", "), "\n")
pw <- intersect(pw, names(gene_list))
gene_list <- gene_list[pw]

# ---- 2. Score each file (AUCell) and aggregate to subtype ----
bases   <- as.vector(outer(regions, groups, paste, sep = "_"))
present <- bases[file.exists(paste0(bases, ".sub.Rdata"))]
missing <- setdiff(bases, present)
if (length(missing)) cat("Skipping missing files:", paste(missing, collapse = ", "), "\n\n")

score_file <- function(f) {
  e <- new.env(); load(paste0(f, ".sub.Rdata"), envir = e)
  obj <- get("A2_sub", envir = e)
  expr2 <- irGSEA.score(object = obj, assay = "Spatial", slot = "data",
                        seeds = 123, ncores = 1, min.cells = 3, min.feature = 0,
                        custom = TRUE, geneset = gene_list, species = "Homo sapiens",
                        method = "AUCell", aucell.MaxRank = NULL,
                        ucell.MaxRank = NULL, kcdf = "Gaussian")
  sc <- as.matrix(expr2@assays$AUCell@data)                 # pathway x spot
  ct <- factor(obj@meta.data$main_celltype, levels = subtypes)
  list(scores = sc, ct = ct)
}

cache <- list()
for (f in present) {
  s <- score_file(f)
  keep <- !is.na(s$ct)                                      # drop Unknown
  sc   <- s$scores[, keep, drop = FALSE]
  ctv  <- s$ct[keep]
  fmat <- t(apply(sc, 1, function(r) tapply(r, ctv, mean))) # pathway x subtype means = f_c^g
  pvec <- as.numeric(table(ctv) / length(ctv))              # subtype proportions   = p_c^g
  cache[[f]] <- list(f = fmat, p = pvec, n = as.numeric(table(ctv)))
  cat("scored", f, "spots =", sum(keep), "\n")
}

# ---- 3. Oaxaca-Blinder decomposition ----
main <- list(); det <- list()
for (rg in regions) {
  ha <- paste0(rg, "_HA"); la <- paste0(rg, "_LA")
  if (!(ha %in% names(cache) && la %in% names(cache))) next
  HA <- cache[[ha]]; LA <- cache[[la]]
  for (k in seq_along(pw)) {
    fHA <- HA$f[k, ]; fLA <- LA$f[k, ]
    pHA <- HA$p;      pLA <- LA$p
    pbar <- (pHA + pLA) / 2
    fbar <- (fHA + fLA) / 2
    state <- sum(pbar * (fHA - fLA))   # state contribution
    comp  <- sum(fbar * (pHA - pLA))   # composition contribution
    FHA   <- sum(pHA * fHA); FLA <- sum(pLA * fLA)
    dF    <- FHA - FLA                 # total difference in region function

    main[[paste(rg, k)]] <- data.frame(
      brain_region = rg, pathway = pw[k],
      region_function_HA = round(FHA, 4), region_function_LA = round(FLA, 4),
      total_difference_dF = round(dF, 4),
      state_contribution = round(state, 4), composition_contribution = round(comp, 4),
      state_proportion = round(state / dF, 3), composition_proportion = round(comp / dF, 3),
      stringsAsFactors = FALSE)

    # ---- contribution of each subtype ----
    state_c <- pbar * (fHA - fLA)     # state contribution of this subtype (function changed)
    comp_c  <- fbar * (pHA - pLA)     # composition contribution of this subtype (proportion changed)
    total_c <- pHA * fHA - pLA * fLA  # contribution of this subtype to the total difference (Sum_c total_c = dF)
    det[[paste(rg, k)]] <- data.frame(
      brain_region = rg, pathway = pw[k], subtype = subtypes,
      p_HA = round(pHA, 4), p_LA = round(pLA, 4), dp = round(pHA - pLA, 4),
      f_HA = round(fHA, 4), f_LA = round(fLA, 4), df = round(fHA - fLA, 4),
      state_contribution_subtype = round(state_c, 4),
      composition_contribution_subtype = round(comp_c, 4),
      total_contribution_subtype = round(total_c, 4),
      total_contribution_proportion = round(total_c / dF, 3),
      stringsAsFactors = FALSE)
  }
}

if (length(main)) {
  main <- do.call(rbind, main); det <- do.call(rbind, det)
  write_csv_bom(main, "functional_difference_state_vs_composition_decomposition.csv")
  write_csv_bom(det,  "functional_difference_state_vs_composition_decomposition_by_subtype.csv")
  cat("\n=== Decomposition results (state proportion ~ 1 mostly state / ~ 0 mostly composition) ===\n")
  print(main, row.names = FALSE)
  cat("\nWrote: functional_difference_state_vs_composition_decomposition.csv / _by_subtype.csv\n")

  # ---- Summary: subtype with the largest contribution for each pathway (by absolute total contribution proportion) ----
  cat("\n=== Dominant subtype for each pathway (largest |total contribution proportion|) ===\n")
  for (rg in unique(det$brain_region)) {
    for (pw in unique(det$pathway[det$brain_region == rg])) {
      s <- det[det$brain_region == rg & det$pathway == pw, ]
      i <- which.max(abs(s$total_contribution_proportion))
      cat(sprintf("%-6s %-45s dominant subtype %s: state contribution=%.4f composition contribution=%.4f total contribution proportion=%.3f\n",
          rg, pw, s$subtype[i], s$state_contribution_subtype[i], s$composition_contribution_subtype[i],
          s$total_contribution_proportion[i]))
    }
  }
}
