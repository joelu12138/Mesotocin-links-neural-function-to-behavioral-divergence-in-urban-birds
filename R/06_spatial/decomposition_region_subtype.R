setwd("data")
suppressPackageStartupMessages(library(Seurat))

## Oaxaca-Blinder decomposition: for the A1/A2 difference in each brain region, which neuron subtype drives it
## DeltaF = Sum p_bar(f_A2-f_A1) [state] + Sum f_bar(p_A2-p_A1) [composition]

load("irgsea_AUCell_neuron_subtype.Rdata")  # expr2

auc <- GetAssayData(expr2, assay = "AUCell", layer = "data")  # genesets x cells
genesets <- rownames(auc)                                       # "Neuron-NEW1" ...
regions  <- c("HPR", "LSR", "NCL", "PVN")
ct <- as.character(expr2$main_celltype)
grp <- expr2$group
reg <- expr2$region

# Subtype order: Excitatory-1..10, Inhibitory-1..3
ex  <- base::unique(ct[grepl("Excitatory", ct)])
inh <- base::unique(ct[grepl("Inhibitory", ct)])
ex  <- ex[order(as.integer(sub("Excitatory-", "", ex)))]
inh <- inh[order(as.integer(sub("Inhibitory-", "", inh)))]
subtypes <- c(ex, inh)

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

main <- list(); det <- list()
for (rg in regions) {
  for (g in genesets) {
    x  <- as.numeric(auc[g, ])
    fA1 <- sapply(subtypes, function(s) mean(x[reg == rg & ct == s & grp == "A1"]))
    fA2 <- sapply(subtypes, function(s) mean(x[reg == rg & ct == s & grp == "A2"]))
    nA1 <- sapply(subtypes, function(s) sum(reg == rg & ct == s & grp == "A1"))
    nA2 <- sapply(subtypes, function(s) sum(reg == rg & ct == s & grp == "A2"))
    # subtypes with no cells in this brain region -> contribution set to 0
    fA1[is.na(fA1)] <- 0; fA2[is.na(fA2)] <- 0
    pA1 <- nA1 / sum(nA1); pA2 <- nA2 / sum(nA2)
    pbar <- (pA1 + pA2) / 2
    fbar <- (fA1 + fA2) / 2

    state <- sum(pbar * (fA2 - fA1))
    comp  <- sum(fbar * (pA2 - pA1))
    FA1 <- sum(pA1 * fA1); FA2 <- sum(pA2 * fA2)
    dF  <- FA2 - FA1

    gs <- sub("Neuron-", "", g)
    main[[paste(rg, gs)]] <- data.frame(
      brain_region = rg, geneset = gs,
      region_function_A1 = round(FA1, 5), region_function_A2 = round(FA2, 5),
      total_difference_dF = round(dF, 5),
      state_contribution = round(state, 5), composition_contribution = round(comp, 5),
      state_proportion = round(state / dF, 3), composition_proportion = round(comp / dF, 3),
      stringsAsFactors = FALSE)

    state_c <- pbar * (fA2 - fA1)
    comp_c  <- fbar * (pA2 - pA1)
    total_c <- pA2 * fA2 - pA1 * fA1   # Sum total_c = dF (exact)
    det[[paste(rg, gs)]] <- data.frame(
      brain_region = rg, geneset = gs, subtype = subtypes,
      p_A1 = round(pA1, 5), p_A2 = round(pA2, 5), dp = round(pA2 - pA1, 5),
      f_A1 = round(fA1, 5), f_A2 = round(fA2, 5), df = round(fA2 - fA1, 5),
      state_contribution_subtype = round(state_c, 5),
      composition_contribution_subtype = round(comp_c, 5),
      total_contribution_subtype = round(total_c, 5),
      total_contribution_proportion = round(total_c / dF, 3),
      stringsAsFactors = FALSE)
  }
}
main <- do.call(rbind, main); det <- do.call(rbind, det)
write_csv_bom(main, "functional_difference_state_vs_composition_decomposition.csv")
write_csv_bom(det,  "functional_difference_state_vs_composition_decomposition_by_subtype.csv")

# Summary: dominant subtype (largest |total contribution|) for each brain region x geneset
cat("=== Dominant subtype (largest |total contribution|) for each brain region x geneset ===\n")
for (rg in regions) {
  for (gs in base::unique(det$geneset)) {
    s <- det[det$brain_region == rg & det$geneset == gs, ]
    if (nrow(s) == 0) next
    dF <- sum(s$total_contribution_subtype)
    if (abs(dF) < 1e-9) next
    i <- which.max(abs(s$total_contribution_subtype))
    cat(sprintf("%-4s %-7s dF=%+.5f dominant subtype %-15s state contribution=%+.5f composition contribution=%+.5f total contribution proportion=%+.2f\n",
        rg, gs, dF, s$subtype[i],
        s$state_contribution_subtype[i], s$composition_contribution_subtype[i], s$total_contribution_proportion[i]))
  }
}
cat("\nWrote: functional_difference_state_vs_composition_decomposition.csv / _by_subtype.csv\n")
