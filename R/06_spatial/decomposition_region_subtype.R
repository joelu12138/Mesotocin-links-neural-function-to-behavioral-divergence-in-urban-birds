setwd("data")
suppressPackageStartupMessages(library(Seurat))

## Oaxaca–Blinder 分解：每个脑区的 A1/A2 差异，由哪种神经元亚型导致
## ΔF = Σ p̄(f_A2−f_A1) [状态] + Σ f̄(p_A2−p_A1) [组成]

load("irgsea_AUCell_neuron_subtype.Rdata")  # expr2

auc <- GetAssayData(expr2, assay = "AUCell", layer = "data")  # genesets x cells
genesets <- rownames(auc)                                       # "Neuron-NEW1" ...
regions  <- c("HPR", "LSR", "NCL", "PVN")
ct <- as.character(expr2$main_celltype)
grp <- expr2$group
reg <- expr2$region

# 亚型排序：Excitatory-1..10, Inhibitory-1..3
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
    # 该脑区无细胞的亚型 → 贡献记为 0
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
      脑区 = rg, 基因集 = gs,
      区域功能_A1 = round(FA1, 5), 区域功能_A2 = round(FA2, 5),
      总差异_dF = round(dF, 5),
      状态贡献 = round(state, 5), 组成贡献 = round(comp, 5),
      状态占比 = round(state / dF, 3), 组成占比 = round(comp / dF, 3),
      stringsAsFactors = FALSE)

    state_c <- pbar * (fA2 - fA1)
    comp_c  <- fbar * (pA2 - pA1)
    total_c <- pA2 * fA2 - pA1 * fA1   # Σ total_c = dF（精确）
    det[[paste(rg, gs)]] <- data.frame(
      脑区 = rg, 基因集 = gs, 亚型 = subtypes,
      p_A1 = round(pA1, 5), p_A2 = round(pA2, 5), dp = round(pA2 - pA1, 5),
      f_A1 = round(fA1, 5), f_A2 = round(fA2, 5), df = round(fA2 - fA1, 5),
      状态贡献_亚型 = round(state_c, 5),
      组成贡献_亚型 = round(comp_c, 5),
      总贡献_亚型 = round(total_c, 5),
      总贡献占比 = round(total_c / dF, 3),
      stringsAsFactors = FALSE)
  }
}
main <- do.call(rbind, main); det <- do.call(rbind, det)
write_csv_bom(main, "功能差异_状态vs组成_分解.csv")
write_csv_bom(det,  "功能差异_状态vs组成_分解_按亚型.csv")

# 摘要：每个脑区×基因集 的主导亚型（|总贡献|最大）
cat("=== 各脑区×基因集 主导亚型（|总贡献|最大） ===\n")
for (rg in regions) {
  for (gs in base::unique(det$基因集)) {
    s <- det[det$脑区 == rg & det$基因集 == gs, ]
    if (nrow(s) == 0) next
    dF <- sum(s$总贡献_亚型)
    if (abs(dF) < 1e-9) next
    i <- which.max(abs(s$总贡献_亚型))
    cat(sprintf("%-4s %-7s dF=%+.5f 主导亚型 %-15s 状态贡献=%+.5f 组成贡献=%+.5f 总贡献占比=%+.2f\n",
        rg, gs, dF, s$亚型[i],
        s$状态贡献_亚型[i], s$组成贡献_亚型[i], s$总贡献占比[i]))
  }
}
cat("\n已写出: 功能差异_状态vs组成_分解.csv / _按亚型.csv\n")
