## ============================================================
## 03_状态vs组成_分解.R  (核心)
## 问题：脑区功能差异(ΔF = F_HA − F_LA)究竟是
##       ① 细胞类型内在功能改变（"状态"），还是
##       ② 细胞组成/比例改变（"组成"）导致？
## 方法：Oaxaca–Blinder 均值分解（精确恒等式，无残差）。
##
##   区域总功能  F^g = Σ_c p_c^g · f_c^g
##   p_c^g : 分组 g 中细胞类型 c 的比例
##   f_c^g : 分组 g 中细胞类型 c 的该通路 AUCell 均值
##
##   ΔF = Σ_c p̄_c (f_HA − f_LA)   ← 状态贡献（细胞本身功能变了）
##      + Σ_c f̄_c (p_HA − p_LA)   ← 组成贡献（比例变了）
##
## 判读：状态占比≈1 → 细胞内在功能是主因；
##       组成占比≈1 → 细胞比例是主因；
##       中间值   → 混合，报各自占比。
##
## 输入：*.sub.Rdata（Seurat 对象 A2_sub，Spatial assay）
##       marker.csv（基因集）、KEGG_pathways_categories_WGCNA.csv（重点通路）
## 输出：功能差异_状态vs组成_分解.csv
##       功能差异_状态vs组成_分解_按亚型.csv
## ============================================================

suppressMessages(library(irGSEA))
suppressMessages(library(UCell))

regions  <- c("CMM","HVC","LAMN","NCM","RA")
groups   <- c("HA","LA")
subtypes <- paste0("Ex_", 1:6)

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("﻿", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

# ---- 1. 读入重点通路与基因集 ----
focus <- read.csv("KEGG_pathways_categories_WGCNA.csv", stringsAsFactors = FALSE)
pw    <- focus$Pathway                                    # 关注通路（63 条）
df <- read.csv("marker.csv", header = TRUE, stringsAsFactors = FALSE)
gene_list <- split(df$SYMBOL, df$names)
miss <- setdiff(pw, names(gene_list))                     # 无对应基因集的通路（如 RNA transport）
if (length(miss)) cat("跳过无基因集通路:", paste(miss, collapse = ", "), "\n")
pw <- intersect(pw, names(gene_list))
gene_list <- gene_list[pw]

# ---- 2. 对每个文件打分(AUCell) 并聚合到亚型 ----
bases   <- as.vector(outer(regions, groups, paste, sep = "_"))
present <- bases[file.exists(paste0(bases, ".sub.Rdata"))]
missing <- setdiff(bases, present)
if (length(missing)) cat("跳过缺失文件:", paste(missing, collapse = ", "), "\n\n")

score_file <- function(f) {
  e <- new.env(); load(paste0(f, ".sub.Rdata"), envir = e)
  obj <- get("A2_sub", envir = e)
  expr2 <- irGSEA.score(object = obj, assay = "Spatial", slot = "data",
                        seeds = 123, ncores = 1, min.cells = 3, min.feature = 0,
                        custom = TRUE, geneset = gene_list, species = "Homo sapiens",
                        method = "AUCell", aucell.MaxRank = NULL,
                        ucell.MaxRank = NULL, kcdf = "Gaussian")
  sc <- as.matrix(expr2@assays$AUCell@data)                 # 通路 × spot
  ct <- factor(obj@meta.data$main_celltype, levels = subtypes)
  list(scores = sc, ct = ct)
}

cache <- list()
for (f in present) {
  s <- score_file(f)
  keep <- !is.na(s$ct)                                      # 去掉 Unknown
  sc   <- s$scores[, keep, drop = FALSE]
  ctv  <- s$ct[keep]
  fmat <- t(apply(sc, 1, function(r) tapply(r, ctv, mean))) # 通路 × 亚型 均值 = f_c^g
  pvec <- as.numeric(table(ctv) / length(ctv))              # 亚型比例   = p_c^g
  cache[[f]] <- list(f = fmat, p = pvec, n = as.numeric(table(ctv)))
  cat("scored", f, "spots =", sum(keep), "\n")
}

# ---- 3. Oaxaca–Blinder 分解 ----
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
    state <- sum(pbar * (fHA - fLA))   # 状态贡献
    comp  <- sum(fbar * (pHA - pLA))   # 组成贡献
    FHA   <- sum(pHA * fHA); FLA <- sum(pLA * fLA)
    dF    <- FHA - FLA                 # 区域功能总差异

    main[[paste(rg, k)]] <- data.frame(
      脑区 = rg, 通路 = pw[k],
      区域功能_HA = round(FHA, 4), 区域功能_LA = round(FLA, 4),
      总差异_dF = round(dF, 4),
      状态贡献 = round(state, 4), 组成贡献 = round(comp, 4),
      状态占比 = round(state / dF, 3), 组成占比 = round(comp / dF, 3),
      stringsAsFactors = FALSE)

    # ---- 每个亚型的贡献度 ----
    state_c <- pbar * (fHA - fLA)     # 该亚型的状态贡献（功能变了）
    comp_c  <- fbar * (pHA - pLA)     # 该亚型的组成贡献（比例变了）
    total_c <- pHA * fHA - pLA * fLA  # 该亚型对总差异的贡献（Σ_c total_c = dF）
    det[[paste(rg, k)]] <- data.frame(
      脑区 = rg, 通路 = pw[k], 亚型 = subtypes,
      p_HA = round(pHA, 4), p_LA = round(pLA, 4), dp = round(pHA - pLA, 4),
      f_HA = round(fHA, 4), f_LA = round(fLA, 4), df = round(fHA - fLA, 4),
      状态贡献_亚型 = round(state_c, 4),
      组成贡献_亚型 = round(comp_c, 4),
      总贡献_亚型 = round(total_c, 4),
      总贡献占比 = round(total_c / dF, 3),
      stringsAsFactors = FALSE)
  }
}

if (length(main)) {
  main <- do.call(rbind, main); det <- do.call(rbind, det)
  write_csv_bom(main, "功能差异_状态vs组成_分解.csv")
  write_csv_bom(det,  "功能差异_状态vs组成_分解_按亚型.csv")
  cat("\n=== 分解结果（状态占比≈1 状态为主 / ≈0 组成为主） ===\n")
  print(main, row.names = FALSE)
  cat("\n已写出: 功能差异_状态vs组成_分解.csv / _按亚型.csv\n")

  # ---- 摘要：每个通路 贡献最大的亚型（按总贡献占比绝对值） ----
  cat("\n=== 每个通路 主导亚型（按 |总贡献占比| 最大者） ===\n")
  for (rg in unique(det$脑区)) {
    for (pw in unique(det$通路[det$脑区 == rg])) {
      s <- det[det$脑区 == rg & det$通路 == pw, ]
      i <- which.max(abs(s$总贡献占比))
      cat(sprintf("%-6s %-45s 主导亚型 %s：状态贡献=%.4f 组成贡献=%.4f 总贡献占比=%.3f\n",
          rg, pw, s$亚型[i], s$状态贡献_亚型[i], s$组成贡献_亚型[i],
          s$总贡献占比[i]))
    }
  }
}
