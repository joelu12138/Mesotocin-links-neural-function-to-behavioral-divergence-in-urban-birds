# 单因素重复度量方差分析 (Condition 组内) + Bonferroni 多重比较
# Urban / Rural 完全分开分析, 各生成 2 个表 (方差分析 + 事后检验), 共 4 个 CSV
# Urban 处理列 = OTA; Rural 处理列 = MT; Bonferroni 按每组 3 个比较单独校正

library(dplyr)
library(tidyr)

sig <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("﻿", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}

groups <- list(
  Urban = list(file = "Urban.csv", conds = c("Baseline", "OTA", "Washout")),
  Rural = list(file = "Rural.csv", conds = c("Baseline", "MT",  "Washout"))
)

for (g in names(groups)) {
  dat   <- read.csv(groups[[g]]$file)
  conds <- groups[[g]]$conds
  k     <- length(conds)

  long <- dat %>% pivot_longer(all_of(conds), names_to = "Condition", values_to = "Value")
  long$Subject   <- factor(long$Indvi)
  long$Condition <- factor(long$Condition, levels = conds)

  # ---- 单因素重复度量方差分析 (主效应) + GG 校正 ----
  a <- aov(Value ~ Condition + Error(Subject), data = long)
  within <- summary(a)[["Error: Within"]][[1]]
  Fval <- within["Condition", "F value"]; df1 <- within["Condition", "Df"]; df2 <- within["Residuals", "Df"]
  p    <- within["Condition", "Pr(>F)"]

  wide <- long %>% select(Indvi, Condition, Value) %>% pivot_wider(names_from = Condition, values_from = Value)
  Y <- as.matrix(wide[, conds])
  C <- contr.helmert(k); Q <- qr.Q(qr(C)); S <- cov(Y %*% Q)
  GG <- sum(diag(S))^2 / ((k - 1) * sum(S^2))
  p_gg <- pf(Fval, (k - 1) * GG, df2 * GG, lower.tail = FALSE)

  anova_tab <- data.frame(效应 = "Condition (组内)",
                          F = Fval, 分子df = df1, 分母df = df2, p = p,
                          GGε = GG, GG校正p = p_gg)
  anova_tab$显著性 <- sig(anova_tab$GG校正p)
  anova_tab$F     <- round(anova_tab$F, 3)
  anova_tab$GGε   <- round(anova_tab$GGε, 4)
  anova_tab$p     <- format.pval(anova_tab$p,     digits = 3, eps = 0.001)
  anova_tab$GG校正p <- format.pval(anova_tab$GG校正p, digits = 3, eps = 0.001)

  # ---- 事后多重比较: 配对 t + Bonferroni (每组 3 个比较) ----
  post <- do.call(rbind, lapply(1:(k - 1), function(i) {
    do.call(rbind, lapply((i + 1):k, function(j) {
      tt <- t.test(long$Value[long$Condition == conds[i]], long$Value[long$Condition == conds[j]], paired = TRUE)
      data.frame(比较 = paste0(conds[i], " vs ", conds[j]),
                 t = tt$statistic, df = tt$parameter, p = tt$p.value)
    }))
  }))
  post$Bonferroni <- p.adjust(post$p, "bonferroni")
  post$显著性 <- sig(post$Bonferroni)
  post$t   <- round(post$t, 3); post$df <- round(post$df, 2)
  post$p   <- format.pval(post$p,   digits = 3, eps = 0.001)
  post$Bonferroni <- format.pval(post$Bonferroni, digits = 3, eps = 0.001)

  cat(sprintf("\n========== %s ==========\n", g))
  print(anova_tab, row.names = FALSE); print(post, row.names = FALSE)

  wr(anova_tab, paste0(g, "_rm_anova_results.csv"))
  wr(post,      paste0(g, "_posthoc_results.csv"))
}
