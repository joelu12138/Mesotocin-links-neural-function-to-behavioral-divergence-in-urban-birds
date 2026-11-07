# One-way repeated measures ANOVA (Condition within-group) + Bonferroni multiple comparisons
# Urban / Rural are analysed completely separately, each producing 2 tables (ANOVA + post hoc test), 4 CSVs in total
# Urban treatment columns = OTA; Rural treatment columns = MT; Bonferroni corrected separately for 3 comparisons per group

library(dplyr)
library(tidyr)

sig <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
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

  # ---- One-way repeated measures ANOVA (main effect) + GG correction ----
  a <- aov(Value ~ Condition + Error(Subject), data = long)
  within <- summary(a)[["Error: Within"]][[1]]
  Fval <- within["Condition", "F value"]; df1 <- within["Condition", "Df"]; df2 <- within["Residuals", "Df"]
  p    <- within["Condition", "Pr(>F)"]

  wide <- long %>% select(Indvi, Condition, Value) %>% pivot_wider(names_from = Condition, values_from = Value)
  Y <- as.matrix(wide[, conds])
  C <- contr.helmert(k); Q <- qr.Q(qr(C)); S <- cov(Y %*% Q)
  GG <- sum(diag(S))^2 / ((k - 1) * sum(S^2))
  p_gg <- pf(Fval, (k - 1) * GG, df2 * GG, lower.tail = FALSE)

  anova_tab <- data.frame(Effect = "Condition (within-group)",
                          F = Fval, df1 = df1, df2 = df2, p = p,
                          GG_epsilon = GG, GG_corrected_p = p_gg)
  anova_tab$Significance <- sig(anova_tab$GG_corrected_p)
  anova_tab$F     <- round(anova_tab$F, 3)
  anova_tab$GG_epsilon   <- round(anova_tab$GG_epsilon, 4)
  anova_tab$p     <- format.pval(anova_tab$p,     digits = 3, eps = 0.001)
  anova_tab$GG_corrected_p <- format.pval(anova_tab$GG_corrected_p, digits = 3, eps = 0.001)

  # ---- Post hoc multiple comparisons: paired t + Bonferroni (3 comparisons per group) ----
  post <- do.call(rbind, lapply(1:(k - 1), function(i) {
    do.call(rbind, lapply((i + 1):k, function(j) {
      tt <- t.test(long$Value[long$Condition == conds[i]], long$Value[long$Condition == conds[j]], paired = TRUE)
      data.frame(Comparison = paste0(conds[i], " vs ", conds[j]),
                 t = tt$statistic, df = tt$parameter, p = tt$p.value)
    }))
  }))
  post$Bonferroni <- p.adjust(post$p, "bonferroni")
  post$Significance <- sig(post$Bonferroni)
  post$t   <- round(post$t, 3); post$df <- round(post$df, 2)
  post$p   <- format.pval(post$p,   digits = 3, eps = 0.001)
  post$Bonferroni <- format.pval(post$Bonferroni, digits = 3, eps = 0.001)

  cat(sprintf("\n========== %s ==========\n", g))
  print(anova_tab, row.names = FALSE); print(post, row.names = FALSE)

  wr(anova_tab, paste0(g, "_rm_anova_results.csv"))
  wr(post,      paste0(g, "_posthoc_results.csv"))
}
