# Following synapse_analysis.R: two-way ANOVA GroupxSex + Fisher LSD for HP / NCL separately
regions <- c("HP", "NCL")
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}

for (reg in regions) {
  d <- read.csv(sprintf("data/%s.csv", reg), stringsAsFactors = FALSE)
  d$Group <- factor(d$Group); d$Sex <- factor(d$Sex)

  cat("\n==========", reg, "==========\n")
  cat("## Descriptive statistics (Mean+/-SEM) ##\n")
  for (s in c("Female", "Male")) for (g in c("Urban", "Rural")) {
    v <- d[d$Sex == s & d$Group == g, "MT"]
    cat(sprintf("%s %s: %.3f +/- %.3f (n=%d)\n", g, s, mean(v), sem(v), length(v)))
  }

  fit <- aov(MT ~ Group * Sex, data = d)
  a <- summary(fit)[[1]]
  rownames(a) <- trimws(rownames(a))
  cat("\n## Two-way ANOVA ##\n"); print(a)

  G <- c("Group", "Sex", "Group:Sex")
  anova_tab <- data.frame(
    Effect = c("Group (Urban vs Rural)", "Sex (biological sex)", "Group x Sex (interaction)"),
    F = round(a[G, "F value"], 3),
    df1 = a[G, "Df"], df2 = a["Residuals", "Df"],
    p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
    Significance = sig(a[G, "Pr(>F)"]), check.names = FALSE)
  cat("\n## Two-way ANOVA three-line table ##\n"); print(anova_tab, row.names = FALSE)

  # Fisher LSD: Urban vs Rural comparison within each sex (pooled error MSE and residual df from the ANOVA)
  MSE <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
  lsd_line <- function(s) {
    su <- d[d$Sex == s & d$Group == "Urban", "MT"]
    sr <- d[d$Sex == s & d$Group == "Rural", "MT"]
    se <- sqrt(MSE * (1/length(su) + 1/length(sr)))
    tval <- (mean(su) - mean(sr)) / se
    list(diff = mean(su) - mean(sr), se = se, t = tval, df = dfe,
         p = 2 * pt(-abs(tval), dfe), lsd = qt(0.975, dfe) * se)
  }
  L <- lapply(c("Female", "Male"), lsd_line)
  se_tab <- data.frame(
    Sex = c("Female", "Male"),
    Urban = sapply(c("Female", "Male"), function(s) {
      v <- d[d$Sex == s & d$Group == "Urban", "MT"]; sprintf("%.3f +/- %.3f", mean(v), sem(v)) }),
    Rural = sapply(c("Female", "Male"), function(s) {
      v <- d[d$Sex == s & d$Group == "Rural", "MT"]; sprintf("%.3f +/- %.3f", mean(v), sem(v)) }),
    Difference = round(sapply(L, `[[`, "diff"), 3),
    LSD = round(sapply(L, `[[`, "lsd"), 3),
    t = round(sapply(L, `[[`, "t"), 3),
    df = sapply(L, `[[`, "df"),
    p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
    Significance = sig(sapply(L, `[[`, "p")), check.names = FALSE)
  cat("\n## Three-line table of Urban vs Rural comparisons within each sex (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

  wr(anova_tab, sprintf("data/%s_MT_two_way_anova_three_line.csv", reg))
  wr(se_tab,    sprintf("data/%s_MT_urban_rural_comparison_three_line.csv", reg))
}
