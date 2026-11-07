# 参考 synapse_analysis.R：对 HP / NCL 分别做 Group×Sex 双因素方差分析 + Fisher LSD
regions <- c("HP", "NCL")
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}

for (reg in regions) {
  d <- read.csv(sprintf("data/%s.csv", reg), stringsAsFactors = FALSE)
  d$Group <- factor(d$Group); d$Sex <- factor(d$Sex)

  cat("\n==========", reg, "==========\n")
  cat("## 描述统计 (Mean±SEM) ##\n")
  for (s in c("Female", "Male")) for (g in c("Urban", "Rural")) {
    v <- d[d$Sex == s & d$Group == g, "MT"]
    cat(sprintf("%s %s: %.3f ± %.3f (n=%d)\n", g, s, mean(v), sem(v), length(v)))
  }

  fit <- aov(MT ~ Group * Sex, data = d)
  a <- summary(fit)[[1]]
  rownames(a) <- trimws(rownames(a))
  cat("\n## 双因素方差分析 ##\n"); print(a)

  G <- c("Group", "Sex", "Group:Sex")
  anova_tab <- data.frame(
    效应 = c("Group (城乡)", "Sex (性别)", "Group × Sex (交互)"),
    F = round(a[G, "F value"], 3),
    分子df = a[G, "Df"], 分母df = a["Residuals", "Df"],
    p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
    显著性 = sig(a[G, "Pr(>F)"]), check.names = FALSE)
  cat("\n## 双因素方差分析三线表 ##\n"); print(anova_tab, row.names = FALSE)

  # Fisher LSD: 同一性别内城乡比较（用 ANOVA 的合并误差 MSE 与残差 df）
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
    性别 = c("Female", "Male"),
    Urban = sapply(c("Female", "Male"), function(s) {
      v <- d[d$Sex == s & d$Group == "Urban", "MT"]; sprintf("%.3f ± %.3f", mean(v), sem(v)) }),
    Rural = sapply(c("Female", "Male"), function(s) {
      v <- d[d$Sex == s & d$Group == "Rural", "MT"]; sprintf("%.3f ± %.3f", mean(v), sem(v)) }),
    差值 = round(sapply(L, `[[`, "diff"), 3),
    LSD = round(sapply(L, `[[`, "lsd"), 3),
    t = round(sapply(L, `[[`, "t"), 3),
    df = sapply(L, `[[`, "df"),
    p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
    显著性 = sig(sapply(L, `[[`, "p")), check.names = FALSE)
  cat("\n## 同一性别城乡比较三线表 (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

  wr(anova_tab, sprintf("data/%s_MT_双因素方差分析_三线表.csv", reg))
  wr(se_tab,    sprintf("data/%s_MT_城乡比较_三线表.csv", reg))
}
