raw <- read.csv("data/AA.csv", check.names = FALSE)
metric_name <- names(raw)[3]                       # 原始指标名 (可能含 - 等特殊字符)
d <- data.frame(Group = factor(raw[[1]]), Region = factor(raw[[2]]), y = as.numeric(raw[[3]]))
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cat("指标:", metric_name, "\n\n## 描述统计 (Mean±SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  s <- d[d$Region == r & d$Group == g, "y"]
  cat(sprintf("%s %s: %.3f±%.3f (n=%d)\n", g, r, mean(s), sem(s), length(s)))
}

fit <- aov(y ~ Group * Region, data = d)
a <- summary(fit)[[1]]
rownames(a) <- trimws(rownames(a))
cat("\n## 双因素方差分析 ##\n"); print(a)

G <- c("Group", "Region", "Group:Region")
anova_tab <- data.frame(
  效应 = c("Group (城乡)", "Region (脑区)", "Group × Region (交互)"),
  F = round(a[G, "F value"], 3),
  分子df = a[G, "Df"], 分母df = a["Residuals", "Df"],
  p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
  显著性 = sig(a[G, "Pr(>F)"]), check.names = FALSE)
cat("\n## 双因素方差分析三线表 ##\n"); print(anova_tab, row.names = FALSE)

# Fisher LSD: 用双因素ANOVA的合并误差(MSE)与残差df做不校正的 t 检验
MSE <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
lsd_line <- function(r) {
  su <- d[d$Region == r & d$Group == "Urban", "y"]
  sr <- d[d$Region == r & d$Group == "Rural", "y"]
  se <- sqrt(MSE * (1/length(su) + 1/length(sr)))
  tval <- (mean(su) - mean(sr)) / se
  list(diff = mean(su) - mean(sr), se = se, t = tval, df = dfe,
       p = 2 * pt(-abs(tval), dfe), lsd = qt(0.975, dfe) * se)
}
L <- lapply(c("HP", "NCL"), lsd_line)
se_tab <- data.frame(
  脑区 = c("HP", "NCL"),
  Urban = sapply(c("HP", "NCL"), function(r) {
    s <- d[d$Region == r & d$Group == "Urban", "y"]; sprintf("%.3f±%.3f", mean(s), sem(s)) }),
  Rural = sapply(c("HP", "NCL"), function(r) {
    s <- d[d$Region == r & d$Group == "Rural", "y"]; sprintf("%.3f±%.3f", mean(s), sem(s)) }),
  差值 = round(sapply(L, `[[`, "diff"), 3),
  LSD = round(sapply(L, `[[`, "lsd"), 3),
  t = round(sapply(L, `[[`, "t"), 3),
  df = sapply(L, `[[`, "df"),
  p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
  显著性 = sig(sapply(L, `[[`, "p")), check.names = FALSE)
cat("\n## 同一脑区城乡比较三线表 (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(anova_tab, paste0("data/", metric_name, "_双因素方差分析_三线表.csv"))
wr(se_tab,    paste0("data/", metric_name, "_脑区城乡比较_三线表.csv"))
