raw <- read.csv("data/AA.csv", check.names = FALSE)
metric_name <- names(raw)[3]
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
MSE <- a["Residuals", "Mean Sq"]; df_resid <- a["Residuals", "Df"]
cat("\n## 双因素方差分析 ##\n"); print(a)

G <- c("Group", "Region", "Group:Region")
anova_tab <- data.frame(
  效应 = c("Group (城乡)", "Region (脑区)", "Group × Region (交互)"),
  F = round(a[G, "F value"], 3),
  分子df = a[G, "Df"], 分母df = df_resid,
  p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
  显著性 = sig(a[G, "Pr(>F)"]), check.names = FALSE)
cat("\n## 双因素方差分析三线表 ##\n"); print(anova_tab, row.names = FALSE)

# 简单效应 (城乡比较, 用合并误差 MSE, 不校正)
se_tab <- do.call(rbind, lapply(c("HP", "NCL"), function(r) {
  u <- d[d$Region == r & d$Group == "Urban", "y"]
  rv <- d[d$Region == r & d$Group == "Rural", "y"]
  diff <- mean(u) - mean(rv)
  se <- sqrt(MSE * (1/length(u) + 1/length(rv)))
  tt <- diff / se
  p <- 2 * pt(-abs(tt), df_resid)
  data.frame(脑区 = r,
             Urban = sprintf("%.3f±%.3f", mean(u), sem(u)),
             Rural = sprintf("%.3f±%.3f", mean(rv), sem(rv)),
             t = round(tt, 3), df = df_resid,
             p = format.pval(p, digits = 3, eps = 0.001),
             显著性 = sig(p), check.names = FALSE)
}))
cat("\n## 同一脑区城乡比较三线表 (不校正) ##\n"); print(se_tab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(anova_tab, paste0("data/", metric_name, "_双因素方差分析_三线表.csv"))
wr(se_tab,    paste0("data/", metric_name, "_脑区城乡比较_三线表.csv"))
