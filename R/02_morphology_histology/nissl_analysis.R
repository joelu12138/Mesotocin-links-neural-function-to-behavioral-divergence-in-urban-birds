# 尼氏小体 双因素方差分析 + Fisher's LSD 多重比较 (复现 GraphPad Prism)
d <- read.csv("data/Nissl bodies.csv", stringsAsFactors = FALSE)
d$Group  <- factor(d$Group,  levels = c("Urban", "Rural"))
d$Region <- factor(d$Region, levels = c("HP", "NCL"))

sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                       ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
pfmt <- function(p) ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))

cat("## 描述统计 (Mean±SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  s <- d[d$Region == r & d$Group == g, "Niss"]
  cat(sprintf("%s %s: %.3f ± %.3f (n=%d)\n", g, r, mean(s), sem(s), length(s)))
}

# ---- 双因素方差分析 ----
fit <- aov(Niss ~ Group * Region, data = d)
a <- summary(fit)[[1]]
rownames(a) <- trimws(rownames(a))
cat("\n## 双因素方差分析 ##\n"); print(a)

# Fisher's LSD 用全模型合并误差项: 残差均方 MSE 与残差自由度 df
MSE    <- a["Residuals", "Mean Sq"]
df_res <- a["Residuals", "Df"]

# ---- Fisher's LSD (未校正): 每行(脑区)内 Urban - Rural ----
alpha <- 0.05
tcrit <- qt(1 - alpha/2, df_res)   # 未校正临界 t, 用于 95% CI

res <- do.call(rbind, lapply(c("HP", "NCL"), function(r) {
  mu <- mean(d$Niss[d$Region == r & d$Group == "Urban"])
  mr <- mean(d$Niss[d$Region == r & d$Group == "Rural"])
  nu <- sum(d$Region == r & d$Group == "Urban")
  nr <- sum(d$Region == r & d$Group == "Rural")
  diff <- mu - mr
  se   <- sqrt(MSE * (1/nu + 1/nr))
  tval <- diff / se
  p    <- 2 * pt(-abs(tval), df_res)   # 未校正双侧 p
  data.frame(Region = r, Mean1 = mu, Mean2 = mr, Diff = diff, SE = se,
             N1 = nu, N2 = nr, t = tval, DF = df_res,
             CI.lo = diff - tcrit*se, CI.hi = diff + tcrit*se, p = p)
}))

cat("\nCompare each cell mean with the other cell mean in that row\n")
cat(sprintf("Number of families: 1\nNumber of comparisons per family: 2\nAlpha: %.2f\n", alpha))

tab1 <- data.frame(
  Region = res$Region,
  `Mean Diff.` = round(res$Diff, 2),
  `95.00% CI of diff.` = sprintf("%.3f to %.3f", res$CI.lo, res$CI.hi),
  `Significant?` = ifelse(res$p < alpha, "Yes", "No"),
  Summary = sig(res$p),
  `Individual P Value` = pfmt(res$p),
  check.names = FALSE)
cat("\nUncorrected Fisher's LSD\n  Urban - Rural\n")
print(tab1, row.names = FALSE)

tab2 <- data.frame(
  Region = res$Region,
  `Mean 1` = round(res$Mean1, 2),
  `Mean 2` = round(res$Mean2, 2),
  `Mean Diff.` = round(res$Diff, 2),
  `SE of diff.` = round(res$SE, 3),
  N1 = res$N1, N2 = res$N2,
  t = round(res$t, 3),
  DF = sprintf("%.2f", res$DF),
  check.names = FALSE)
cat("\nTest details\n  Urban - Rural\n")
print(tab2, row.names = FALSE)

# ---- 写出 CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(tab1, "data/尼氏小体_FisherLSD多重比较.csv")
wr(tab2, "data/尼氏小体_FisherLSD检验细节.csv")
