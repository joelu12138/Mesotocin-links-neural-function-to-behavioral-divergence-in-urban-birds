options(digits = 6)

d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

d$PI <- (d$Six - d$One) / (d$Six + d$One)   # 偏好指数, >0 = 偏好 6 鸟笼

## ---- 1. 组内成对比较：One vs Six（配对 Wilcoxon 符号秩） ----
cat("## 1. 组内成对比较 (One vs Six, 配对 Wilcoxon 符号秩) ##\n")
pair <- do.call(rbind, lapply(c("总体", levels(d$Group)), function(g) {
  if (g == "总体") { x <- d$One; y <- d$Six } else { x <- d$One[d$Group == g]; y <- d$Six[d$Group == g] }
  wt <- wilcox.test(x, y, paired = TRUE)
  data.frame(组 = g, n = length(x),
             One中位 = median(x), Six中位 = median(y),
             V = unname(wt$statistic),
             p = format.pval(wt$p.value, digits = 3, eps = 0.001),
             显著性 = sig(wt$p.value), check.names = FALSE)
}))
print(pair, row.names = FALSE)

## ---- 2. 组间比较：Mann-Whitney U ----
cat("\n## 2. 组间比较 (Urban vs Rural, Mann-Whitney U) ##\n")
comp <- do.call(rbind, lapply(
  list(c("One 频次", "One", "%.1f"),
       c("Six 频次", "Six", "%.1f"),
       c("总频次",   "Tot", "%.1f"),
       c("偏好指数 PI", "PI", "%.3f")),
  function(spec) {
    lab <- spec[1]; v <- spec[2]; fmt <- spec[3]
    x <- if (v == "Tot") d$One + d$Six else d[[v]]
    u <- x[d$Group == "Urban"]; r <- x[d$Group == "Rural"]
    wt <- wilcox.test(u, r)
    data.frame(指标 = lab,
               Urban = sprintf(paste0(fmt, "±", fmt), mean(u), sem(u)),
               Rural = sprintf(paste0(fmt, "±", fmt), mean(r), sem(r)),
               W = unname(wt$statistic),
               p = format.pval(wt$p.value, digits = 3, eps = 0.001),
               显著性 = sig(wt$p.value), check.names = FALSE)
  }
))
print(comp, row.names = FALSE)

## ---- 写 CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(pair, "data/聚群实验_组内配对Wilcoxon.csv")
wr(comp, "data/聚群实验_组间MannWhitneyU.csv")
