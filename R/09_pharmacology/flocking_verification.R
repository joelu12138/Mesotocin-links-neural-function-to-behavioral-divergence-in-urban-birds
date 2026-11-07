options(contrasts = c("contr.sum", "contr.poly"))
options(digits = 8)

d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
colnames(d)[1] <- "Group"          # 首列无表头, 实为 Group
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))

cat("One + Six 范围:", range(d$One + d$Six), "(总时长不固定)\n\n")

cat("## 描述统计 (Mean±SEM, 秒) ##\n")
for (g in as.character(unique(d$Group))) {
  s <- d[d$Group == g, ]
  cat(sprintf("%s: One = %.3f ± %.3f (n=%d);  Six = %.3f ± %.3f (n=%d)\n",
      g, mean(s$One), sem(s$One), nrow(s), mean(s$Six), sem(s$Six), nrow(s)))
}

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

## ---- 1. 重复度量方差分析 ----
mlm <- lm(cbind(One, Six) ~ Group, data = d)
idata <- data.frame(Condition = factor(c("One", "Six")))
aov.out <- car::Anova(mlm, idata = idata, idesign = ~Condition, type = "III")
s <- summary(aov.out, multivariate = FALSE)
ut <- s$univariate.tests

rows <- c("Group", "Condition", "Group:Condition")
anova_tab <- data.frame(
  效应 = c("Group (组间)", "Condition (组内)", "Group × Condition (组间差异)"),
  F = round(ut[rows, "F value"], 3),
  分子df = ut[rows, "num Df"], 分母df = ut[rows, "den Df"],
  p = format.pval(ut[rows, "Pr(>F)"], digits = 3, eps = 0.001),
  显著性 = sig(ut[rows, "Pr(>F)"]), check.names = FALSE)
cat("\n## 重复度量方差分析 (组间 Group × 组内 Condition) ##\n")
print(anova_tab, row.names = FALSE)

## ---- 2. 分组配对 t 检验 ----
ptab <- do.call(rbind, lapply(as.character(unique(d$Group)), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$Six, paired = TRUE)
  data.frame(组 = g,
             One = sprintf("%.3f ± %.3f", mean(s2$One), sem(s2$One)),
             Six = sprintf("%.3f ± %.3f", mean(s2$Six), sem(s2$Six)),
             差值One减Six = sprintf("%.3f ± %.3f", mean(s2$One - s2$Six), sem(s2$One - s2$Six)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             显著性 = sig(t$p.value), check.names = FALSE)
}))
cat("\n## 分组配对 t 检验 (One vs Six, 不校正) ##\n")
print(ptab, row.names = FALSE)

## ---- 写 CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(anova_tab, "data/聚群实验_重复度量方差分析_三线表.csv")
wr(ptab,      "data/聚群实验_配对t检验_三线表.csv")
