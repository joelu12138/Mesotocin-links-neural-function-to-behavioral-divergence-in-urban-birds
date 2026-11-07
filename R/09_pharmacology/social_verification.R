options(contrasts = c("contr.sum", "contr.poly"))
options(digits = 8)

d <- read.csv("data/Total exploration time.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))

# 检查总时长是否固定
cat("One + No 范围:", range(d$One + d$No), "(每只鸟固定为", d$One[1] + d$No[1], "秒)\n\n")

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

## ---- 1. 重复度量方差分析 ----
mlm <- lm(cbind(One, No) ~ Group, data = d)
idata <- data.frame(Condition = factor(c("One", "No")))
aov.out <- car::Anova(mlm, idata = idata, idesign = ~Condition, type = "III")
s <- summary(aov.out, multivariate = FALSE)
ut <- s$univariate.tests

anova_tab <- data.frame(
  效应 = c("Condition (组内)", "Group × Condition (组间差异)"),
  F = round(c(ut["Condition", "F value"], ut["Group:Condition", "F value"]), 3),
  分子df = c(ut["Condition", "num Df"], ut["Group:Condition", "num Df"]),
  分母df = c(ut["Condition", "den Df"], ut["Group:Condition", "den Df"]),
  p = c(format.pval(ut["Condition", "Pr(>F)"], digits = 3, eps = 0.001),
        format.pval(ut["Group:Condition", "Pr(>F)"], digits = 3, eps = 0.001)),
  显著性 = c(sig(ut["Condition", "Pr(>F)"]), sig(ut["Group:Condition", "Pr(>F)"])),
  check.names = FALSE)
cat("## 重复度量方差分析 (组间 Group × 组内 Condition) ##\n")
print(anova_tab, row.names = FALSE)
cat("注: Group 主效应退化(One+No 固定为 300s, 组间总时长无方差), 组间差异由交互项体现。\n\n")

## ---- 2. 分组配对 t 检验 ----
ptab <- do.call(rbind, lapply(as.character(unique(d$Group)), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$No, paired = TRUE)
  data.frame(组 = g,
             One = sprintf("%.3f ± %.3f", mean(s2$One), sem(s2$One)),
             No = sprintf("%.3f ± %.3f", mean(s2$No), sem(s2$No)),
             差值One减No = sprintf("%.3f ± %.3f", mean(s2$One - s2$No), sem(s2$One - s2$No)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             显著性 = sig(t$p.value), check.names = FALSE)
}))
cat("## 分组配对 t 检验 (One vs No, 不校正) ##\n")
print(ptab, row.names = FALSE)

# 组间差异补充: One 的组间独立 t 检验 (与交互项等价)
w <- t.test(d$One[d$Group == as.character(unique(d$Group))[1]], d$One[d$Group == as.character(unique(d$Group))[2]])
cat(sprintf("\n组间 One 独立 t 检验: t=%.3f, df=%.2f, p=%s (与交互项等价)\n",
            w$statistic, w$parameter, format.pval(w$p.value, digits = 3)))

## ---- 写 CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(anova_tab, "data/社会实验_重复度量方差分析_三线表.csv")
wr(ptab,      "data/社会实验_配对t检验_三线表.csv")
