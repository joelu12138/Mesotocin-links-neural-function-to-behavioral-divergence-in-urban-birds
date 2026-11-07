d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
n <- nrow(d)
sem <- function(x) sd(x)/sqrt(length(x))

cat("## 描述统计 (Mean±SEM, 频次) ##\n")
for (g in c("Urban", "Rural")) {
  s <- d[d$Group == g, ]
  cat(sprintf("%s: One = %.2f±%.2f (n=%d);  Six = %.2f±%.2f (n=%d)\n",
      g, mean(s$One), sem(s$One), nrow(s), mean(s$Six), sem(s$Six), nrow(s)))
}

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# 长格式 + 混合方差分析 (base R aov, 组内 2 水平 -> 无需球形校正)
long <- data.frame(
  y = c(d$One, d$Six),
  Group = rep(d$Group, 2),
  Cond = factor(rep(c("One", "Six"), each = n)),
  ID = factor(rep(seq_len(n), 2))
)
fit <- aov(y ~ Group * Cond + Error(ID / Cond), data = long)
s <- summary(fit)
bet <- s[["Error: ID"]][[1]]       # Group (组间)
win <- s[["Error: ID:Cond"]][[1]]  # Cond, Group:Cond (组内)

anova_tab <- data.frame(
  效应 = c("Group (组间)", "Condition (组内)", "Group × Condition (组间差异)"),
  F = round(c(bet["Group", "F value"], win["Cond", "F value"], win["Group:Cond", "F value"]), 3),
  分子df = c(bet["Group", "Df"], win["Cond", "Df"], win["Group:Cond", "Df"]),
  分母df = c(bet["Residuals", "Df"], win["Residuals", "Df"], win["Residuals", "Df"]),
  p = format.pval(c(bet["Group", "Pr(>F)"], win["Cond", "Pr(>F)"], win["Group:Cond", "Pr(>F)"]),
                  digits = 3, eps = 0.001),
  显著性 = sig(c(bet["Group", "Pr(>F)"], win["Cond", "Pr(>F)"], win["Group:Cond", "Pr(>F)"])),
  check.names = FALSE)
cat("\n## 重复度量方差分析 (组间 Group × 组内 Condition) ##\n")
print(anova_tab, row.names = FALSE)

# 分组配对 t 检验
ptab <- do.call(rbind, lapply(c("Urban", "Rural"), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$Six, paired = TRUE)
  data.frame(组 = g,
             差值One减Six = sprintf("%.2f±%.2f", mean(s2$One - s2$Six), sem(s2$One - s2$Six)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             显著性 = sig(t$p.value), check.names = FALSE)
}))
cat("\n## 分组配对 t 检验 (One vs Six, 不校正) ##\n")
print(ptab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(anova_tab, "data/聚群实验_重复度量方差分析_三线表.csv")
wr(ptab,      "data/聚群实验_配对t检验_三线表.csv")
