options(contrasts = c("contr.sum", "contr.poly"))
options(digits = 8)

d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
colnames(d)[1] <- "Group"          # first column has no header, it is actually Group
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))

cat("One + Six range:", range(d$One + d$Six), "(total duration is not fixed)\n\n")

cat("## Descriptive statistics (Mean+/-SEM, s) ##\n")
for (g in as.character(unique(d$Group))) {
  s <- d[d$Group == g, ]
  cat(sprintf("%s: One = %.3f +/- %.3f (n=%d);  Six = %.3f +/- %.3f (n=%d)\n",
      g, mean(s$One), sem(s$One), nrow(s), mean(s$Six), sem(s$Six), nrow(s)))
}

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

## ---- 1. Repeated-measures ANOVA ----
mlm <- lm(cbind(One, Six) ~ Group, data = d)
idata <- data.frame(Condition = factor(c("One", "Six")))
aov.out <- car::Anova(mlm, idata = idata, idesign = ~Condition, type = "III")
s <- summary(aov.out, multivariate = FALSE)
ut <- s$univariate.tests

rows <- c("Group", "Condition", "Group:Condition")
anova_tab <- data.frame(
  Effect = c("Group (between-group)", "Condition (within-group)", "Group x Condition (between-group difference)"),
  F = round(ut[rows, "F value"], 3),
  df1 = ut[rows, "num Df"], df2 = ut[rows, "den Df"],
  p = format.pval(ut[rows, "Pr(>F)"], digits = 3, eps = 0.001),
  Significance = sig(ut[rows, "Pr(>F)"]), check.names = FALSE)
cat("\n## Repeated-measures ANOVA (between-group Group x within-group Condition) ##\n")
print(anova_tab, row.names = FALSE)

## ---- 2. Paired t-tests by group ----
ptab <- do.call(rbind, lapply(as.character(unique(d$Group)), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$Six, paired = TRUE)
  data.frame(Group = g,
             One = sprintf("%.3f +/- %.3f", mean(s2$One), sem(s2$One)),
             Six = sprintf("%.3f +/- %.3f", mean(s2$Six), sem(s2$Six)),
             Difference_One_minus_Six = sprintf("%.3f +/- %.3f", mean(s2$One - s2$Six), sem(s2$One - s2$Six)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             Significance = sig(t$p.value), check.names = FALSE)
}))
cat("\n## Paired t-test by group (One vs Six, unadjusted) ##\n")
print(ptab, row.names = FALSE)

## ---- Write CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(anova_tab, "data/flock_size_rm_anova_three_line_table.csv")
wr(ptab,      "data/flock_size_paired_t_test_three_line_table.csv")
