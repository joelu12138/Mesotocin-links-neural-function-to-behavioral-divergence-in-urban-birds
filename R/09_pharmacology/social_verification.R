options(contrasts = c("contr.sum", "contr.poly"))
options(digits = 8)

d <- read.csv("data/Total exploration time.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))

# Check whether total duration is fixed
cat("One + No range:", range(d$One + d$No), "(fixed at", d$One[1] + d$No[1], "s per bird)\n\n")

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

## ---- 1. Repeated-measures ANOVA ----
mlm <- lm(cbind(One, No) ~ Group, data = d)
idata <- data.frame(Condition = factor(c("One", "No")))
aov.out <- car::Anova(mlm, idata = idata, idesign = ~Condition, type = "III")
s <- summary(aov.out, multivariate = FALSE)
ut <- s$univariate.tests

anova_tab <- data.frame(
  Effect = c("Condition (within-group)", "Group x Condition (between-group difference)"),
  F = round(c(ut["Condition", "F value"], ut["Group:Condition", "F value"]), 3),
  df1 = c(ut["Condition", "num Df"], ut["Group:Condition", "num Df"]),
  df2 = c(ut["Condition", "den Df"], ut["Group:Condition", "den Df"]),
  p = c(format.pval(ut["Condition", "Pr(>F)"], digits = 3, eps = 0.001),
        format.pval(ut["Group:Condition", "Pr(>F)"], digits = 3, eps = 0.001)),
  Significance = c(sig(ut["Condition", "Pr(>F)"]), sig(ut["Group:Condition", "Pr(>F)"])),
  check.names = FALSE)
cat("## Repeated-measures ANOVA (between-group Group x within-group Condition) ##\n")
print(anova_tab, row.names = FALSE)
cat("Note: the Group main effect is degenerate (One+No fixed at 300s, no between-group variance in total duration); the between-group difference is captured by the interaction.\n\n")

## ---- 2. Paired t-tests by group ----
ptab <- do.call(rbind, lapply(as.character(unique(d$Group)), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$No, paired = TRUE)
  data.frame(Group = g,
             One = sprintf("%.3f +/- %.3f", mean(s2$One), sem(s2$One)),
             No = sprintf("%.3f +/- %.3f", mean(s2$No), sem(s2$No)),
             Difference_One_minus_No = sprintf("%.3f +/- %.3f", mean(s2$One - s2$No), sem(s2$One - s2$No)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             Significance = sig(t$p.value), check.names = FALSE)
}))
cat("## Paired t-test by group (One vs No, unadjusted) ##\n")
print(ptab, row.names = FALSE)

# Between-group supplement: independent t-test on One (equivalent to the interaction)
w <- t.test(d$One[d$Group == as.character(unique(d$Group))[1]], d$One[d$Group == as.character(unique(d$Group))[2]])
cat(sprintf("\nBetween-group independent t-test on One: t=%.3f, df=%.2f, p=%s (equivalent to the interaction)\n",
            w$statistic, w$parameter, format.pval(w$p.value, digits = 3)))

## ---- Write CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(anova_tab, "data/social_rm_anova_three_line_table.csv")
wr(ptab,      "data/social_paired_t_test_three_line_table.csv")
