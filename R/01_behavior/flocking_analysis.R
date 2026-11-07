d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
n <- nrow(d)
sem <- function(x) sd(x)/sqrt(length(x))

cat("## Descriptive statistics (Mean+/-SEM, frequency) ##\n")
for (g in c("Urban", "Rural")) {
  s <- d[d$Group == g, ]
  cat(sprintf("%s: One = %.2f+/-%.2f (n=%d);  Six = %.2f+/-%.2f (n=%d)\n",
      g, mean(s$One), sem(s$One), nrow(s), mean(s$Six), sem(s$Six), nrow(s)))
}

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# long format + mixed ANOVA (base R aov, within-group 2 levels -> no sphericity correction needed)
long <- data.frame(
  y = c(d$One, d$Six),
  Group = rep(d$Group, 2),
  Cond = factor(rep(c("One", "Six"), each = n)),
  ID = factor(rep(seq_len(n), 2))
)
fit <- aov(y ~ Group * Cond + Error(ID / Cond), data = long)
s <- summary(fit)
bet <- s[["Error: ID"]][[1]]       # Group (between-group)
win <- s[["Error: ID:Cond"]][[1]]  # Cond, Group:Cond (within-group)

anova_tab <- data.frame(
  Effect = c("Group (between-group)", "Condition (within-group)", "Group x Condition (between-group difference)"),
  F = round(c(bet["Group", "F value"], win["Cond", "F value"], win["Group:Cond", "F value"]), 3),
  df1 = c(bet["Group", "Df"], win["Cond", "Df"], win["Group:Cond", "Df"]),
  df2 = c(bet["Residuals", "Df"], win["Residuals", "Df"], win["Residuals", "Df"]),
  p = format.pval(c(bet["Group", "Pr(>F)"], win["Cond", "Pr(>F)"], win["Group:Cond", "Pr(>F)"]),
                  digits = 3, eps = 0.001),
  Significance = sig(c(bet["Group", "Pr(>F)"], win["Cond", "Pr(>F)"], win["Group:Cond", "Pr(>F)"])),
  check.names = FALSE)
cat("\n## Repeated-measures ANOVA (between-group Group x within-group Condition) ##\n")
print(anova_tab, row.names = FALSE)

# paired t-test by group
ptab <- do.call(rbind, lapply(c("Urban", "Rural"), function(g) {
  s2 <- d[d$Group == g, ]
  t <- t.test(s2$One, s2$Six, paired = TRUE)
  data.frame(Group = g,
             DifferenceOneMinusSix = sprintf("%.2f+/-%.2f", mean(s2$One - s2$Six), sem(s2$One - s2$Six)),
             t = round(t$statistic, 3), df = t$parameter,
             p = format.pval(t$p.value, digits = 3, eps = 0.001),
             Significance = sig(t$p.value), check.names = FALSE)
}))
cat("\n## Paired t-test by group (One vs Six, uncorrected) ##\n")
print(ptab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(anova_tab, "data/flock_size_rm_anova_three_line.csv")
wr(ptab,      "data/flock_size_paired_t_three_line.csv")
