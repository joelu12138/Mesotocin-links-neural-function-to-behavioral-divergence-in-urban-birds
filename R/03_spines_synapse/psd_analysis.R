d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
d$Group <- factor(d$Group); d$Region <- factor(d$Region)
y <- setdiff(names(d), c("Group", "Region"))[1]   # automatically take the metric column
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cat("Metric:", y, "\n\n## Descriptive statistics (Mean+/-SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  s <- d[d$Region == r & d$Group == g, y]
  cat(sprintf("%s %s: %.3f+/-%.3f (n=%d)\n", g, r, mean(s), sem(s), length(s)))
}

fit <- aov(as.formula(paste(y, "~ Group * Region")), data = d)
a <- summary(fit)[[1]]
rownames(a) <- trimws(rownames(a))
cat("\n## Two-way ANOVA ##\n"); print(a)

G <- c("Group", "Region", "Group:Region")
anova_tab <- data.frame(
  Effect = c("Group (Urban vs Rural)", "Region (brain region)", "Group x Region (interaction)"),
  F = round(a[G, "F value"], 3),
  df1 = a[G, "Df"], df2 = a["Residuals", "Df"],
  p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
  Significance = sig(a[G, "Pr(>F)"]), check.names = FALSE)
cat("\n## Two-way ANOVA three-line table ##\n"); print(anova_tab, row.names = FALSE)

# Fisher LSD: uncorrected t tests using the pooled error (MSE) and residual df of the two-way ANOVA
MSE <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
lsd_line <- function(r) {
  su <- d[d$Region == r & d$Group == "Urban", y]
  sr <- d[d$Region == r & d$Group == "Rural", y]
  se <- sqrt(MSE * (1/length(su) + 1/length(sr)))
  tval <- (mean(su) - mean(sr)) / se
  list(diff = mean(su) - mean(sr), se = se, t = tval, df = dfe,
       p = 2 * pt(-abs(tval), dfe), lsd = qt(0.975, dfe) * se)
}
L <- lapply(c("HP", "NCL"), lsd_line)
se_tab <- data.frame(
  brain_region = c("HP", "NCL"),
  Urban = sapply(c("HP", "NCL"), function(r) {
    s <- d[d$Region == r & d$Group == "Urban", y]; sprintf("%.3f+/-%.3f", mean(s), sem(s)) }),
  Rural = sapply(c("HP", "NCL"), function(r) {
    s <- d[d$Region == r & d$Group == "Rural", y]; sprintf("%.3f+/-%.3f", mean(s), sem(s)) }),
  Difference = round(sapply(L, `[[`, "diff"), 3),
  LSD = round(sapply(L, `[[`, "lsd"), 3),
  t = round(sapply(L, `[[`, "t"), 3),
  df = sapply(L, `[[`, "df"),
  p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
  Significance = sig(sapply(L, `[[`, "p")), check.names = FALSE)
cat("\n## Three-line table of Urban vs Rural comparisons within each brain region (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(anova_tab, paste0("data/", y, "_two_way_anova_three_line.csv"))
wr(se_tab,    paste0("data/", y, "_brain_region_urban_rural_comparison_three_line.csv"))
