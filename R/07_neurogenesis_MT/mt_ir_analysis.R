sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}

## ================= AA_HPNCL.csv: two-way ANOVA + Fisher LSD =================
d1 <- read.csv("data/AA_HPNCL.csv", stringsAsFactors = FALSE)
y <- setdiff(names(d1), c("Group", "Region"))[1]
names(d1)[names(d1) == y] <- "MT"
d1$Group <- factor(d1$Group); d1$Region <- factor(d1$Region)

cat("========== AA_HPNCL (MT-ir_cell) ==========\n")
cat("## Descriptive statistics (Mean+/-SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  v <- d1[d1$Region == r & d1$Group == g, "MT"]
  cat(sprintf("%s %s: %.3f +/- %.3f (n=%d)\n", g, r, mean(v), sem(v), length(v)))
}

fit <- aov(MT ~ Group * Region, data = d1)
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

# Fisher LSD: Urban vs Rural comparison within the same brain region
MSE <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
lsd_line <- function(r) {
  su <- d1[d1$Region == r & d1$Group == "Urban", "MT"]
  sr <- d1[d1$Region == r & d1$Group == "Rural", "MT"]
  se <- sqrt(MSE * (1/length(su) + 1/length(sr)))
  tval <- (mean(su) - mean(sr)) / se
  list(diff = mean(su) - mean(sr), se = se, t = tval, df = dfe,
       p = 2 * pt(-abs(tval), dfe), lsd = qt(0.975, dfe) * se)
}
L <- lapply(c("HP", "NCL"), lsd_line)
se_tab <- data.frame(
  brain_region = c("HP", "NCL"),
  Urban = sapply(c("HP", "NCL"), function(r) {
    v <- d1[d1$Region == r & d1$Group == "Urban", "MT"]; sprintf("%.3f +/- %.3f", mean(v), sem(v)) }),
  Rural = sapply(c("HP", "NCL"), function(r) {
    v <- d1[d1$Region == r & d1$Group == "Rural", "MT"]; sprintf("%.3f +/- %.3f", mean(v), sem(v)) }),
  Difference = round(sapply(L, `[[`, "diff"), 3),
  LSD = round(sapply(L, `[[`, "lsd"), 3),
  t = round(sapply(L, `[[`, "t"), 3),
  df = sapply(L, `[[`, "df"),
  p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
  Significance = sig(sapply(L, `[[`, "p")), check.names = FALSE)
cat("\n## Three-line table of Urban vs Rural comparisons within each brain region (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

wr(anova_tab, "data/AA_HPNCL_two_way_anova_three_line.csv")
wr(se_tab,    "data/AA_HPNCL_urban_rural_comparison_three_line.csv")

## ================= AA_PVN.csv: independent-samples t test =================
d2 <- read.csv("data/AA_PVN.csv", stringsAsFactors = FALSE)
tt <- t.test(d2$Urban, d2$Rural)   # Welch (default, does not assume equal variances)

desc <- data.frame(
  Group = c("Urban", "Rural"),
  n = c(length(d2$Urban), length(d2$Rural)),
  `Mean+/-SEM` = c(sprintf("%.3f +/- %.3f", mean(d2$Urban), sem(d2$Urban)),
                 sprintf("%.3f +/- %.3f", mean(d2$Rural), sem(d2$Rural))),
  check.names = FALSE)

tt_tab <- data.frame(
  Comparison = "Urban vs Rural",
  t = round(as.numeric(tt$statistic), 3),
  df = round(as.numeric(tt$parameter), 3),
  p = format.pval(tt$p.value, digits = 3, eps = 0.001),
  Significance = sig(tt$p.value), check.names = FALSE)

cat("\n========== AA_PVN: independent-samples t test ==========\n")
cat("## Descriptive statistics (Mean+/-SEM) ##\n"); print(desc, row.names = FALSE)
cat("\n## t test three-line table ##\n"); print(tt_tab, row.names = FALSE)

wr(desc,   "data/AA_PVN_descriptive_statistics_three_line.csv")
wr(tt_tab, "data/AA_PVN_t_test_three_line.csv")
