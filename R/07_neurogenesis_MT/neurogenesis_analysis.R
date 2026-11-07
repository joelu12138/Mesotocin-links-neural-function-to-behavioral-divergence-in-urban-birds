raw <- read.csv("data/AA.csv", check.names = FALSE)
metric_name <- names(raw)[3]
d <- data.frame(Group = factor(raw[[1]]), Region = factor(raw[[2]]), y = as.numeric(raw[[3]]))
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cat("Metric:", metric_name, "\n\n## Descriptive statistics (Mean+/-SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  s <- d[d$Region == r & d$Group == g, "y"]
  cat(sprintf("%s %s: %.3f+/-%.3f (n=%d)\n", g, r, mean(s), sem(s), length(s)))
}

fit <- aov(y ~ Group * Region, data = d)
a <- summary(fit)[[1]]
rownames(a) <- trimws(rownames(a))
MSE <- a["Residuals", "Mean Sq"]; df_resid <- a["Residuals", "Df"]
cat("\n## Two-way ANOVA ##\n"); print(a)

G <- c("Group", "Region", "Group:Region")
anova_tab <- data.frame(
  Effect = c("Group (Urban vs Rural)", "Region (brain region)", "Group x Region (interaction)"),
  F = round(a[G, "F value"], 3),
  df1 = a[G, "Df"], df2 = df_resid,
  p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
  Significance = sig(a[G, "Pr(>F)"]), check.names = FALSE)
cat("\n## Two-way ANOVA three-line table ##\n"); print(anova_tab, row.names = FALSE)

# Simple effects (Urban vs Rural comparison, pooled error MSE, uncorrected)
se_tab <- do.call(rbind, lapply(c("HP", "NCL"), function(r) {
  u <- d[d$Region == r & d$Group == "Urban", "y"]
  rv <- d[d$Region == r & d$Group == "Rural", "y"]
  diff <- mean(u) - mean(rv)
  se <- sqrt(MSE * (1/length(u) + 1/length(rv)))
  tt <- diff / se
  p <- 2 * pt(-abs(tt), df_resid)
  data.frame(brain_region = r,
             Urban = sprintf("%.3f+/-%.3f", mean(u), sem(u)),
             Rural = sprintf("%.3f+/-%.3f", mean(rv), sem(rv)),
             t = round(tt, 3), df = df_resid,
             p = format.pval(p, digits = 3, eps = 0.001),
             Significance = sig(p), check.names = FALSE)
}))
cat("\n## Three-line table of Urban vs Rural comparisons within each brain region (uncorrected) ##\n"); print(se_tab, row.names = FALSE)

wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(anova_tab, paste0("data/", metric_name, "_two_way_anova_three_line.csv"))
wr(se_tab,    paste0("data/", metric_name, "_brain_region_urban_rural_comparison_three_line.csv"))
