options(digits = 6)

d <- read.csv("data/AA.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Group) & d$Group != "", ]
d$Group <- factor(d$Group)
sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

d$PI <- (d$Six - d$One) / (d$Six + d$One)   # Preference index, >0 = prefers the 6-bird cage

## ---- 1. Within-group paired comparison: One vs Six (paired Wilcoxon signed-rank) ----
cat("## 1. Within-group paired comparison (One vs Six, paired Wilcoxon signed-rank) ##\n")
pair <- do.call(rbind, lapply(c("Overall", levels(d$Group)), function(g) {
  if (g == "Overall") { x <- d$One; y <- d$Six } else { x <- d$One[d$Group == g]; y <- d$Six[d$Group == g] }
  wt <- wilcox.test(x, y, paired = TRUE)
  data.frame(Group = g, n = length(x),
             OneMedian = median(x), SixMedian = median(y),
             V = unname(wt$statistic),
             p = format.pval(wt$p.value, digits = 3, eps = 0.001),
             Significance = sig(wt$p.value), check.names = FALSE)
}))
print(pair, row.names = FALSE)

## ---- 2. Between-group comparison: Mann-Whitney U ----
cat("\n## 2. Between-group comparison (Urban vs Rural, Mann-Whitney U) ##\n")
comp <- do.call(rbind, lapply(
  list(c("One frequency", "One", "%.1f"),
       c("Six frequency", "Six", "%.1f"),
       c("Total frequency", "Tot", "%.1f"),
       c("Preference index PI", "PI", "%.3f")),
  function(spec) {
    lab <- spec[1]; v <- spec[2]; fmt <- spec[3]
    x <- if (v == "Tot") d$One + d$Six else d[[v]]
    u <- x[d$Group == "Urban"]; r <- x[d$Group == "Rural"]
    wt <- wilcox.test(u, r)
    data.frame(Metric = lab,
               Urban = sprintf(paste0(fmt, "+/-", fmt), mean(u), sem(u)),
               Rural = sprintf(paste0(fmt, "+/-", fmt), mean(r), sem(r)),
               W = unname(wt$statistic),
               p = format.pval(wt$p.value, digits = 3, eps = 0.001),
               Significance = sig(wt$p.value), check.names = FALSE)
  }
))
print(comp, row.names = FALSE)

## ---- Write CSV ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}
wr(pair, "data/flock_size_paired_wilcoxon.csv")
wr(comp, "data/flock_size_between_group_mann_whitney_u.csv")
