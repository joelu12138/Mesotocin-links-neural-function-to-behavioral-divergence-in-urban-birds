suppressMessages({library(afex); library(emmeans)})
options(digits = 6)

d <- read.csv("data/data_long.csv", stringsAsFactors = FALSE)
d$Group  <- factor(d$Group,  levels = c("Urban", "Rural"))
d$Region <- factor(d$Region, levels = c("HP", "NCL"))
d$Metric <- factor(d$Metric, levels = c("m", "s", "t"))

d$Subject <- NA
for (g in levels(d$Group)) {
  gi <- which(d$Group == g)
  d$Subject[gi] <- rep(paste0(substr(g, 1, 1), seq_len(length(gi) / 6)), each = 6)
}
d$Subject <- factor(d$Subject)

sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
fmtp <- function(p) format.pval(p, digits = 3, eps = 0.001)
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}

w <- reshape(d, idvar = c("Subject", "Group", "Region"), timevar = "Metric", direction = "wide")
w$total <- w$Exp.m + w$Exp.s + w$Exp.t

## ---- Table 1: three-way mixed ANOVA (GG correction) ----
m3 <- aov_car(Exp ~ Group * Region * Metric + Error(Subject / (Region * Metric)), data = d)
at <- m3$anova_table
aov_tab <- data.frame(
  Effect = rownames(at),
  F = round(at[["F"]], 3),
  df1 = round(at[["num Df"]], 2),
  df2 = round(at[["den Df"]], 2),
  p = fmtp(at[["Pr(>F)"]]),
  Significance = sig(at[["Pr(>F)"]]),
  check.names = FALSE)
cat("## Table 1: three-way mixed ANOVA (GG correction) ##\n")
print(aov_tab, row.names = FALSE)
wr(aov_tab, "data/spine_three_way_mixed_anova_three_line.csv")

## ---- Table 2: multiple comparisons (emmeans pooled error + Bonferroni) ----
desc_cells <- function(df, cell_cols, value_col) {
  cells <- unique(df[, cell_cols, drop = FALSE])
  cells <- cells[do.call(order, cells), , drop = FALSE]
  res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    mask <- rep(TRUE, nrow(df))
    for (cc in cell_cols) mask <- mask & (df[[cc]] == cells[[cc]][i])
    xu <- df[[value_col]][mask & df$Group == "Urban"]
    xr <- df[[value_col]][mask & df$Group == "Rural"]
    data.frame(Urban = sprintf("%.3f +/- %.3f", mean(xu), sd(xu) / sqrt(length(xu))),
               Rural = sprintf("%.3f +/- %.3f", mean(xr), sd(xr) / sqrt(length(xr))),
               check.names = FALSE)
  }))
  cbind(cells, res)
}

mk_tab <- function(pr, desc, cell_cols, cn) {
  pr$p_Bonf_raw <- p.adjust(pr$p.value, "bonferroni")
  m <- merge(desc, pr[, c(cell_cols, "t.ratio", "df", "p.value", "p_Bonf_raw")],
             by = cell_cols, sort = FALSE)
  cell_part <- m[, cell_cols, drop = FALSE]
  names(cell_part) <- cn
  out <- data.frame(
    Urban = m$Urban, Rural = m$Rural,
    t = round(m$t.ratio, 3), df = round(m$df, 2),
    p = fmtp(m$p.value), p_Bonf = fmtp(m$p_Bonf_raw),
    Significance = sig(m$p_Bonf_raw), check.names = FALSE)
  cbind(cell_part, out)[c(cn, "Urban", "Rural", "t", "df", "p", "p_Bonf", "Significance")]
}

# m/s/t density (6 comparisons per family)
pr_d   <- as.data.frame(pairs(emmeans(m3, ~ Group | Region * Metric), adjust = "none"))
desc_d <- desc_cells(d, c("Region", "Metric"), "Exp")
tab_d  <- mk_tab(pr_d, desc_d, c("Region", "Metric"), c("brain_region", "Metric"))

# total density (2 comparisons per family)
m_tot  <- aov_car(total ~ Group * Region + Error(Subject / Region), data = w)
pr_t   <- as.data.frame(pairs(emmeans(m_tot, ~ Group | Region), adjust = "none"))
desc_t <- desc_cells(w, "Region", "total")
tab_t  <- mk_tab(pr_t, desc_t, "Region", "brain_region")
tab_t$Metric <- "Total"
tab_t <- tab_t[, c("brain_region", "Metric", "Urban", "Rural", "t", "df", "p", "p_Bonf", "Significance")]

tab_d$Metric <- as.character(tab_d$Metric)
tab_mc <- rbind(tab_d, tab_t)
tab_mc$Metric <- factor(tab_mc$Metric, levels = c("m", "s", "t", "Total"))
tab_mc <- tab_mc[order(tab_mc$brain_region, tab_mc$Metric), ]
rownames(tab_mc) <- NULL

cat("\n## Table 2: multiple comparisons (emmeans pooled error + Bonferroni) ##\n")
print(tab_mc, row.names = FALSE)
wr(tab_mc, "data/spine_multiple_comparisons_three_line.csv")
