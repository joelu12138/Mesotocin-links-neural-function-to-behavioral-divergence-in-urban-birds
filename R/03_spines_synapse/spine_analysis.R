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
  cat("已写出:", path, "\n")
}

w <- reshape(d, idvar = c("Subject", "Group", "Region"), timevar = "Metric", direction = "wide")
w$total <- w$Exp.m + w$Exp.s + w$Exp.t

## ---- 表 1: 三因素混合方差分析 (GG 校正) ----
m3 <- aov_car(Exp ~ Group * Region * Metric + Error(Subject / (Region * Metric)), data = d)
at <- m3$anova_table
aov_tab <- data.frame(
  效应 = rownames(at),
  F = round(at[["F"]], 3),
  分子df = round(at[["num Df"]], 2),
  分母df = round(at[["den Df"]], 2),
  p = fmtp(at[["Pr(>F)"]]),
  显著性 = sig(at[["Pr(>F)"]]),
  check.names = FALSE)
cat("## 表 1: 三因素混合方差分析 (GG 校正) ##\n")
print(aov_tab, row.names = FALSE)
wr(aov_tab, "data/树突棘_三因素混合方差分析_三线表.csv")

## ---- 表 2: 多重比较 (emmeans 模型合并误差 + Bonferroni) ----
desc_cells <- function(df, cell_cols, value_col) {
  cells <- unique(df[, cell_cols, drop = FALSE])
  cells <- cells[do.call(order, cells), , drop = FALSE]
  res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    mask <- rep(TRUE, nrow(df))
    for (cc in cell_cols) mask <- mask & (df[[cc]] == cells[[cc]][i])
    xu <- df[[value_col]][mask & df$Group == "Urban"]
    xr <- df[[value_col]][mask & df$Group == "Rural"]
    data.frame(Urban = sprintf("%.3f ± %.3f", mean(xu), sd(xu) / sqrt(length(xu))),
               Rural = sprintf("%.3f ± %.3f", mean(xr), sd(xr) / sqrt(length(xr))),
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
    显著性 = sig(m$p_Bonf_raw), check.names = FALSE)
  cbind(cell_part, out)[c(cn, "Urban", "Rural", "t", "df", "p", "p_Bonf", "显著性")]
}

# m/s/t 密度（6 次比较一个家族）
pr_d   <- as.data.frame(pairs(emmeans(m3, ~ Group | Region * Metric), adjust = "none"))
desc_d <- desc_cells(d, c("Region", "Metric"), "Exp")
tab_d  <- mk_tab(pr_d, desc_d, c("Region", "Metric"), c("脑区", "指标"))

# 总密度（2 次比较一个家族）
m_tot  <- aov_car(total ~ Group * Region + Error(Subject / Region), data = w)
pr_t   <- as.data.frame(pairs(emmeans(m_tot, ~ Group | Region), adjust = "none"))
desc_t <- desc_cells(w, "Region", "total")
tab_t  <- mk_tab(pr_t, desc_t, "Region", "脑区")
tab_t$指标 <- "总"
tab_t <- tab_t[, c("脑区", "指标", "Urban", "Rural", "t", "df", "p", "p_Bonf", "显著性")]

tab_d$指标 <- as.character(tab_d$指标)
tab_mc <- rbind(tab_d, tab_t)
tab_mc$指标 <- factor(tab_mc$指标, levels = c("m", "s", "t", "总"))
tab_mc <- tab_mc[order(tab_mc$脑区, tab_mc$指标), ]
rownames(tab_mc) <- NULL

cat("\n## 表 2: 多重比较 (emmeans 合并误差 + Bonferroni) ##\n")
print(tab_mc, row.names = FALSE)
wr(tab_mc, "data/树突棘_多重比较_三线表.csv")
