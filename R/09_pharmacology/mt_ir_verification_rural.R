# Reproduce GraphPad two-way ANOVA (Type III SS) + Bonferroni multiple comparisons
# MT-ir_cell ~ Group x Region; key points: contr.sum contrasts + emmeans Bonferroni
suppressMessages({
  library(car)
  library(emmeans)
})

d <- read.csv("AA_HPNCL.csv", stringsAsFactors = FALSE, check.names = FALSE)
d$Group  <- factor(d$Group,  levels = c("Saline", "MT", "NMDA"))
d$Region <- factor(d$Region)

fit <- aov(`MT-ir_cell` ~ Group * Region, data = d,
           contrasts = list(Group = contr.sum, Region = contr.sum))

# ---- 1) Two-way ANOVA (Type III SS) ----
a3 <- car::Anova(fit, type = 3)
ss_int <- a3["Group:Region", "Sum Sq"]; df_int <- a3["Group:Region", "Df"]
ss_reg <- a3["Region", "Sum Sq"];        df_reg <- a3["Region", "Df"]
ss_grp <- a3["Group", "Sum Sq"];         df_grp <- a3["Group", "Df"]
ss_res <- a3["Residuals", "Sum Sq"];     df_res <- a3["Residuals", "Df"]
ms_int <- ss_int/df_int; ms_reg <- ss_reg/df_reg; ms_grp <- ss_grp/df_grp; ms_res <- ss_res/df_res
f_int  <- ms_int/ms_res; f_reg  <- ms_reg/ms_res;  f_grp  <- ms_grp/ms_res
p_int  <- pf(f_int, df_int, df_res, lower.tail = FALSE)
p_reg  <- pf(f_reg, df_reg, df_res, lower.tail = FALSE)
p_grp  <- pf(f_grp, df_grp, df_res, lower.tail = FALSE)
ss_total <- sum((d$`MT-ir_cell` - mean(d$`MT-ir_cell`))^2)
pct <- c(Interaction = ss_int, Row = ss_reg, Column = ss_grp, Residual = ss_res) / ss_total * 100

anova_tab <- data.frame(
  Source = c("Interaction", "Row Factor (Region)", "Column Factor (Group)", "Residual"),
  SS     = round(c(ss_int, ss_reg, ss_grp, ss_res), 1),
  DF     = c(df_int, df_reg, df_grp, df_res),
  MS     = round(c(ms_int, ms_reg, ms_grp, ms_res), 1),
  F      = round(c(f_int, f_reg, f_grp, NA), 4),
  P      = c(format.pval(c(p_int, p_reg, p_grp), digits = 4, eps = 0.0001), NA),
  Pct    = round(pct, 2),
  check.names = FALSE
)

# ---- 2) Bonferroni multiple comparisons (pairwise Group comparisons within each brain region, corrected for 3 comparisons) ----
n_grp <- setNames(as.numeric(table(d$Group, d$Region)[, "HP"]),
                  rownames(table(d$Group, d$Region)))  # n per group within each brain region = 5
cell  <- as.data.frame(emmeans(fit, ~ Group | Region))
pair  <- pairs(emmeans(fit, ~ Group | Region), adjust = "bonferroni")
tuk   <- as.data.frame(pair)
ci    <- as.data.frame(confint(pair))

parts <- strsplit(as.character(tuk$contrast), " - ")
gA <- sapply(parts, `[`, 1); gB <- sapply(parts, `[`, 2)
fmtp <- function(p) ifelse(p > 0.9999, ">0.9999",
                    ifelse(p < 0.0001, "<0.0001",
                           format(round(p, 4), nsmall = 4, scientific = FALSE, trim = TRUE)))

tuk_tab <- data.frame(
  Region     = tuk$Region,
  Comparison = tuk$contrast,
  Mean1  = signif(sapply(seq_len(nrow(tuk)), function(i)
             cell$emmean[cell$Group == gA[i] & cell$Region == tuk$Region[i]]), 4),
  Mean2  = signif(sapply(seq_len(nrow(tuk)), function(i)
             cell$emmean[cell$Group == gB[i] & cell$Region == tuk$Region[i]]), 4),
  Diff   = signif(tuk$estimate, 4),
  SE     = signif(tuk$SE, 4),
  N1     = n_grp[gA],
  N2     = n_grp[gB],
  t      = signif(abs(tuk$t.ratio), 4),
  df     = tuk$df,
  P_adj  = fmtp(tuk$p.value),
  CI_lo  = signif(ci$lower.CL, 4),
  CI_hi  = signif(ci$upper.CL, 4),
  check.names = FALSE
)

cat("===== 1) Two-way ANOVA (Type III SS) =====\n")
print(anova_tab, row.names = FALSE)
cat("\n===== 2) Bonferroni multiple comparisons =====\n")
print(tuk_tab, row.names = FALSE)

wr <- function(x, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(x, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("\nWrote:", path, "\n")
}
wr(anova_tab, "AA_HPNCL_two_way_anova.csv")
wr(tuk_tab,   "AA_HPNCL_bonferroni_posthoc.csv")
