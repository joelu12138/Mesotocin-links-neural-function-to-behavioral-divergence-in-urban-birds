# Reproduce GraphPad two-way repeated-measures ANOVA + multiple comparisons (Bonferroni correction)
#
# Usage: Rscript rm_anova_graphpad.R <input.csv> [time_column]
#   Default time_column = total_exploration_time
# Input is a long table; required columns: Index (individual), Group (between-group factor), Day (within-group factor), <time> (measured value)
#
# Output two tables (ready for a three-line table):
#   1) Two-way RM ANOVA table: Type III SS/DF/MS/F/P, matching GraphPad
#      (Time main effect gets Geisser-Greenhouse correction; Interaction keeps the unadjusted df)
#   2) Multiple comparison table: per-day between-group comparisons (emmeans simple effects) + Bonferroni correction

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "Index.csv"
time_col    <- if (length(args) >= 2) args[2] else "total_exploration_time"

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(afex)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

run_graphpad_rm_anova <- function(input_file, time_col = "total_exploration_time") {
  dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

  long <- dat %>%
    rename(Time = all_of(time_col)) %>%
    mutate(
      Index  = factor(Index),
      Group  = factor(Group, levels = sort(unique(Group))),
      DayNum = as.integer(gsub("[^0-9]", "", Day)),
      Day    = factor(DayNum, levels = sort(unique(DayNum)))
    )

  g      <- levels(long$Group)
  k      <- nlevels(long$Day)
  n_subj <- nlevels(long$Index)

  res <- aov_ez(id = "Index", dv = "Time", data = long,
                between = "Group", within = "Day", type = 3)
  s   <- summary(res)
  u   <- s$univariate.tests
  adj <- s$pval.adjustments

  # ---- Type III SS / df / MS ----
  ss_int <- u["Group:Day", "Sum Sq"]; df_int <- u["Group:Day", "num Df"]
  ss_day <- u["Day", "Sum Sq"];       df_day <- u["Day", "num Df"]
  ss_grp <- u["Group", "Sum Sq"];     df_grp <- u["Group", "num Df"]
  ss_sub <- u["Group", "Error SS"];   df_sub <- u["Group", "den Df"]
  ss_res <- u["Day", "Error SS"];     df_res <- u["Day", "den Df"]

  ms_int <- ss_int / df_int; ms_day <- ss_day / df_day; ms_grp <- ss_grp / df_grp
  ms_sub <- ss_sub / df_sub; ms_res <- ss_res / df_res

  f_int <- u["Group:Day", "F value"]; p_int <- u["Group:Day", "Pr(>F)"]
  f_day <- u["Day", "F value"];       p_day <- u["Day", "Pr(>F)"]
  f_grp <- u["Group", "F value"];     p_grp <- u["Group", "Pr(>F)"]

  f_sub <- ms_sub / ms_res
  p_sub <- pf(f_sub, df_sub, df_res, lower.tail = FALSE)

  eps_gg      <- adj["Day", "GG eps"]
  df_day_gg_1 <- eps_gg * df_day
  df_day_gg_2 <- eps_gg * df_res
  p_day_gg    <- adj["Day", "Pr(>F[GG])"]

  ss_total <- sum((long$Time - mean(long$Time))^2)
  pct <- c(Interaction = ss_int, Time = ss_day, Group = ss_grp,
           Subject = ss_sub, Residual = ss_res) / ss_total * 100

  # ---- Table 1: two-way RM ANOVA table ----
  tab <- data.frame(
    Source = c("Group", "Time", "Group x Time", "Subject", "Residual"),
    SS     = round(c(ss_grp, ss_day, ss_int, ss_sub, ss_res)),
    DF     = c(df_grp, df_day, df_int, df_sub, df_res),
    MS     = round(c(ms_grp, ms_day, ms_int, ms_sub, ms_res)),
    F      = round(c(f_grp, f_day, f_int, f_sub, NA), 3),
    P      = c(p_grp, p_day_gg, p_int, p_sub, NA),
    check.names = FALSE
  )

  # ---- Table 2: multiple comparisons (per-day independent-samples Welch t-test + Bonferroni, matching GraphPad) ----
  days <- sort(unique(long$DayNum))
  mc_rows <- list()
  for (d in days) {
    dd <- long %>% filter(DayNum == d)
    x1 <- dd$Time[dd$Group == g[1]]; x2 <- dd$Time[dd$Group == g[2]]
    n1 <- length(x1); n2 <- length(x2)
    s1 <- sd(x1); s2 <- sd(x2)
    if (n1 < 2 || n2 < 2 || (s1 == 0 && s2 == 0)) next   # zero variance (e.g. Day1) cannot be tested
    m1 <- mean(x1); m2 <- mean(x2)
    diff <- m2 - m1                        # g[2] - g[1]
    se <- sqrt(s1^2 / n1 + s2^2 / n2)
    tval <- diff / se
    df <- (s1^2 / n1 + s2^2 / n2)^2 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
    p <- 2 * pt(-abs(tval), df)
    ci_lo <- diff - qt(0.975, df) * se
    ci_hi <- diff + qt(0.975, df) * se
    mc_rows[[length(mc_rows) + 1]] <- data.frame(
      DayNum = d, Mean1 = m1, Mean2 = m2, Diff = diff, SE = se,
      N1 = n1, N2 = n2, t = tval, df = df, p = p, CI_lo = ci_lo, CI_hi = ci_hi)
  }
  mc <- bind_rows(mc_rows)
  mc$P_Bonf <- p.adjust(mc$p, "bonferroni")

  # Daily Mean +/- SEM per group (raw time, for the three-line table)
  sem <- long %>%
    group_by(Group, DayNum) %>%
    summarise(MeanSEM = sprintf("%.3f +/- %.3f", mean(Time), sd(Time) / sqrt(n())),
              .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = MeanSEM)

  mc <- mc %>%
    left_join(sem, by = "DayNum") %>%
    select(DayNum, all_of(g), Diff, SE, t, df, p, P_Bonf, CI_lo, CI_hi)

  # ---- Output ----
  cat("\n================= Table 1: two-way RM ANOVA =================\n")
  cat("Input:", input_file, "  between-group:", paste(g, collapse = " vs "),
      "  within-group levels (k):", k, "  individuals:", n_subj, "\n\n")
  print(tab, row.names = FALSE)
  cat("\nGeisser-Greenhouse epsilon =", round(eps_gg, 4), "\n")

  cat("\n================= Table 2: multiple comparisons (Bonferroni) =================\n")
  print(mc, row.names = FALSE)

  # ---- Write CSV (UTF-8 BOM, to avoid garbled text in Excel) ----
  wr <- function(x, path) {
    con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
    write.table(x, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
                quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  }
  out1 <- sub("\\.csv$", "_anova_table.csv", input_file)
  out2 <- sub("\\.csv$", "_multiple_comparison.csv", input_file)
  wr(tab, out1)
  wr(mc, out2)
  cat("\nWrote:", out1, "\nWrote:", out2, "\n")

  invisible(list(anova = tab, multiple_comparison = mc, epsilon = eps_gg))
}

run_graphpad_rm_anova(input_file, time_col)
