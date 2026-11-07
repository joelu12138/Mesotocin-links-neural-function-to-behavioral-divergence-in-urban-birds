# 复现 GraphPad 双向重复测量 ANOVA + 多重比较（Bonferroni 校正）
#
# 用法: Rscript rm_anova_graphpad.R <input.csv> [time_column]
#   默认 time_column = total_exploration_time
# 输入为长表，必需列: Index(个体), Group(组间因子), Day(组内因子), <time>(测量值)
#
# 输出两张表（三线表可直接用）:
#   1) 双因素 RM ANOVA 表: Type III SS/DF/MS/F/P，对齐 GraphPad
#      (Time 主效应做 Geisser-Greenhouse 校正; Interaction 保留未校正 df)
#   2) 多重比较表: 逐日组间比较 (emmeans 简单效应) + Bonferroni 校正

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

  # 自动识别宽表(D1..D6)/长表(Day 列)
  if ("Day" %in% names(dat)) {
    long <- dat %>%
      rename(Time = all_of(time_col)) %>%
      mutate(DayNum = as.integer(gsub("[^0-9]", "", Day)))
  } else {
    day_cols <- grep("^D[0-9]+$", names(dat), value = TRUE)
    long <- dat %>%
      pivot_longer(all_of(day_cols), names_to = "Day", values_to = "Time") %>%
      mutate(DayNum = as.integer(gsub("[^0-9]", "", Day)))
  }

  long <- long %>%
    filter(DayNum != 1) %>%   # Day 1 = adaptation, excluded per Methods
    mutate(
      Index  = factor(Index),
      Group  = factor(Group, levels = sort(unique(Group))),
      Day    = factor(DayNum, levels = sort(unique(DayNum)))
    ) %>%
    arrange(Index, DayNum)

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

  # ---- 表 1: 双因素 RM ANOVA 表 ----
  tab <- data.frame(
    Source = c("Group", "Time", "Group × Time", "Subject", "Residual"),
    SS     = round(c(ss_grp, ss_day, ss_int, ss_sub, ss_res)),
    DF     = c(df_grp, df_day, df_int, df_sub, df_res),
    MS     = round(c(ms_grp, ms_day, ms_int, ms_sub, ms_res)),
    F      = round(c(f_grp, f_day, f_int, f_sub, NA), 3),
    P      = c(p_grp, p_day_gg, p_int, p_sub, NA),
    check.names = FALSE
  )

  # ---- 表 2: 多重比较（逐日独立样本 Welch t 检验 + Bonferroni，对齐 GraphPad）----
  days <- sort(unique(long$DayNum))
  mc_rows <- list()
  for (d in days) {
    dd <- long %>% filter(DayNum == d)
    x1 <- dd$Time[dd$Group == g[1]]; x2 <- dd$Time[dd$Group == g[2]]
    n1 <- length(x1); n2 <- length(x2)
    s1 <- sd(x1); s2 <- sd(x2)
    if (n1 < 2 || n2 < 2 || (s1 == 0 && s2 == 0)) next   # 零方差(如 Day1)无法检验
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

  # 各组每日 Mean ± SEM（原始时间，用于三线表）
  sem <- long %>%
    group_by(Group, DayNum) %>%
    summarise(MeanSEM = sprintf("%.1f ± %.1f", mean(Time), sd(Time) / sqrt(n())),
              .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = MeanSEM)

  mc <- mc %>%
    left_join(sem, by = "DayNum") %>%
    select(DayNum, all_of(g), Diff, SE, t, df, p, P_Bonf, CI_lo, CI_hi)

  # ---- 输出 ----
  cat("\n================= 表 1: 双因素 RM ANOVA =================\n")
  cat("输入:", input_file, "  组间:", paste(g, collapse = " vs "),
      "  组内水平(k):", k, "  个体:", n_subj, "\n\n")
  print(tab, row.names = FALSE)
  cat("\nGeisser-Greenhouse epsilon =", round(eps_gg, 4), "\n")

  cat("\n================= 表 2: 多重比较 (Bonferroni) =================\n")
  print(mc, row.names = FALSE)

  # ---- 写 CSV (UTF-8 BOM，避免 Excel 乱码) ----
  # C locale 下 write.table 的 connection 会把 ×/± 转成 <U+00D7> 或抹掉高位字节，
  # 这里按字节手工序列化，用 writeChar(useBytes=TRUE) 原样落盘。
  wr <- function(x, path) {
    q <- function(v) paste0('"', gsub('"', '""', as.character(v), fixed = TRUE), '"')
    fmt <- function(v) if (is.numeric(v)) { y <- as.character(v); y[is.na(v)] <- "NA"; y } else q(v)
    body <- vapply(seq_len(nrow(x)), function(i)
      paste(vapply(x, function(col) fmt(col[i]), character(1)), collapse = ","), character(1))
    txt <- c(paste(vapply(names(x), function(nm) q(nm), character(1)), collapse = ","), body)
    con <- file(path, "wb")
    writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
    writeChar(paste0(paste(txt, collapse = "\r\n"), "\r\n"), con, eos = NULL, useBytes = TRUE)
    close(con)
  }
  out1 <- sub("\\.csv$", "_anova_table.csv", input_file)
  out2 <- sub("\\.csv$", "_multiple_comparison.csv", input_file)
  wr(tab, out1)
  wr(mc, out2)
  cat("\n已写出:", out1, "\n已写出:", out2, "\n")

  invisible(list(anova = tab, multiple_comparison = mc, epsilon = eps_gg))
}

run_graphpad_rm_anova(input_file, time_col)
