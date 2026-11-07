# Three-way repeated-measures ANOVA: Site x Sex (between-group) x Day (within-group) + GG correction
# Input: Relative exploration time.csv (Site, Group[=Sex], 1d..6d)
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(afex)
})

sig <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

wr <- function(x, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(x, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("Wrote:", path, "\n")
}

dat <- read.csv("Relative exploration time.csv", check.names = FALSE)
dat$ID <- factor(seq_len(nrow(dat)))

day_cols <- grep("^[0-9]d$", names(dat), value = TRUE)
day_cols <- day_cols[day_cols != "1d"]   # Day 1 = adaptation, excluded per Methods

long <- dat %>%
  rename(Sex = Group) %>%
  pivot_longer(all_of(day_cols), names_to = "Day", values_to = "Time") %>%
  mutate(
    Site   = factor(Site),
    Sex    = factor(Sex),
    DayNum = as.integer(gsub("[^0-9]", "", Day)),
    Day    = factor(DayNum, levels = sort(unique(DayNum)))
  ) %>%
  arrange(ID, DayNum)

res <- aov_ez(id = "ID", dv = "Time", data = long,
              between = c("Site", "Sex"), within = "Day", type = 3)
s   <- summary(res)
u   <- s$univariate.tests
adj <- s$pval.adjustments

eps_gg <- adj["Day", "GG eps"]          # GG epsilon for the within-group factor Day

eff  <- c("Site", "Sex", "Site:Sex", "Day", "Site:Day", "Sex:Day", "Site:Sex:Day")
lbl  <- c("Site (Urban vs Rural)", "Sex (Female vs Male)", "Site x Sex",
          "Day", "Site x Day", "Sex x Day", "Site x Sex x Day")
within <- grepl("Day", eff)             # terms containing Day require GG correction

tab <- data.frame(
  Effect = lbl,
  SS     = round(u[eff, "Sum Sq"]),
  df1    = u[eff, "num Df"],
  df2    = u[eff, "den Df"],
  MS     = round(u[eff, "Sum Sq"] / u[eff, "num Df"]),
  F      = round(u[eff, "F value"], 3),
  P      = u[eff, "Pr(>F)"],
  check.names = FALSE
)

# within-group terms: use GG-corrected P and df
tab$P[within] <- adj[eff[within], "Pr(>F[GG])"]
tab$GGdf1 <- NA_real_
tab$GGdf2 <- NA_real_
tab$GGdf1[within] <- u[eff[within], "num Df"] * eps_gg
tab$GGdf2[within] <- u[eff[within], "den Df"] * eps_gg

tab$Significance <- sig(tab$P)
tab$P <- format.pval(tab$P, digits = 3, eps = 0.001)

cat("===== Three-way repeated-measures ANOVA (Site x Sex x Day) =====\n")
cat(sprintf("between-group factors: Site, Sex   within-group factor: Day   individuals: %d\n", nrow(dat)))
cat(sprintf("Geisser-Greenhouse epsilon = %.4f\n\n", eps_gg))
print(tab, row.names = FALSE)

wr(tab, "three_way_rm_anova_results.csv")
