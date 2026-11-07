# !!! SUPERSEDED !!!
# This script (readout: Relative_exploration_time) is obsolete. The current
# analysis uses total_exploration_time (latest version) in
# spatial_cognition_rm_anova.R. Kept here only for provenance; do not use for
# the reported results.
#
# Repeated-measures ANOVA (Group x Time) + per-timepoint Welch t-tests, Bonferroni.
# Same algorithm as rm_anova.R / rm_anova_table.R.
# Plain-loop style (no lapply/rbind/sprintf) to dodge intermittent Rscript
# native segfault on this Windows box.
options(contrasts = c("contr.sum", "contr.poly"))
options(digits = 6)

input_file <- "Index.csv"   # EDIT: relative path (repo-external absolute path removed)
d <- read.csv(input_file, stringsAsFactors = FALSE)
d$Day <- factor(d$Day, levels = c("1d","2d","3d","4d","5d","6d"))

wide <- reshape(d, idvar = c("Index", "Group"), timevar = "Day",
                direction = "wide", v.names = "Relative_exploration_time")
names(wide)[-(1:2)] <- paste0("X", 1:6, "d")
wide$Group <- factor(wide$Group)
times <- c("2d","3d","4d","5d","6d")

# ---- RM-ANOVA ----
mlm <- lm(cbind(X2d,X3d,X4d,X5d,X6d) ~ Group, data = wide)
idata <- data.frame(Time = factor(1:5, labels = times))
aov.out <- car::Anova(mlm, idata = idata, idesign = ~Time, type = "III")
s <- summary(aov.out, multivariate = FALSE)
print(s)

# ---- three-line table (Group / Time / Group:Time), GG-corrected p ----
ut <- s$univariate.tests
gg <- s$pval.adjustments
mauchly <- s$sphericity.tests
eff <- c("Group", "Time", "Group:Time")
Fv <- numeric(0); ndf <- integer(0); ddf <- integer(0)
praw <- character(0); pggv <- character(0); ssig <- character(0)
for (nm in eff) {
  Fv  <- c(Fv,  as.numeric(ut[nm, "F value"]))
  ndf <- c(ndf, as.integer(ut[nm, "num Df"]))
  ddf <- c(ddf, as.integer(ut[nm, "den Df"]))
  pr  <- as.numeric(ut[nm, "Pr(>F)"])
  praw <- c(praw, format.pval(pr, digits = 3, eps = 0.001))
  pg  <- if (nm %in% rownames(gg)) as.numeric(gg[nm, "Pr(>F[GG])"]) else NA
  pggv <- c(pggv, ifelse(is.na(pg), "-", format.pval(pg, digits = 3, eps = 0.001)))
  puse <- ifelse(is.na(pg), pr, pg)
  ssig <- c(ssig, ifelse(puse < 0.001, "***", ifelse(puse < 0.01, "**",
                    ifelse(puse < 0.05, "*", "ns"))))
}
tab <- data.frame(effect = eff, F = round(Fv, 3), numDf = ndf, denDf = ddf,
                  p = praw, pGG = pggv, sig = ssig,
                  stringsAsFactors = FALSE, check.names = FALSE)
print(tab, row.names = FALSE)
cat("Mauchly: W =", round(mauchly["Time", "Test statistic"], 4),
    ", p =", format.pval(mauchly["Time", "p-value"], digits = 3, eps = 0.001), "\n")

# ---- per-timepoint Welch t-tests, Bonferroni (Methods: day-specific contrasts) ----
glev <- levels(wide$Group)
day <- character(0); tval <- numeric(0); dfval <- numeric(0); pval <- numeric(0)
m1 <- numeric(0); s1 <- numeric(0); n1 <- integer(0)
m2 <- numeric(0); s2 <- numeric(0); n2 <- integer(0)
for (tp in times) {
  col <- paste0("X", tp)
  a <- wide[[col]][wide$Group == glev[1]]
  b <- wide[[col]][wide$Group == glev[2]]
  w <- t.test(a, b)
  day  <- c(day, tp)
  tval <- c(tval, as.numeric(w$statistic))
  dfval<- c(dfval, as.numeric(w$parameter))
  pval <- c(pval, as.numeric(w$p.value))
  m1 <- c(m1, mean(a)); s1 <- c(s1, sd(a)/sqrt(length(a))); n1 <- c(n1, length(a))
  m2 <- c(m2, mean(b)); s2 <- c(s2, sd(b)/sqrt(length(b))); n2 <- c(n2, length(b))
}
pbon <- p.adjust(pval, method = "bonferroni")
sig <- ifelse(pbon < 0.001, "***", ifelse(pbon < 0.01, "**", ifelse(pbon < 0.05, "*", "ns")))
t2 <- data.frame(day = day,
                 g1 = sprintf("%.3f+-%.3f (n=%d)", m1, s1, n1),
                 g2 = sprintf("%.3f+-%.3f (n=%d)", m2, s2, n2),
                 t = round(tval, 3), df = round(dfval, 2),
                 p = signif(pval, 6), P_Bonf = signif(pbon, 6), sig = sig,
                 stringsAsFactors = FALSE, check.names = FALSE)
names(t2)[2:3] <- glev
print(t2, row.names = FALSE)

out1 <- sub("\\.csv$", "_rm_anova_three_line_table.csv", input_file)
out2 <- sub("\\.csv$", "_timepoint_between_group_comparison.csv", input_file)
write.csv(tab, out1, row.names = FALSE)
write.csv(t2, out2, row.names = FALSE)
cat("wrote", out1, ",", out2, "\n")
