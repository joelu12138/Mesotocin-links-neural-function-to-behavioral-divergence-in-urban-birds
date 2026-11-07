suppressMessages(library(openxlsx))
pfmt <- function(p) format.pval(p, digits = 3, eps = 0.001)
sem <- function(x) sd(x)/sqrt(length(x))
fmt_sem <- function(x, digits = 2) sprintf(paste0("%.", digits, "f +/- %.", digits, "f (n=%d)"), mean(x), sem(x), length(x))

## ================= Supplementary Table 1. Boldness =================
ad  <- read.csv("alerting_distance/Alerting distance.csv")
fid <- read.csv("alerting_distance/Flight initiation distance.csv")
b_test <- function(dat, label) {
  u <- na.omit(dat$Urban); r <- na.omit(dat$Rural)
  w <- t.test(u, r, var.equal = TRUE)
  c(label, fmt_sem(u), fmt_sem(r), sprintf("%.2f", mean(u) - mean(r)),
    round(w$statistic, 3), round(w$parameter, 2), pfmt(w$p.value))
}
bold <- as.data.frame(rbind(b_test(ad, "Alerting distance (m)"),
                            b_test(fid, "Flight initiation distance (m)")), stringsAsFactors = FALSE)
colnames(bold) <- c("Measure", "Urban (mean +/- SEM, n)", "Rural (mean +/- SEM, n)", "Difference", "t", "df", "P")

## ================= Supplementary Tables 2 & 3. Spatial cognition =================
re <- read.csv("relative_exploration_time/Relative exploration time.csv", stringsAsFactors = FALSE)
re <- re[!is.na(re$Index) & re$Index != "", ]; re$Group <- factor(re$Group)
times <- c("2d","3d","4d","5d","6d"); cols <- paste0("X", times)
Y <- as.matrix(re[, cols]); n <- nrow(re); k <- 5
long <- data.frame(y = as.vector(t(Y)),
                   Group = factor(rep(re$Group, each = k)),
                   Time = factor(rep(times, times = n), levels = times),
                   ID = factor(rep(1:n, each = k)))
fit <- aov(y ~ Group * Time + Error(ID / Time), data = long)
s <- summary(fit)
bet <- s[["Error: ID"]][[1]]; rownames(bet) <- trimws(rownames(bet))
win <- s[["Error: ID:Time"]][[1]]; rownames(win) <- trimws(rownames(win))
C <- contr.helmert(k); Q <- qr.Q(qr(C)); Z <- Y %*% Q
S <- matrix(0, k-1, k-1)
for (g in levels(re$Group)) { Zg <- Z[re$Group == g, , drop = FALSE]; S <- S + (nrow(Zg)-1)*cov(Zg) }
S <- S / (n - nlevels(re$Group))
GG <- sum(diag(S))^2 / ((k-1)*sum(S^2))
den_df <- win["Residuals","Df"]
sp_anova <- data.frame(
  Effect = c("Group (urban vs rural)", "Time (trial/day)", "Group x Time"),
  F = round(c(bet["Group","F value"], win["Time","F value"], win["Group:Time","F value"]), 3),
  num_df = c(bet["Group","Df"], win["Time","Df"], win["Group:Time","Df"]),
  den_df = c(bet["Residuals","Df"], den_df, den_df),
  P = c(pfmt(bet["Group","Pr(>F)"]), pfmt(win["Time","Pr(>F)"]), pfmt(win["Group:Time","Pr(>F)"])),
  GG_P = c("-", pfmt(pf(win["Time","F value"], 4*GG, den_df*GG, lower.tail=FALSE)),
                 pfmt(pf(win["Group:Time","F value"], 4*GG, den_df*GG, lower.tail=FALSE)))
)
sp_day <- do.call(rbind, lapply(times, function(tp) {
  u <- re[[paste0("X", tp)]][re$Group == "Urban"]; r <- re[[paste0("X", tp)]][re$Group == "Rural"]
  w <- t.test(u, r)
  data.frame(Day = paste0("Day ", substr(tp, 1, 1)),
             Urban = sprintf("%.3f +/- %.3f", mean(u), sem(u)),
             Rural = sprintf("%.3f +/- %.3f", mean(r), sem(r)),
             Difference = sprintf("%.3f", mean(u) - mean(r)),
             t = round(w$statistic, 3), df = round(w$parameter, 2),
             P = w$p.value, stringsAsFactors = FALSE)
}))
sp_day$FDR_P <- p.adjust(sp_day$P, "BH")
sp_day$P <- pfmt(sp_day$P); sp_day$FDR_P <- pfmt(sp_day$FDR_P)

## ================= Supplementary Table 6. Sociality =================
so <- read.csv("social_test/Total exploration time.csv"); so$Group <- factor(so$Group)
long_so <- data.frame(y = c(so$One, so$No), Group = rep(so$Group, 2),
                      Cond = factor(rep(c("One","No"), each = nrow(so))),
                      ID = factor(rep(1:nrow(so), 2)))
fit_so <- aov(y ~ Group * Cond + Error(ID / Cond), data = long_so)
s_so <- summary(fit_so)
win_so <- s_so[["Error: ID:Cond"]][[1]]; rownames(win_so) <- trimws(rownames(win_so))
so_anova <- data.frame(
  Effect = c("Condition (conspecific vs empty)", "Group x Condition"),
  F = round(c(win_so["Cond","F value"], win_so["Group:Cond","F value"]), 3),
  num_df = c(win_so["Cond","Df"], win_so["Group:Cond","Df"]),
  den_df = c(win_so["Residuals","Df"], win_so["Residuals","Df"]),
  P = c(pfmt(win_so["Cond","Pr(>F)"]), pfmt(win_so["Group:Cond","Pr(>F)"]))
)
so_pt <- do.call(rbind, lapply(c("Urban","Rural"), function(g) {
  s2 <- so[so$Group == g, ]; t <- t.test(s2$One, s2$No, paired = TRUE)
  data.frame(Group = g, Diff = sprintf("%.2f +/- %.2f", mean(s2$One - s2$No), sem(s2$One - s2$No)),
             t = round(t$statistic, 3), df = t$parameter, P = pfmt(t$p.value))
}))

## ================= Supplementary Table 7. Flocking =================
flock <- function(path) {
  d <- read.csv(path)
  if (!"Group" %in% names(d)) colnames(d)[1] <- "Group"
  d$Group <- factor(d$Group)
  long_d <- data.frame(y = c(d$One, d$Six), Group = rep(d$Group, 2),
                       Cond = factor(rep(c("One","Six"), each = nrow(d))),
                       ID = factor(rep(1:nrow(d), 2)))
  f <- aov(y ~ Group * Cond + Error(ID / Cond), data = long_d)
  s <- summary(f); b <- s[["Error: ID"]][[1]]; rownames(b) <- trimws(rownames(b))
  w <- s[["Error: ID:Cond"]][[1]]; rownames(w) <- trimws(rownames(w))
  anova <- data.frame(
    Effect = c("Group", "Condition (1 vs 6)", "Group x Condition"),
    F = round(c(b["Group","F value"], w["Cond","F value"], w["Group:Cond","F value"]), 3),
    num_df = c(b["Group","Df"], w["Cond","Df"], w["Group:Cond","Df"]),
    den_df = c(b["Residuals","Df"], w["Residuals","Df"], w["Residuals","Df"]),
    P = c(pfmt(b["Group","Pr(>F)"]), pfmt(w["Cond","Pr(>F)"]), pfmt(w["Group:Cond","Pr(>F)"]))
  )
  pt <- do.call(rbind, lapply(c("Urban","Rural"), function(g) {
    s2 <- d[d$Group == g, ]; t <- t.test(s2$One, s2$Six, paired = TRUE)
    data.frame(Group = g, Diff = sprintf("%.2f +/- %.2f", mean(s2$One - s2$Six), sem(s2$One - s2$Six)),
               t = round(t$statistic, 3), df = t$parameter, P = pfmt(t$p.value))
  }))
  list(anova = anova, pt = pt)
}
fl_t <- flock("flocking/AA.csv")
fl_f <- flock("flocking_frequency/AA.csv")

## ================= Supplementary Table 4. Sex effects on spatial learning =================
d4 <- read.csv("Relative exploration time.csv", check.names = FALSE)
colnames(d4)[1:2] <- c("Hab", "Sex")
d4$Hab <- factor(d4$Hab); d4$Sex <- factor(d4$Sex)
times4 <- c("2d","3d","4d","5d","6d")
Y4 <- as.matrix(d4[, times4]); n4 <- nrow(d4); k4 <- 5
long4 <- data.frame(y = as.vector(t(Y4)),
                    Hab = factor(rep(d4$Hab, each = k4)),
                    Sex = factor(rep(d4$Sex, each = k4)),
                    Time = factor(rep(times4, times = n4), levels = times4),
                    ID = factor(rep(1:n4, each = k4)))
fit4 <- aov(y ~ Hab * Sex * Time + Error(ID / Time), data = long4)
s4 <- summary(fit4)
bet4 <- s4[["Error: ID"]][[1]]; rownames(bet4) <- trimws(rownames(bet4))
win4 <- s4[["Error: ID:Time"]][[1]]; rownames(win4) <- trimws(rownames(win4))
C4 <- contr.helmert(k4); Q4 <- qr.Q(qr(C4)); Z4 <- Y4 %*% Q4
cells4 <- interaction(d4$Hab, d4$Sex)
S4 <- matrix(0, k4-1, k4-1)
for (cl in levels(cells4)) { Zg <- Z4[cells4 == cl, , drop = FALSE]; S4 <- S4 + (nrow(Zg)-1)*cov(Zg) }
S4 <- S4 / (n4 - length(levels(cells4)))
GG4 <- sum(diag(S4))^2 / ((k4-1)*sum(S4^2))
den4 <- win4["Residuals","Df"]
sex_tab <- data.frame(
  Effect = c("Sex", "Sex x Habitat", "Sex x Time", "Sex x Habitat x Time"),
  F = round(c(bet4["Sex","F value"], bet4["Hab:Sex","F value"],
              win4["Sex:Time","F value"], win4["Hab:Sex:Time","F value"]), 3),
  num_df = c(bet4["Sex","Df"], bet4["Hab:Sex","Df"], win4["Sex:Time","Df"], win4["Hab:Sex:Time","Df"]),
  den_df = c(bet4["Residuals","Df"], bet4["Residuals","Df"], den4, den4),
  P = c(pfmt(bet4["Sex","Pr(>F)"]), pfmt(bet4["Hab:Sex","Pr(>F)"]),
        pfmt(win4["Sex:Time","Pr(>F)"]), pfmt(win4["Hab:Sex:Time","Pr(>F)"])),
  GG_P = c("-", "-",
           pfmt(pf(win4["Sex:Time","F value"], win4["Sex:Time","Df"]*GG4, den4*GG4, lower.tail=FALSE)),
           pfmt(pf(win4["Hab:Sex:Time","F value"], win4["Hab:Sex:Time","Df"]*GG4, den4*GG4, lower.tail=FALSE)))
)
grid <- expand.grid(Site = c("Urban","Rural"), Day = times4, stringsAsFactors = FALSE)
sex_mc <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  site <- grid$Site[i]; tp <- grid$Day[i]
  f <- d4[[tp]][d4$Hab == site & d4$Sex == "Female"]
  m <- d4[[tp]][d4$Hab == site & d4$Sex == "Male"]
  w <- t.test(f, m)
  data.frame(Site = site, Day = paste0("Day ", substr(tp, 1, 1)),
             Female = sprintf("%.1f +/- %.1f (n=%d)", mean(f), sem(f), length(f)),
             Male = sprintf("%.1f +/- %.1f (n=%d)", mean(m), sem(m), length(m)),
             Difference = sprintf("%.1f", mean(f) - mean(m)),
             t = round(w$statistic, 3), df = round(w$parameter, 2),
             P = w$p.value, stringsAsFactors = FALSE)
}))
sex_mc$FDR_P <- p.adjust(sex_mc$P, "BH")
sex_mc$P <- pfmt(sex_mc$P); sex_mc$FDR_P <- pfmt(sex_mc$FDR_P)

## ================= Write Excel =================
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", border = "TopBottom", borderStyle = "thin")
tt <- createStyle(fontSize = 12, textDecoration = "bold")
lg <- createStyle(fontSize = 9, fontColour = "#333333")

add_block <- function(wb, sheet, start, title, tbl, legend) {
  writeData(wb, sheet, title, startRow = start, startCol = 1)
  addStyle(wb, sheet, tt, rows = start, cols = 1)
  hr <- start + 2
  writeData(wb, sheet, tbl, startRow = hr, startCol = 1, colNames = TRUE)
  addStyle(wb, sheet, hs, rows = hr, cols = 1:ncol(tbl))
  lr <- hr + nrow(tbl)
  addStyle(wb, sheet, createStyle(border = "bottom", borderStyle = "thin"), rows = lr, cols = 1:ncol(tbl))
  writeData(wb, sheet, legend, startRow = lr + 2, startCol = 1)
  addStyle(wb, sheet, lg, rows = lr + 2, cols = 1)
  invisible(lr + 3)
}

addWorksheet(wb, "Table1_Boldness")
add_block(wb, "Table1_Boldness", 1,
  "Supplementary Table 1 | Alerting and flight initiation distances (m) of urban and rural tree sparrows.",
  bold,
  "Data are mean +/- SEM (n). Difference = Urban - Rural. P values are two-sided from independent-samples Student's t-tests (equal variances). *P < 0.05, **P < 0.01, ***P < 0.001.")

addWorksheet(wb, "Table2_SpatialCognition")
add_block(wb, "Table2_SpatialCognition", 1,
  "Supplementary Table 2 | Repeated-measures ANOVA of relative exploration time in the spatial-foraging task.",
  sp_anova,
  "Day 1 is the baseline and was excluded. GG_P, Greenhouse-Geisser corrected P (sphericity violated). *P < 0.05, **P < 0.01, ***P < 0.001.")

addWorksheet(wb, "Table3_SpatialCognition_PerDay")
add_block(wb, "Table3_SpatialCognition_PerDay", 1,
  "Supplementary Table 3 | Between-population comparison of relative exploration time at each trial (day).",
  sp_day,
  "Values are mean +/- SEM (Urban n = 23, Rural n = 26). Difference = Urban - Rural. P values are two-sided from Welch's t-tests; FDR_P, Benjamini-Hochberg adjusted P. More negative values indicate faster food localization relative to baseline.")

addWorksheet(wb, "Table4_SexSpatial")
add_block(wb, "Table4_SexSpatial", 1,
  "Supplementary Table 4 | Effects of sex on relative exploration time in the spatial-foraging task.",
  sex_tab,
  "Three-way repeated-measures ANOVA (Habitat x Sex x Time; Habitat = urban/rural, Sex = female/male, Time = 5 trials, day-1 baseline excluded). GG_P, Greenhouse-Geisser corrected P. Sample sizes: 4 urban females, 19 urban males, 10 rural females, 17 rural males. No significant sex effects were detected.")
add_block(wb, "Table4_SexSpatial", 10,
  "Between-sex comparison at each site and trial (Welch's t-test, FDR-adjusted).",
  sex_mc,
  "Female vs Male. Values are mean +/- SEM (n). Difference = Female - Male. FDR_P, Benjamini-Hochberg adjusted P. No significant sex differences were detected.")

addWorksheet(wb, "Table6_Sociality")
add_block(wb, "Table6_Sociality", 1,
  "Supplementary Table 6 | Sociality test: exploration time (s) near a conspecific (One) vs an empty compartment (No).",
  so_anova,
  "Repeated-measures ANOVA (Group x Condition). Total exploration time was fixed at 300 s per bird, so the Group main effect is degenerate; the between-population difference is given by the Group x Condition interaction.")
add_block(wb, "Table6_Sociality", 9,
  "Paired t-test (One vs No) within each population.",
  so_pt,
  "Diff = One - No (s), mean +/- SEM. *P < 0.05, **P < 0.01, ***P < 0.001.")

addWorksheet(wb, "Table7_Flocking")
add_block(wb, "Table7_Flocking", 1,
  "Supplementary Table 7 | Flocking test: exploration time (s) and frequency near 1 vs 6 conspecifics.",
  fl_t$anova,
  "Repeated-measures ANOVA (Group x Condition). Exploration time (s).")
add_block(wb, "Table7_Flocking", 8,
  "Paired t-test (One vs Six), time (s).",
  fl_t$pt,
  "Diff = One - Six (s), mean +/- SEM.")
add_block(wb, "Table7_Flocking", 14,
  "Repeated-measures ANOVA (Group x Condition). Exploration frequency (counts).",
  fl_f$anova,
  "Exploration frequency (counts).")
add_block(wb, "Table7_Flocking", 21,
  "Paired t-test (One vs Six), frequency (counts).",
  fl_f$pt,
  "Diff = One - Six (frequency), mean +/- SEM. *P < 0.05, **P < 0.01, ***P < 0.001.")

saveWorkbook(wb, "Supplementary_Tables.xlsx", overwrite = TRUE)
cat("Wrote: Supplementary_Tables.xlsx (6 sheets: Table 1, 2, 3, 4, 6, 7)\n")
