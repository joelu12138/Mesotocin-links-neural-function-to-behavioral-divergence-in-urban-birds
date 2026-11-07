# 城乡树麻雀脑组织不同电流强度下放电峰值数
# 双因素重复度量方差分析 (Two-way RM ANOVA, 不假设球形 → Geisser-Greenhouse 校正)
#   + Bonferroni 多重比较 (复现 GraphPad Prism)

raw <- read.csv("data/AA.csv", check.names = FALSE)
labels <- c("0A","20A","40A","60A","80A","100A","120A","140A","160A","180A","200A")
d <- raw[, c("Group", labels)]
d$Group <- factor(d$Group)
Y <- as.matrix(d[, labels])
n <- nrow(d); k <- length(labels)

sem <- function(x) sd(x)/sqrt(length(x))
sig  <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                        ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))
pfmt <- function(p) ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))
fmt  <- function(x, d = 4) sprintf("%g", signif(x, d))
# 4 位有效数字 + 四舍五入(向远离 0)，保留尾随 0，对齐 Prism
fmt4 <- function(x) {
  if (x == 0) return("0.0000")
  ni <- floor(log10(abs(x))) + 1
  dd <- max(0, 4 - ni)
  v  <- sign(x) * trunc(abs(x) * 10^dd + 0.5) / 10^dd
  sprintf(paste0("%.", dd, "f"), v)
}

# ---- 长格式 + aov (混合模型: 组间 Group × 组内 Current) ----
long <- data.frame(Subject = factor(rep(1:n, each = k)),
                   Group   = factor(rep(d$Group, each = k)),
                   Current = factor(rep(labels, times = n), levels = labels),
                   y = as.vector(t(Y)))
fit <- aov(y ~ Group * Current + Error(Subject / Current), data = long)
s   <- summary(fit)
bet <- s[["Error: Subject"]][[1]];         rownames(bet) <- trimws(rownames(bet))
win <- s[["Error: Subject:Current"]][[1]]; rownames(win) <- trimws(rownames(win))

# ---- 提取 SS / DF / MS ----
SS_A  <- bet["Group","Sum Sq"];        MS_A  <- bet["Group","Mean Sq"]
SS_S  <- bet["Residuals","Sum Sq"];    MS_S  <- bet["Residuals","Mean Sq"]
SS_B  <- win["Current","Sum Sq"];      MS_B  <- win["Current","Mean Sq"]
SS_AB <- win["Group:Current","Sum Sq"];MS_AB <- win["Group:Current","Mean Sq"]
SS_R  <- win["Residuals","Sum Sq"];    MS_R  <- win["Residuals","Mean Sq"]
dfd_B <- win["Residuals","Df"]

# ---- F 值与自由度 ----
F_A  <- MS_A / MS_S    # Column Factor, df (1, 14)
F_B  <- MS_B / MS_R    # Time,           df (10, 140)
F_AB <- MS_AB / MS_R   # Interaction,    df (10, 140)
F_S  <- MS_S / MS_R    # Subject,        df (14, 140)

p_A  <- pf(F_A,  bet["Group","Df"],       bet["Residuals","Df"], lower.tail = FALSE)
p_B  <- pf(F_B,  win["Current","Df"],     dfd_B, lower.tail = FALSE)
p_AB <- pf(F_AB, win["Group:Current","Df"],dfd_B, lower.tail = FALSE)
p_S  <- pf(F_S,  bet["Residuals","Df"],   dfd_B, lower.tail = FALSE)

# ---- Geisser-Greenhouse epsilon (组内因子 11 水平) ----
C <- contr.helmert(k); Q <- qr.Q(qr(C)); Z <- Y %*% Q
Smat <- matrix(0, k - 1, k - 1)
for (g in levels(d$Group)) { Zg <- Z[d$Group == g, , drop = FALSE]; Smat <- Smat + (nrow(Zg) - 1) * cov(Zg) }
Smat <- Smat / (n - nlevels(d$Group))
GG <- sum(diag(Smat))^2 / ((k - 1) * sum(Smat^2))

dfn_B <- win["Current","Df"]
dfn_Bgg <- dfn_B * GG; dfd_Bgg <- dfd_B * GG
p_Bgg <- pf(F_B, dfn_Bgg, dfd_Bgg, lower.tail = FALSE)

SS_tot <- SS_A + SS_S + SS_B + SS_AB + SS_R

cat("Two-way RM ANOVA    Matching: Stacked\n")
cat("Assume sphericity?    No\nAlpha    0.05\n\n")

# ---- Source of Variation 表 ----
src <- data.frame(
  `Source of Variation` = c("Interaction", "Time", "Column Factor", "Subject"),
  `% of total variation` = fmt(c(SS_AB, SS_B, SS_A, SS_S) / SS_tot * 100),
  `P value`              = pfmt(c(p_AB, p_Bgg, p_A, p_S)),
  `P value summary`      = sig(c(p_AB, p_Bgg, p_A, p_S)),
  `Significant?`         = ifelse(c(p_AB, p_Bgg, p_A, p_S) < 0.05, "Yes", "No"),
  `Geisser-Greenhouse's epsilon` = c("", sprintf("%.4f", GG), "", ""),
  check.names = FALSE)
cat("Source of Variation    % of total variation    P value    P value summary    Significant?    Geisser-Greenhouse's epsilon\n")
print(src, row.names = FALSE)

# ---- ANOVA table ----
anv <- data.frame(
  `ANOVA table` = c("Interaction", "Time", "Column Factor", "Subject", "Residual"),
  SS  = fmt(c(SS_AB, SS_B, SS_A, SS_S, SS_R)),
  DF  = c(win["Group:Current","Df"], win["Current","Df"], bet["Group","Df"],
          bet["Residuals","Df"], dfd_B),
  MS  = fmt(c(MS_AB, MS_B, MS_A, MS_S, MS_R)),
  `F (DFn, DFd)` = c(sprintf("F (%d, %d) = %.2f", as.integer(win["Group:Current","Df"]), as.integer(dfd_B), F_AB),
                     sprintf("F (%.3f, %.2f) = %.1f", dfn_Bgg, dfd_Bgg, F_B),
                     sprintf("F (%d, %d) = %.2f",   as.integer(bet["Group","Df"]), as.integer(bet["Residuals","Df"]), F_A),
                     sprintf("F (%d, %d) = %.2f", as.integer(bet["Residuals","Df"]), as.integer(dfd_B), F_S),
                     ""),
  `P value` = c(sprintf("P%s", pfmt(c(p_AB, p_Bgg, p_A, p_S))), ""),
  check.names = FALSE)
cat("\nANOVA table    SS    DF    MS    F (DFn, DFd)    P value\n")
print(anv, row.names = FALSE)

# ---- Bonferroni 多重比较: 每行(电流)内 R-Bird - U-Bird ----
# HP 与 NCL 分别运行（换输入文件即可，Methods 要求两脑区独立分析）。
cmp_idx <- 2:k                       # 跳过 0 pA（无诱发放电）
m_cmp   <- length(cmp_idx)           # Bonferroni 比较次数由数据推导
alpha   <- 0.05
nU <- sum(d$Group == "Urban"); nR <- sum(d$Group == "Rural")
SE_pool <- sqrt(MS_R * (1/nU + 1/nR))   # 合并 ANOVA 残差 MSE 的 SE（GraphPad Bonferroni）

bonf <- do.call(rbind, lapply(cmp_idx, function(i) {
  cur <- labels[i]
  u  <- d[d$Group == "Urban", cur]
  rv <- d[d$Group == "Rural", cur]
  md <- mean(rv) - mean(u)                    # R - U
  tcrit <- qt(0.975, dfd_B)
  lo <- md - tcrit * SE_pool                  # 合并 MSE 的 95% CI（与 p 一致）
  hi <- md + tcrit * SE_pool
  t_p   <- md / SE_pool                        # 合并 SE 的 t
  p_raw <- 2 * pt(-abs(t_p), dfd_B)
  p_bon <- min(1, p_raw * m_cmp)
  data.frame(Row = sprintf("Row %d", i),
             `Mean Diff.` = fmt4(md),
             `95.00% CI of diff.` = paste0(fmt4(lo), " to ", fmt4(hi)),
             `Significant?` = ifelse(p_bon < alpha, "Yes", "No"),
             check.names = FALSE)
}))

cat("\nCompare each cell mean with the other cell mean in that row\n")
cat(sprintf("Number of families: 1\nNumber of comparisons per family: %d\nAlpha: %.2f\n", m_cmp, alpha))
cat("\nBonferroni's multiple comparisons test    Mean Diff.    95.00% CI of diff.    Significant?\n")
cat("  R-Bird - U-Bird\n")
print(bonf, row.names = FALSE)

# ---- 写出 CSV (UTF-8 BOM) ----
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}
wr(src,  "data/放电峰值_重复度量方差分析_变异来源表.csv")
wr(anv,  "data/放电峰值_重复度量方差分析_ANOVA表.csv")
wr(bonf, "data/放电峰值_Bonferroni多重比较.csv")
