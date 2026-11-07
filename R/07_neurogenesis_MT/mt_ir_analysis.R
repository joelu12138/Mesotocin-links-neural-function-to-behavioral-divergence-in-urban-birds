sem <- function(x) sd(x)/sqrt(length(x))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
wr <- function(tab, path) {
  con <- file(path, "w", encoding = "UTF-8"); writeLines("\ufeff", con, sep = ""); close(con)
  write.table(tab, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = TRUE, append = TRUE, fileEncoding = "UTF-8")
  cat("已写出:", path, "\n")
}

## ================= AA_HPNCL.csv：双因素方差分析 + Fisher LSD =================
d1 <- read.csv("data/AA_HPNCL.csv", stringsAsFactors = FALSE)
y <- setdiff(names(d1), c("Group", "Region"))[1]
names(d1)[names(d1) == y] <- "MT"
d1$Group <- factor(d1$Group); d1$Region <- factor(d1$Region)

cat("========== AA_HPNCL (MT-ir_cell) ==========\n")
cat("## 描述统计 (Mean±SEM) ##\n")
for (r in c("HP", "NCL")) for (g in c("Urban", "Rural")) {
  v <- d1[d1$Region == r & d1$Group == g, "MT"]
  cat(sprintf("%s %s: %.3f ± %.3f (n=%d)\n", g, r, mean(v), sem(v), length(v)))
}

fit <- aov(MT ~ Group * Region, data = d1)
a <- summary(fit)[[1]]
rownames(a) <- trimws(rownames(a))
cat("\n## 双因素方差分析 ##\n"); print(a)

G <- c("Group", "Region", "Group:Region")
anova_tab <- data.frame(
  效应 = c("Group (城乡)", "Region (脑区)", "Group × Region (交互)"),
  F = round(a[G, "F value"], 3),
  分子df = a[G, "Df"], 分母df = a["Residuals", "Df"],
  p = format.pval(a[G, "Pr(>F)"], digits = 3, eps = 0.001),
  显著性 = sig(a[G, "Pr(>F)"]), check.names = FALSE)
cat("\n## 双因素方差分析三线表 ##\n"); print(anova_tab, row.names = FALSE)

# Fisher LSD: 同一脑区城乡比较
MSE <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
lsd_line <- function(r) {
  su <- d1[d1$Region == r & d1$Group == "Urban", "MT"]
  sr <- d1[d1$Region == r & d1$Group == "Rural", "MT"]
  se <- sqrt(MSE * (1/length(su) + 1/length(sr)))
  tval <- (mean(su) - mean(sr)) / se
  list(diff = mean(su) - mean(sr), se = se, t = tval, df = dfe,
       p = 2 * pt(-abs(tval), dfe), lsd = qt(0.975, dfe) * se)
}
L <- lapply(c("HP", "NCL"), lsd_line)
se_tab <- data.frame(
  脑区 = c("HP", "NCL"),
  Urban = sapply(c("HP", "NCL"), function(r) {
    v <- d1[d1$Region == r & d1$Group == "Urban", "MT"]; sprintf("%.3f ± %.3f", mean(v), sem(v)) }),
  Rural = sapply(c("HP", "NCL"), function(r) {
    v <- d1[d1$Region == r & d1$Group == "Rural", "MT"]; sprintf("%.3f ± %.3f", mean(v), sem(v)) }),
  差值 = round(sapply(L, `[[`, "diff"), 3),
  LSD = round(sapply(L, `[[`, "lsd"), 3),
  t = round(sapply(L, `[[`, "t"), 3),
  df = sapply(L, `[[`, "df"),
  p = format.pval(sapply(L, `[[`, "p"), digits = 3, eps = 0.001),
  显著性 = sig(sapply(L, `[[`, "p")), check.names = FALSE)
cat("\n## 同一脑区城乡比较三线表 (Fisher LSD) ##\n"); print(se_tab, row.names = FALSE)

wr(anova_tab, "data/AA_HPNCL_双因素方差分析_三线表.csv")
wr(se_tab,    "data/AA_HPNCL_城乡比较_三线表.csv")

## ================= AA_PVN.csv：独立样本 t 检验 =================
d2 <- read.csv("data/AA_PVN.csv", stringsAsFactors = FALSE)
tt <- t.test(d2$Urban, d2$Rural)   # Welch (默认，不假设方差齐)

desc <- data.frame(
  组 = c("Urban", "Rural"),
  n = c(length(d2$Urban), length(d2$Rural)),
  `Mean±SEM` = c(sprintf("%.3f ± %.3f", mean(d2$Urban), sem(d2$Urban)),
                 sprintf("%.3f ± %.3f", mean(d2$Rural), sem(d2$Rural))),
  check.names = FALSE)

tt_tab <- data.frame(
  比较 = "Urban vs Rural",
  t = round(as.numeric(tt$statistic), 3),
  df = round(as.numeric(tt$parameter), 3),
  p = format.pval(tt$p.value, digits = 3, eps = 0.001),
  显著性 = sig(tt$p.value), check.names = FALSE)

cat("\n========== AA_PVN：独立样本 t 检验 ==========\n")
cat("## 描述统计 (Mean±SEM) ##\n"); print(desc, row.names = FALSE)
cat("\n## t 检验三线表 ##\n"); print(tt_tab, row.names = FALSE)

wr(desc,   "data/AA_PVN_描述统计_三线表.csv")
wr(tt_tab, "data/AA_PVN_t检验_三线表.csv")
