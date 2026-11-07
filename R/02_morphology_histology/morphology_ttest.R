library(dplyr)

# 城乡两组的形态学特征文件（每个文件三列：Index, Group, <特征>）
files <- c("Tarsus.csv", "Wing.csv", "Beak.csv", "Brain.csv",
           "HWI.csv", "Ratio.csv", "mass.csv", "HP.csv")
files <- files[file.exists(files)]   # 只跑目录里实际存在的文件

equal_var <- FALSE   # FALSE = Welch t 检验（不假设方差齐，稳健）；TRUE = 经典 Student t 检验

res <- lapply(files, function(f) {
  d <- read.csv(f)
  feature <- setdiff(names(d), c("Index", "Group"))
  d$Value <- d[[feature]]
  d$Group <- factor(d$Group, levels = c("Urban", "Rural"))

  tt <- t.test(Value ~ Group, data = d, var.equal = equal_var)

  data.frame(
    Feature = feature,
    N_Urban = sum(d$Group == "Urban"),
    N_Rural = sum(d$Group == "Rural"),
    Mean_Urban = mean(d$Value[d$Group == "Urban"]),
    Mean_Rural = mean(d$Value[d$Group == "Rural"]),
    SD_Urban   = sd(d$Value[d$Group == "Urban"]),
    SD_Rural   = sd(d$Value[d$Group == "Rural"]),
    t  = as.numeric(tt$statistic),
    df = as.numeric(tt$parameter),
    p  = tt$p.value
  )
}) %>% bind_rows()

# 多重检验校正：Bonferroni
res$p_Bonferroni <- p.adjust(res$p, method = "bonferroni")
res$sig_Bonferroni <- cut(res$p_Bonferroni, breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c("***", "**", "*", "ns"), right = FALSE)

out <- res %>%
  mutate(across(c(Mean_Urban, Mean_Rural, SD_Urban, SD_Rural), ~ round(.x, 3)),
         t  = round(t, 3),
         df = round(df, 1),
         across(c(p, p_Bonferroni), ~ round(.x, 3)))

print(as.data.frame(out))
write.csv(out, "morphology_ttest_results.csv", row.names = FALSE)
cat("\nsaved morphology_ttest_results.csv\n")
