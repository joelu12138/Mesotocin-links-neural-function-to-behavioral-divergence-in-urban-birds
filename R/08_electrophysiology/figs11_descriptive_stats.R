files <- list.files("data", pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("descriptive_statistics", files)]
sem <- function(x) sd(x)/sqrt(length(x))

out <- do.call(rbind, lapply(files, function(f) {
  base <- sub("\\.csv$", "", basename(f))
  region <- sub("^.*-", "", base)          # text after the last hyphen = HP / NCL
  metric <- sub("-(HP|NCL)$", "", base)
  d <- read.csv(f)
  data.frame(
    Variable = metric,
    `Brain region` = region,
    Urban = sprintf("%.3f +/- %.3f", mean(d$Urban), sem(d$Urban)),
    Rural = sprintf("%.3f +/- %.3f", mean(d$Rural), sem(d$Rural)),
    check.names = FALSE)
}))

out <- out[order(out$Variable, out$`Brain region`), ]
out$Variable[duplicated(out$Variable)] <- ""   # the variable name is shown only once, on the HP row
print(out, row.names = FALSE)

con <- file("data/descriptive_statistics.csv", "w", encoding = "UTF-8")
writeLines("\ufeff", con, sep = ""); close(con)
write.table(out, file = "data/descriptive_statistics.csv", sep = ",",
            row.names = FALSE, col.names = TRUE, quote = TRUE, append = TRUE,
            fileEncoding = "UTF-8")
cat("\nWrote: data/descriptive_statistics.csv\n")
