files <- list.files("data", pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("descriptive_statistics", files)]
sem <- function(x) sd(x)/sqrt(length(x))

display <- c("Event frequency" = "Event frequency (Hz)",
             "Membrane potential" = "Membrane potential (mV)",
             "Spikes" = "Spikes")

out <- do.call(rbind, lapply(files, function(f) {
  base <- sub("\\.csv$", "", basename(f))
  region <- sub("^.*-", "", base)          # text after the last hyphen = HP / NCL
  metric <- sub("-(HP|NCL)$", "", base)
  d <- read.csv(f)
  tt <- t.test(d$Urban, d$Rural, var.equal = TRUE)   # independent-samples t-test (Student)
  data.frame(
    Variable = unname(display[metric]),
    `Brain region` = region,
    Urban = sprintf("%.3f +/- %.3f", mean(d$Urban), sem(d$Urban)),
    Rural = sprintf("%.3f +/- %.3f", mean(d$Rural), sem(d$Rural)),
    t = round(as.numeric(tt$statistic), 3),
    df = as.integer(tt$parameter),
    P = format.pval(tt$p.value, digits = 3, eps = 0.001),
    check.names = FALSE)
}))

out <- out[order(out$Variable, out$`Brain region`), ]
out$Variable[duplicated(out$Variable)] <- ""   # the variable name is shown only once, on the HP row
print(out, row.names = FALSE)

con <- file("data/t_test_results.csv", "w", encoding = "UTF-8")
writeLines("\ufeff", con, sep = ""); close(con)
write.table(out, file = "data/t_test_results.csv", sep = ",",
            row.names = FALSE, col.names = TRUE, quote = TRUE, append = TRUE,
            fileEncoding = "UTF-8")
cat("\nWrote: data/t_test_results.csv\n")
