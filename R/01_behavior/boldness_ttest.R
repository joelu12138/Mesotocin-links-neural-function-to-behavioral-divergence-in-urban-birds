files <- c("Alerting distance", "Flight initiation distance")

rows <- lapply(files, function(f) {
  d <- read.csv(sprintf("data/%s.csv", f), stringsAsFactors = FALSE)
  u <- na.omit(as.numeric(d$Urban)); r <- na.omit(as.numeric(d$Rural))
  w <- t.test(u, r, var.equal = TRUE)
  sig <- if (w$p.value < 0.001) "***" else if (w$p.value < 0.01) "**" else if (w$p.value < 0.05) "*" else "ns"
  p   <- if (w$p.value < 0.001) "<0.001" else sprintf("%.3f", w$p.value)
  data.frame(Metric = f,
             Urban = sprintf("%.3f+/-%.3f (n=%d)", mean(u), sd(u)/sqrt(length(u)), length(u)),
             Rural = sprintf("%.3f+/-%.3f (n=%d)", mean(r), sd(r)/sqrt(length(r)), length(r)),
             t = round(w$statistic, 3), df = round(w$parameter, 2), p = p, Significance = sig,
             check.names = FALSE)
})

out <- do.call(rbind, rows)

path <- "data/three_line_table_stats.csv"
con <- file(path, "w", encoding = "UTF-8")
writeLines("\ufeff", con, sep = "")   # UTF-8 BOM so Excel opens the CSV cleanly
close(con)
write.table(out, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, append = TRUE, fileEncoding = "UTF-8")

cat("Wrote:", path, "\n\n")
print(out, row.names = FALSE)
