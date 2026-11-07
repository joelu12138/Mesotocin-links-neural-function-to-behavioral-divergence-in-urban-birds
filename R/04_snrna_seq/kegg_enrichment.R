# KEGG enrichment of differential genes (snRNA-seq)
# Methods: Q < 0.1 used as the enrichment threshold.
# TERM2GENE / TERM2NAME come from a user-curated pathway-gene table
# (sparrow genes mapped to KEGG orthologs), not from AnnotationHub.

library(clusterProfiler)

# ---- input paths (edit before running) ----
pathway_table_file <- "inputs/PATH_ID_NAME.csv"  # cols: pathway_id, gene, pathway_name
deg_file           <- "inputs/DEGs.csv"          # must contain column 'Index' (gene symbols)

PATH_ID_NAME <- read.csv(pathway_table_file, header = TRUE)
head(PATH_ID_NAME)

go2gene <- PATH_ID_NAME[, c(1, 2)]
go2name <- PATH_ID_NAME[, c(1, 3)]
data2 <- read.csv(deg_file, header = TRUE)
head(data2)

x <- enricher(
  data2$Index,
  TERM2GENE = go2gene,
  TERM2NAME = go2name,
  pvalueCutoff = 1,
  qvalueCutoff = 0.1)
x

dotplot(x)
a <- x@result
write.csv(a, "KEGG.csv", row.names = FALSE)