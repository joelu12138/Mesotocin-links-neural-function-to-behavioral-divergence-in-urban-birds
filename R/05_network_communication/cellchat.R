# Cell-cell communication inference with CellChat v2.2.0 (urban vs rural).
# Database: CellChatDB.treeSparrow_conservative (1-to-1 ortholog map, 1757
# ligand-receptor pairs; matches Methods), built from CellChatDB.human via
# OrthoFinder ortholog mapping by build_sparrow_cellchat_db.R (which also
# produces a full ortholog variant, CellChatDB.treeSparrow, 12763 pairs).
# truncatedMean probability estimation; groups with <10 cells excluded;
# 100 permutations; interactions retained at P < 0.05; objects merged to
# compare interaction numbers, aggregate strengths and differential networks.
#
# NOTE: the DB and the expression matrix (subobj) share the same gene-symbol
# namespace. The ortholog TSV inputs of the build script are upstream artifacts
# not included here; the prebuilt .rda is shipped ready to use.
#
# Input: neuronal subtype object saved by 04_snrna_seq/neuronal_subtype_diff.R.

library(CellChat)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

data_dir <- "data"   # run from repository root; inputs/outputs live in data/
setwd(data_dir)

load("neuronal_subtypes.Rdata")   # produces subobj ($celltype, $group Urban/Rural)

# Species-adapted ligand-receptor database (conservative 1-to-1, 1757 pairs).
load("../R/05_network_communication/CellChatDB.treeSparrow_conservative.rda")  # -> CellChatDB.treeSparrow_conservative
CellChatDB <- CellChatDB.treeSparrow_conservative
showDatabaseCategory(CellChatDB)
write.csv(CellChatDB$interaction, "database.csv", row.names = FALSE)

# ---- helper: build a CellChat object from one habitat subset ----
run_cellchat <- function(obj, db) {
  DefaultAssay(obj) <- "RNA"
  data.input <- GetAssayData(obj, assay = "RNA", slot = "data")
  data.input <- data.input[!grepl("^LOC", rownames(data.input)), ]  # drop uncharacterized LOC genes
  labels <- obj@meta.data$celltype                                 # labels from THIS subset (aligned)
  meta <- data.frame(group = labels, row.names = colnames(data.input))
  meta$group <- droplevels(factor(meta$group))
  stopifnot(all(colnames(data.input) == rownames(meta)))

  cc <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # truncatedMean + 100 permutations (Methods)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                          raw.use = TRUE, population.size = TRUE, nboot = 100)
  cc <- filterCommunication(cc, min.cells = 10)   # exclude groups with <10 cells
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# ---- helper: standard per-habitat visualizations + tables ----
plot_cellchat <- function(cc, prefix) {
  write.csv(cc@LR$LRsig, paste0("cellchat-LR-", prefix, ".csv"), row.names = FALSE)
  df.net <- subsetCommunication(cc)               # default thresh = 0.05 (P < 0.05)
  write.csv(df.net, paste0(prefix, "-cell_communications_all.csv"), row.names = FALSE)

  groupSize <- as.numeric(table(cc@idents))
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cc@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cc@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction weights/strength")

  gg1 <- netVisual_heatmap(cc)
  gg2 <- netVisual_heatmap(cc, measure = "weight")
  pdf(paste0("combined_heatmap-", prefix, ".pdf"), width = 12, height = 6)
  draw(gg1 + gg2)
  dev.off()

  netVisual_bubble(cc, remove.isolate = FALSE)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  for (pattern in c("outgoing", "incoming", "all")) {
    h <- netAnalysis_signalingRole_heatmap(cc, pattern = pattern, width = 10, height = 20)
    pdf(paste0(pattern, "_", prefix, ".pdf"), width = 6, height = 12)
    draw(h)
    dev.off()
  }
  cc
}

# ---- Rural and Urban analyzed separately ----
dt_rural <- subset(subobj, group == "Rural")
Idents(dt_rural) <- dt_rural@meta.data$celltype
cellchat.rural <- run_cellchat(dt_rural, CellChatDB)
cellchat.rural <- plot_cellchat(cellchat.rural, "Rural")

dt_urban <- subset(subobj, group == "Urban")
Idents(dt_urban) <- dt_urban@meta.data$celltype
cellchat.urban <- run_cellchat(dt_urban, CellChatDB)
cellchat.urban <- plot_cellchat(cellchat.urban, "Urban")

# ---- merge to compare urban vs rural ----
object.list <- list(Rural = cellchat.rural, Urban = cellchat.urban)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
pdf("combined_heatmap-diff.pdf", width = 6, height = 8)
print(gg1 + gg2)
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("count", "weight"))
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = weight.max[1], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE, measure = "weight")
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE, measure = "weight",
                          sources.use = "Neuron")

gg1 <- netVisual_heatmap(cellchat.merged)
gg2 <- netVisual_heatmap(cellchat.merged, measure = "weight")
pdf("Ser_heatmap.pdf", height = 6, width = 12)
draw(gg1 + gg2)
dev.off()

p5 <- rankNet(cellchat.merged, mode = "comparison", stacked = TRUE,
              do.stat = TRUE, paired.test = FALSE)
p5 <- p5 + scale_fill_manual(values = c("Rural" = "#4492A3", "Urban" = "#DF5E21"))
pdf("cellchat-pathway-diff.pdf", height = 8, width = 6)
print(p5)
dev.off()

p5 <- rankNet(cellchat.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
write.csv(p5$data, "RA_rank.csv", row.names = FALSE)

saveRDS(cellchat.merged, "cellchat_merged.rds")