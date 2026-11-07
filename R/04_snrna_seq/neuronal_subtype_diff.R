# Neuronal subtype analysis: extract neuronal nuclei, recluster (resolution 0.8),
# annotate 10 excitatory + 3 inhibitory subtypes, composition, subtype DE and
# gene-set scoring.
#
# Methods: the 7,555 neuronal nuclei were extracted and reclustered at resolution
# 0.8; excitatory/inhibitory identity assigned with SLC17A6/NRGN/CAMK2A and
# GAD1/GAD2/SLC32A1. Input is the integrated object saved by single_cell_diff.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(irGSEA)
  library(UCell)
})
options(stringsAsFactors = FALSE)

data_dir <- "data"   # adjust to local data directory
setwd(data_dir)

# ---- 1. Extract neurons and recluster --------------------------------------
load("sce_all_integrated.Rdata")          # produces sce.all (from single_cell_diff.R)
subobj <- subset(sce.all, celltype == "Neuron")

f <- "neuronal_subtypes.Rdata"
if (!file.exists(f)) {
  subobj <- subobj %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(pc.genes = VariableFeatures(.)) %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:30, reduction = "harmony") %>%
    FindClusters(resolution = 0.8) %>%
    RunUMAP(dims = 1:30, reduction = "harmony")
  save(subobj, file = f)
} else {
  load(f)
}

# ---- 2. Annotate 13 neuronal subtypes (manual curation) ---------------------
Idents(subobj) <- "seurat_clusters"
cl <- as.character(subobj$seurat_clusters)
ct <- rep(NA_character_, length(cl))
ct[cl %in% c("17")]                  <- "In1"
ct[cl %in% c("8", "19", "12")]       <- "In2"
ct[cl %in% c("6", "18", "14")]       <- "In3"
ct[cl %in% c("1", "27")]             <- "Ex1"
ct[cl %in% c("28", "23", "3", "30")] <- "Ex2"
ct[cl %in% c("22")]                  <- "Ex3"
ct[cl %in% c("21", "10", "4", "26")] <- "Ex4"
ct[cl %in% c("24", "29", "0")]       <- "Ex5"
ct[cl %in% c("7")]                   <- "Ex6"
ct[cl %in% c("15", "5", "20", "25")] <- "Ex7"
ct[cl %in% c("2", "13")]             <- "Ex8"
ct[cl %in% c("9")]                   <- "Ex9"
ct[cl %in% c("16", "11")]            <- "Ex10"
subobj$celltype <- ct
Idents(subobj) <- "celltype"
table(subobj$celltype)

# ---- 3. Subtype markers -----------------------------------------------------
markers <- FindAllMarkers(subobj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.59)
write.csv(markers, "Markers_neuronal_subtype.csv", row.names = FALSE)

# ---- 4. Urban/rural grouping and subtype composition -----------------------
subobj$group <- ifelse(subobj$orig.ident %in% c("U1", "U2", "U3"), "Urban", "Rural")
cell_counts <- table(Idents(subobj))
cell.num.group <- table(Idents(subobj), subobj$group)
cell.freq.group <- round(prop.table(cell.num.group, margin = 2) * 100, 2)
cell.all <- cbind(cell_counts = cell_counts,
                  cell_Freq = round(prop.table(cell_counts) * 100, 2),
                  cell.num.group, cell.freq.group)
write.csv(cell.all, "neuronal_subtype_composition.csv")

# ---- 5. Within-subtype urban-rural DE (example: Ex3; Wilcoxon, BH) ----------
sub.obj <- subset(subobj, idents = "Ex3")
dfmarkers <- FindMarkers(sub.obj, ident.1 = "Urban", ident.2 = "Rural",
                         group.by = "group", min.pct = 0.25,
                         logfc.threshold = 0.59, verbose = FALSE)
write.csv(dfmarkers, "FindMarkers-Ex3.csv")

# ---- 6. Gene-set scoring (irGSEA / AUCell) ----------------------------------
df <- read.csv("marker.csv", header = TRUE)
gene_list <- split(as.character(df$SYMBOL), as.character(df$names))
expr2 <- irGSEA.score(object = subobj, assay = "RNA", slot = "data",
                      seeds = 123, ncores = 1, min.cells = 3, min.feature = 0,
                      custom = TRUE, geneset = gene_list, species = "Homo sapiens",
                      method = "AUCell", aucell.MaxRank = NULL,
                      ucell.MaxRank = NULL, kcdf = "Gaussian")
DefaultAssay(expr2) <- "AUCell"
result.dge <- irGSEA.integrate(object = expr2, group.by = "celltype",
                               method = "AUCell")

# Urban vs rural gene-set comparison (per habitat)
result.dge.group <- irGSEA.integrate(object = expr2, group.by = "group",
                                     method = "AUCell")
write.csv(result.dge.group$AUCell, "geneset_AUCell_group_subtype.csv", row.names = FALSE)