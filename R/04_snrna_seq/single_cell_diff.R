# Single-nucleus RNA-seq: sample integration, clustering, annotation, urban-rural DE
# and gene-set scoring.
#
# Methods: Seurat v5 + Harmony (sample identity as integration variable), 30 dims,
# Louvain resolution 0.5. QC: 200-8000 detected genes, >=500 UMIs, <5% ribosomal-
# protein and <5% hemoglobin transcripts. Doublets were removed per sample with
# DoubletFinder v2.0.2 (pK = 0.09) upstream of this script (see Methods).
#
# Input/output paths point to the local data directory; adjust before running.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(irGSEA)
  library(UCell)
})
options(stringsAsFactors = FALSE)

data_dir <- "data"   # input/output directory; adjust to local data
setwd(data_dir)

# ---- 1. Read 10X matrices and merge samples ---------------------------------
multi_sample_path <- "data/inputs"   # folder of per-sample 10X outputs
sample_folders <- list.dirs(multi_sample_path, full.names = TRUE, recursive = FALSE)
seurat_list <- lapply(sample_folders, function(folder) {
  m <- Read10X(data.dir = folder)
  s <- basename(folder)
  colnames(m) <- paste0(s, "_", colnames(m))
  CreateSeuratObject(counts = m, project = s)
})
sce.all <- JoinLayers(merge(seurat_list[[1]], seurat_list[-1]))
table(sce.all$orig.ident)

# ---- 2. Quality control (per Methods) --------------------------------------
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")
for (p in c("percent.mt", "percent.rp", "percent.hb")) {
  sce.all@meta.data[[p]][is.na(sce.all@meta.data[[p]])] <- 0
}
sce.all <- subset(sce.all,
                  subset = nFeature_RNA >= 200 & nFeature_RNA <= 8000 &
                           nCount_RNA >= 500 & percent.rp < 5 & percent.hb < 5)
table(sce.all$orig.ident)

# ---- 3. Normalize, Harmony integration, cluster, embed ---------------------
sce.all <- sce.all %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(pc.genes = VariableFeatures(.)) %>%
  RunHarmony("orig.ident") %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:30, reduction = "harmony")

# ---- 4. Manual annotation of major cell classes (markers in Table 13) -------
# 11 major classes; astrocytes resolved into 4 subtypes. Cluster-to-class mapping
# reflects the final curated assignment (canonical markers + DotPlot/VlnPlot).
Idents(sce.all) <- "seurat_clusters"
cl <- as.character(sce.all$seurat_clusters)
ct <- rep(NA_character_, length(cl))
ct[cl %in% c("17", "18", "19", "21")]        <- "Neuron"
ct[cl %in% c("3", "6")]                       <- "Astrocytes-1"
ct[cl %in% c("2")]                            <- "Astrocytes-2"
ct[cl %in% c("7")]                            <- "Astrocytes-3"
ct[cl %in% c("5", "22")]                      <- "Astrocytes-4"
ct[cl %in% c("0", "1", "10", "14", "27", "30")] <- "Oligodendrocytes"
ct[cl %in% c("8", "15", "25", "28")]          <- "Microglial cells"
ct[cl %in% c("24", "9", "31")]                <- "Endothelial cells"
ct[cl %in% c("29")]                           <- "OPCs"
ct[cl %in% c("26", "11")]                     <- "Ependymal cells"
ct[cl %in% c("13")]                           <- "Red blood cells"
ct[cl %in% c("20")]                           <- "T cells"
ct[cl %in% c("16")]                           <- "CPEs"
ct[cl %in% c("12", "32")]                     <- "Heterophils"
sce.all$celltype <- ct
Idents(sce.all) <- "celltype"
table(sce.all$celltype)

# ---- 5. Cluster markers -----------------------------------------------------
markers <- FindAllMarkers(sce.all, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.59)
write.csv(markers, "Markers_cell_type.csv", row.names = FALSE)

# ---- 6. Urban/rural grouping and cell composition --------------------------
sce.all$group <- ifelse(sce.all$orig.ident %in% c("U1", "U2", "U3"), "Urban", "Rural")
cell_counts <- table(Idents(sce.all))
cell.num.group <- table(Idents(sce.all), sce.all$group)
cell.freq.group <- round(prop.table(cell.num.group, margin = 2) * 100, 2)
cell.all <- cbind(cell_counts = cell_counts,
                  cell_Freq = round(prop.table(cell_counts) * 100, 2),
                  cell.num.group, cell.freq.group)
write.csv(cell.all, "cell_composition.csv")

# ---- 7. Neuronal urban-rural differential expression (Wilcoxon, BH) ---------
subobj <- subset(sce.all, idents = "Neuron")
dfmarkers <- FindMarkers(subobj, ident.1 = "Urban", ident.2 = "Rural",
                         group.by = "group", min.pct = 0.25,
                         logfc.threshold = 0.59, verbose = FALSE)
write.csv(dfmarkers, "Group-Neuron.csv")

# ---- 8. Gene-set scoring (irGSEA / AUCell, human-ortholog gene sets) ---------
df <- read.csv("marker.csv", header = TRUE)
gene_list <- split(as.character(df$SYMBOL), as.character(df$names))
expr2 <- irGSEA.score(object = sce.all, assay = "RNA", slot = "data",
                      seeds = 123, ncores = 1, min.cells = 3, min.feature = 0,
                      custom = TRUE, geneset = gene_list, species = "Homo sapiens",
                      method = "AUCell", aucell.MaxRank = NULL,
                      ucell.MaxRank = NULL, kcdf = "Gaussian")
DefaultAssay(expr2) <- "AUCell"
result.dge <- irGSEA.integrate(object = expr2, group.by = "celltype",
                               method = "AUCell")
write.csv(result.dge$AUCell, "geneset_AUCell_celltype.csv", row.names = FALSE)

# Urban vs rural gene-set comparison (per habitat)
result.dge.group <- irGSEA.integrate(object = expr2, group.by = "group",
                                     method = "AUCell")
write.csv(result.dge.group$AUCell, "geneset_AUCell_group.csv", row.names = FALSE)

# ---- 9. Save integrated object for downstream scripts ----------------------
save(sce.all, file = "sce_all_integrated.Rdata")