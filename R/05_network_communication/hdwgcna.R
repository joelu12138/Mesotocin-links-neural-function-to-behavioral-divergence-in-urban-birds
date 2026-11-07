# hdWGCNA v0.4.12 co-expression network on neuronal nuclei.
# Methods: genes in >=5% of neuronal nuclei; unsigned network, soft power 4;
# minModuleSize 50, deepSplit 4, mergeCutHeight 0.2 -> 11 modules (M11 focal).
# Input: neuronal subtype object saved by 04_snrna_seq/neuronal_subtype_diff.R.

rm(list = ls())
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
library(reshape2)

load("neuronal_subtypes.Rdata")   # produces subobj (neuronal subtypes)

set.seed(12345)
table(subobj@meta.data$orig.ident)
seurat_obj <- SetupForWGCNA(
  subobj,
  gene_select = "fraction", # default fraction; other options are variable (genes stored in the Seurat object) and custom (user-defined)
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Neuron" # the name of the hdWGCNA experiment
)
# Build metacells for each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "orig.ident"), # specify the columns in seurat_obj@meta.data used for grouping
  reduction = 'harmony', # the dimensional reduction used for KNN
  k = 25, # nearest neighbor parameter
  max_shared = 10, # maximum number of cells shared between two metacells
  ident.group = 'celltype', # set the identity of the metacell Seurat object
  min_cells = 100 # exclude cell subgroups with fewer than 100 cells
)
# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = 'unsigned'# you can also use "unsigned" or "signed hybrid"
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
# check the following data
power_table <- GetPowerTable(seurat_obj)
head(power_table)

# If no soft power is specified, ConstructNetwork selects one automatically
# construct co-expression network (Methods: soft power 4, minModuleSize 50,
# deepSplit 4, mergeCutHeight 0.2):
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 4,
  networkType = 'unsigned',
  TOMType = 'unsigned',
  minModuleSize = 50,
  deepSplit = 4,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM",
  tom_name = 'Neuron', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE# allow overwriting an existing file with the same name
)

# Visualize the WGCNA dendrogram
# The "grey" module consists of genes that were not assigned to any co-expression module. The grey module should be ignored in all downstream analyses and interpretations.
PlotDendrogram(seurat_obj, main='Neuron hdWGCNA Dendrogram')
# ScaleData must be run first, otherwise harmony reports an error:
# seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# Compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)

# Harmonized module eigengenes:
# Apply Harmony batch correction to the MEs to generate harmonized module eigengenes (hMEs)
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
#MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# Compute eigengene-based connectivity (kME):
# focus on hub genes
seurat_obj <- ModuleConnectivity(
  seurat_obj
)

# Rename modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Neuron_NEW"
)
# Plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=4)
p

# Get within-module information:
# This part is key, since the goal is to obtain the genes of each module
# This step removes the unwanted grey module genes
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey')
# Show the first six columns:
head(modules[,1:6])
write.csv(modules,"modules.csv")

# Obtain hub genes
# Extract a table of the top N hub genes ranked by kME; 10 is used here
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)
write.csv(hub_df,"hub_df.csv")
# Save the data
saveRDS(seurat_obj, 'hdWGCNA_object.rds')

# Compute kME scores of the top 25 hub genes in each module
# Use the UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell' # Seurat method (AddModuleScore)
)

# Feature plots of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# Feature plots of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# Module behavior in different samples (by group Urban/Rural, using all neurons)
ModuleRadarPlot(
  seurat_obj,
  group.by = 'group',
  barcodes = rownames(seurat_obj@meta.data),
  axis.label.size=4,
  grid.label.size=4
)

# View the module correlogram
ModuleCorrelogram(seurat_obj)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'celltype')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p


# Use ModuleNetworkPlot to visualize the top 50 (adjustable) hub genes of each module
ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)


options(future.globals.maxSize = 5 * 1024^3)  # 5GB

# hubgene network (number of genes is adjustable)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other=2,
  edge_prop = 0.75,
  mods = 'all'
)

# select the number of modules
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other= 2,
  edge_prop = 0.75,
  mods = mods[1:5] # only select 5 modules
)


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1# min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2# size of each point based on intramodular connectivity
  ) +
  umap_theme()


ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

#Gene set enrichment analysis
###Gene set enrichment analysis
#Gene set scoring
cran.packages<-c("msigdbr","dplyr","purrr","stringr","magrittr",
                 "RobustRankAggreg","tibble","reshape2","ggsci",
                 "tidyr","aplot","ggfun","ggplotify","ggridges",
                 "gghalves","Seurat","SeuratObject","methods",
                 "devtools","BiocManager","data.table","doParallel",
                 "doRNG")
if(!requireNamespace(cran.packages,quietly=TRUE)){
  install.packages(cran.packages,ask=F,update=F)
}
#??install??packages??from??Bioconductor
bioconductor.packages<-c("GSEABase","AUCell","SummarizedExperiment",
                         "singscore","GSVA","ComplexHeatmap","ggtree",
                         "Nebulosa")
if(!requireNamespace(bioconductor.packages,quietly=TRUE)){
  BiocManager::install(bioconductor.packages,ask=F,update=F)
}
#??install??packages??from??Github
if(!requireNamespace("UCell",quietly=TRUE)){
  devtools::install_github("carmonalab/UCell")
}
if(!requireNamespace("irGSEA",quietly=TRUE)){
  devtools::install_github("chuiqin/irGSEA")
}
ST <- subobj   # AUCell module/gene-set scores across all neuronal subtypes
table(ST@meta.data$celltype)
library(irGSEA)
df<-read.csv("pathway.csv",header = T)
table(df$names)
df$names <- as.character(df$names)
df$SYMBOL <- as.character(df$SYMBOL)
gene_list <- split(df$SYMBOL, df$names)
names(ST@assays)
expr2 <- irGSEA.score(object = ST,assay = "RNA", 
                      slot = "data", seeds = 123, ncores = 1,                    
                      min.cells = 3, min.feature = 0,                    
                      custom = T, geneset = gene_list,                     
                      species = "Homo sapiens",                    
                      method = c("AUCell"),                    
                      aucell.MaxRank = NULL, ucell.MaxRank = NULL,                     
                      kcdf = 'Gaussian')
DefaultAssay(expr2) <- "AUCell"
#Integrate the differential gene sets
head(expr2@meta.data)
result.dge<-irGSEA.integrate(object=expr2,
                             group.by="celltype",
                             metadata=NULL,col.name=NULL,
                             method=c("AUCell"))
write.csv(result.dge$AUCell,"KEGG_celltype.csv")

#######
TOM <- GetTOM(seurat_obj, wgcna_name='Neuron')
hub_df <- GetHubGenes(seurat_obj, n_hubs = 300)
TOM <- TOM[hub_df$gene_name,hub_df$gene_name]
dim(TOM)

adj_mat <- melt(TOM, varnames = c("from", "to"), value.name = "tom")
adj_mat <- adj_mat[adj_mat$tom>0.02,]
edge_file <- adj_mat[,1:2]
write.csv(edge_file, file = 'edge_file.csv', row.names = FALSE)

#Get nodes; the gene and module columns are sufficient
gene_module <- hub_df[,1:2]
gene_module <- gene_module[gene_module$gene_name %in% unique(c(unique(edge_file$from),unique(edge_file$to))),]
#module color (the grey module is excluded from downstream analyses and removed)
mod_color_df <- GetModules(seurat_obj) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>%
  filter(module != 'grey') %>%
  arrange(module)
gene_module <- gene_module %>%
  left_join(mod_color_df, by = "module")
write.csv(gene_module, file = 'gene_module.csv', row.names = FALSE)

