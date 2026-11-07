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
  gene_select = "fraction", # 默认fraction;其他有variable(Seurat对象中存储的基因),custom(自定义)
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Neuron" # the name of the hdWGCNA experiment
)
# 各组构建metacell
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "orig.ident"), #指定seurat_obj@meta.data中要分组的列
  reduction = 'harmony', # 选择要执行KNN的降维
  k = 25, # 最近邻居参数
  max_shared = 10, # 两个metacell之间共享细胞的最大数目
  ident.group = 'celltype', # 设置metacell安全对象的标识
  min_cells = 100 # 排除数量小于100的细胞亚群
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
# check以下数据
power_table <- GetPowerTable(seurat_obj)
head(power_table)

# 如果没有指定软阈值，construcNetwork会自动指定软阈值
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
  overwrite_tom = TRUE# 允许覆盖已存在的同名文件
)

# 可视化WGCNA树状图
# “灰色”模块由那些未被归入任何共表达模块的基因组成。对于所有下游分析和解释，应忽略灰色模块。
PlotDendrogram(seurat_obj, main='Neuron hdWGCNA Dendrogram')
# 需要先运行ScaleData，否则harmony会报错:
# seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# 计算完整单细胞数据集中的所有MEs
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)

# 协调模特征基因:
# 允许用户对MEs应用Harmony批量校正，生成协调模块特征基因(hMEs)
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
#MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# 计算基于特征基因的连接性(kME)：
# 关注枢纽基因
seurat_obj <- ModuleConnectivity(
  seurat_obj
)

# 模块重命名
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Neuron_NEW"
)
# 绘制每个模块按kME排序的基因
p <- PlotKMEs(seurat_obj, ncol=4)
p

# 获取模块内部信息:
# 个人认为这一部分很关键，毕竟咱们就是想要得到不同模块的基因
# 这一步去除了不需要的灰色模块基因
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey')
# 显示前六列:
head(modules[,1:6])
write.csv(modules,"modules.csv")

# 得到枢纽基因
# 可以提取按kME排序的前N个枢纽基因的表格,这里选择了10
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)
write.csv(hub_df,"hub_df.csv")
# 保存数据
saveRDS(seurat_obj, 'hdWGCNA_object.rds')

# 计算每个模块前25个枢纽基因的kME得分
# 使用UCell方法
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell' # Seurat方法(AddModuleScore)
)

# 每个模块制作hMEs的特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# 每个模块制作hub scores的特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

# 每个模块在不同样本中的情况（按 Urban/Rural 分组，使用全部神经元）
ModuleRadarPlot(
  seurat_obj,
  group.by = 'group',
  barcodes = rownames(seurat_obj@meta.data),
  axis.label.size=4,
  grid.label.size=4
)

# 查看模块相关图
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


# 使用ModuleNetworkPlot可视化每个模块前50(数值可自定)的hub gene
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

# hubgene network(基因数可自定)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other=2,
  edge_prop = 0.75,
  mods = 'all'
)

# 可以选择模块数
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

#基因集的富集分析
###基因集的富集分析
#基因集的打分
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
#整合差异基因集
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

#获取node，需要基因及module两列即可
gene_module <- hub_df[,1:2]
gene_module <- gene_module[gene_module$gene_name %in% unique(c(unique(edge_file$from),unique(edge_file$to))),]
#module color（灰色模块不参与下游分析，剔除）
mod_color_df <- GetModules(seurat_obj) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>%
  filter(module != 'grey') %>%
  arrange(module)
gene_module <- gene_module %>%
  left_join(mod_color_df, by = "module")
write.csv(gene_module, file = 'gene_module.csv', row.names = FALSE)

