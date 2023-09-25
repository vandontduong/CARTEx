# extract newborn data
# Vandon Duong

experiment = 'GSE136184'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(gtools)



####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

CART.combined.CD8pos <- readRDS(paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "AgeGroup")

head(CART.combined.CD8pos)

# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical <- VlnPlot(CART.combined.CD8pos, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "AgeGroup2", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical, paste('./plots/', experiment, '_vlnplot_naiveT_canonical', sep = ''))

####################################################################################################
######################################## Extract newborn data ######################################
####################################################################################################


CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "AgeGroup2")
CART.combined.CD8pos.newborn <- subset(CART.combined.CD8pos, idents = c("Newborn"))


# Quality filter
CART.combined.CD8pos.newborn <- subset(CART.combined.CD8pos.newborn, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.MT < 10)

# Examine features
vlnplot_quality <- VlnPlot(object = CART.combined.CD8pos.newborn, features = c('nFeature_RNA','nCount_RNA', "percent.MT"), group.by = 'AgeGroup2', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality_newborn', sep = ''))

# Variable features and initial UMAP analysis
CART.combined.CD8pos.newborn <- FindVariableFeatures(CART.combined.CD8pos.newborn, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(CART.combined.CD8pos.newborn), 10)
varplt <- VariableFeaturePlot(CART.combined.CD8pos.newborn)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled_newborn', sep = ''))

CART.combined.CD8pos.newborn <- ScaleData(CART.combined.CD8pos.newborn, features = all.genes)
CART.combined.CD8pos.newborn <- RunPCA(CART.combined.CD8pos.newborn, features = VariableFeatures(object = CART.combined.CD8pos.newborn), npcs = 40)

CART.combined.CD8pos.newborn <- FindNeighbors(CART.combined.CD8pos.newborn, dims = 1:10)
CART.combined.CD8pos.newborn <- FindClusters(CART.combined.CD8pos.newborn, resolution = 0.5)
CART.combined.CD8pos.newborn <- RunUMAP(CART.combined.CD8pos.newborn, dims = 1:10)

saveRDS(CART.combined.CD8pos.newborn, file = paste('./data/', experiment, '_newborn.rds', sep = ''))


umap_seurat_clusters <- DimPlot(CART.combined.CD8pos.newborn, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters_newborn', sep = ''))

# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical_newborn <- VlnPlot(CART.combined.CD8pos.newborn, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "seurat_clusters", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical_newborn, paste('./plots/', experiment, '_vlnplot_naiveT_canonical_newborn', sep = ''))







