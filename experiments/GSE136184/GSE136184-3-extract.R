# extract newborn data
# Vandon Duong

set.seed(123)
experiment = 'GSE136184'
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(gtools)

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_annotated_cross.rds', sep = ''))
# expt.obj <- SetIdent(expt.obj, value = "AgeGroup")
head(expt.obj)

cellIDs_newborn <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Newborn'])
cellIDs_under30 <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Under 30'])
cellIDs_elderly <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Elderly'])

cellIDs_naiveCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Naive CD8 T cells'])
cellIDs_terminalCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Terminal effector CD8 T cells'])

expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% intersect(union(cellIDs_newborn, cellIDs_under30), cellIDs_naiveCD8T))]
expt.obj.elderly_terminal <- expt.obj[,(colnames(expt.obj) %in% intersect(cellIDs_elderly, cellIDs_terminalCD8T))]

expt.obj.young_naive@meta.data$extract.ident <- "Young_Naive"
expt.obj.elderly_terminal@meta.data$extract.ident <- "Elderly_Terminal"

rm(expt.obj)

expt.obj <- merge(expt.obj.young_naive, y = expt.obj.elderly_terminal, project = "aging_extract")

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_allgenes_aging_extract.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'orig.ident', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_prepare_vlnplot_quality_aging_extract', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled_aging_extract', sep = ''))


expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

saveRDS(expt.obj, file = paste('./data/', experiment, '_aging_extract.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '_aging_extract.rds', sep = ''))
















# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical <- VlnPlot(expt.obj, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "AgeGroup2", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical, paste('./plots/', experiment, '_vlnplot_naiveT_canonical', sep = ''))

####################################################################################################
######################################## Extract newborn data ######################################
####################################################################################################


expt.obj <- SetIdent(expt.obj, value = "AgeGroup2")
expt.obj <- subset(expt.obj, idents = c("Newborn"))


# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.MT < 10)

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.MT"), group.by = 'AgeGroup2', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality_newborn', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled_newborn', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

saveRDS(expt.obj, file = paste('./data/', experiment, '_newborn.rds', sep = ''))


umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters_newborn', sep = ''))

# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical_newborn <- VlnPlot(expt.obj, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "seurat_clusters", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical_newborn, paste('./plots/', experiment, '_vlnplot_naiveT_canonical_newborn', sep = ''))







