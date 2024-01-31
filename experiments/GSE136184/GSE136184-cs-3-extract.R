#!/usr/bin/env Rscript

### Script name: GSE136184-cs-3-extract.R
### Description: extract young naive and old terminal (cross-sectional study)
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
experiment = 'GSE136184'
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_cs_scored.rds', sep = ''))

head(expt.obj)

expt.obj@meta.data$extract.ident <- paste(expt.obj@meta.data$AgeGroup, expt.obj@meta.data$monaco, sep = ", ")

ridgeplt_CARTEx_84 <- RidgePlot(expt.obj, features = c("CARTEx_84"), group.by = "AgeGroup2") + xlim(c(NA, 4))
generate_figs(ridgeplt_CARTEx_84, paste('./plots/', experiment, '_extract_ridgeplt_CARTEx_84', sep = ''), c(8,4))

ridgeplt_CARTEx_200 <- RidgePlot(expt.obj, features = c("CARTEx_200"), group.by = "AgeGroup2") + xlim(c(NA, 4))
generate_figs(ridgeplt_CARTEx_200, paste('./plots/', experiment, '_extract_ridgeplt_CARTEx_200', sep = ''), c(8,4))

ridgeplt_CARTEx_630 <- RidgePlot(expt.obj, features = c("CARTEx_630"), group.by = "AgeGroup2") + xlim(c(NA, 4))
generate_figs(ridgeplt_CARTEx_630, paste('./plots/', experiment, '_extract_ridgeplt_CARTEx_630', sep = ''), c(8,4))



cellIDs_newborn <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Newborn'])
cellIDs_under30 <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Under 30'])
cellIDs_elderly <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Elderly'])

cellIDs_naiveCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Naive CD8 T cells'])
cellIDs_terminalCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Terminal effector CD8 T cells'])

expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% intersect(union(cellIDs_newborn, cellIDs_under30), cellIDs_naiveCD8T))]
# expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% intersect(cellIDs_newborn, cellIDs_naiveCD8T))]
expt.obj.elderly_terminal <- expt.obj[,(colnames(expt.obj) %in% intersect(cellIDs_elderly, cellIDs_terminalCD8T))]

# expt.obj.young_naive@meta.data$extract.ident <- "Young, naive CD8 T cells"
# expt.obj.elderly_terminal@meta.data$extract.ident <- "Elderly, terminal effector CD8 T cells"


SelectSortCells <- function(atlas, module, direction, counts){
  md <- atlas@meta.data
  if (direction == 'high'){cells <- rownames(md[order(-md[module]), ][1:counts,])}
  if (direction == 'low'){cells <- rownames(md[order(md[module]), ][1:counts,])}
  return(cells)
}

cellIDs_YN_low <- SelectSortCells(expt.obj.young_naive, 'CARTEx_84', 'low', 3000)
cellIDs_ET_high <- SelectSortCells(expt.obj.elderly_terminal, 'CARTEx_84', 'high', 3000)
saveRDS(cellIDs_YN_low, file="./data/cellIDs_YN_low.RData")
saveRDS(cellIDs_ET_high, file="./data/cellIDs_ET_high.RData")


# expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% cellIDs_YN_low)]
# expt.obj.elderly_terminal <- expt.obj[,(colnames(expt.obj) %in% cellIDs_ET_high)]

# rm(expt.obj)

# expt.obj <- merge(expt.obj.young_naive, y = expt.obj.elderly_terminal, project = "aging_extract")

umap_extract_highlight_YN_low <- DimPlot(expt.obj, cells.highlight = cellIDs_YN_low, cols.highlight = "blue", cols = "grey", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("Selected from young naive CD8+ T cells") + theme(legend.position="none")
umap_extract_highlight_ET_high <- DimPlot(expt.obj, cells.highlight = cellIDs_ET_high, cols.highlight = "blue", cols = "grey", order = TRUE, pt.size = 0.1, sizes.highlight = 0.1) + ggtitle("Selected from elderly terminal effector CD8+ T cells") + theme(legend.position="none")
generate_figs(umap_extract_highlight_YN_low, paste('./plots/', experiment, '_cs_extract_umap_highlight_YN_low', sep = ''), c(6,5))
generate_figs(umap_extract_highlight_ET_high, paste('./plots/', experiment, '_cs_extract_umap_highlight_ET_high', sep = ''), c(6,5))

umap_extract_highlight_combined <- umap_extract_highlight_YN_low + umap_extract_highlight_ET_high
generate_figs(umap_extract_highlight_combined, paste('./plots/', experiment, '_cs_extract_umap_highlight_combined', sep = ''), c(12,5))


expt.obj <- expt.obj[,(colnames(expt.obj) %in% union(cellIDs_YN_low, cellIDs_ET_high))]

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
top10_varfeats <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_varfeats, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled_aging_extract', sep = ''))


expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

inspect_elbow <- ElbowPlot(expt.obj)
generate_figs(inspect_elbow, paste('./plots/', experiment, '_cs_extract_inspect_elbow', sep = ''), c(5, 4))

# select dimensionality based on elbowplot analysis
dim.max <- ChoosePC(expt.obj)
expt.obj <- FindNeighbors(expt.obj, dims = 1:dim.max)

# Examine a range of clustering resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.2)
expt.obj <- FindClusters(expt.obj, resolution = resolution.range)

inspect_clustering <- clustree(expt.obj, prefix = "RNA_snn_res.")
generate_figs(inspect_clustering, paste('./plots/', experiment, '_cs_extract_inspect_clustering', sep = ''), c(12, 8))

# Select clustering resolution based on clustree analysis
expt.obj@meta.data$seurat_clusters <- expt.obj@meta.data$RNA_snn_res.0.6
expt.obj@meta.data$seurat_clusters <- factor(expt.obj@meta.data$seurat_clusters, levels = SortNumStrList(unique(expt.obj@meta.data$seurat_clusters), shift = TRUE))

print("Calculating UMAPs...")
expt.obj <- RunUMAP(expt.obj, dims = 1:dim.max)

Idents(expt.obj) <- "extract.ident"
expt.obj <- RenameIdents(expt.obj, `Old, Terminal effector CD8 T cells` = "OldTerminal", `Young, Naive CD8 T cells` = "YoungNaive")
table(Idents(expt.obj))
expt.obj[['extract.ident']] <- Idents(expt.obj)
table(expt.obj[['extract.ident']])

expt.obj@meta.data$extract.ident <- factor(expt.obj@meta.data$extract.ident, levels = c('YoungNaive', 'OldTerminal'))

# count how many cells are zero
sum(expt.obj@meta.data$PFSD.CARTEx_84 == 0)

# print("Calculating diffusion maps...")
# diffusion_map <- RunDiffusion(expt.obj, k_int = 10)
# expt.obj <- diffusion_map$atlas
# saveRDS(diffusion_map$dmap, file = paste('./data/', experiment, '_cs_extract_dmap.rds', sep = ''))
# rm(diffusion_map)


saveRDS(expt.obj, file = paste('./data/', experiment, '_aging_extract.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_aging_extract.rds', sep = ''))


vlnplot_CARTEx_84_extract_ident <- VlnPlot(expt.obj, c("CARTEx_84"), group.by = "extract.ident", pt.size = 1)
generate_figs(vlnplot_CARTEx_84_extract_ident, paste('./plots/', experiment, '_extract_vlnplot_CARTEx_84_extract_ident', sep = ''), c(6,5))

vlnplot_CARTEx_84_age_group_2 <- VlnPlot(expt.obj, c("CARTEx_84"), group.by = "AgeGroup2", pt.size = 1)
generate_figs(vlnplot_CARTEx_84_age_group_2, paste('./plots/', experiment, '_extract_vlnplot_CARTEx_84_age_group_2', sep = ''), c(6,5))

vlnplot_CARTEx_630_extract_ident <- VlnPlot(expt.obj, c("CARTEx_630"), group.by = "extract.ident", pt.size = 0, cols=c("royalblue", "orchid")) + theme(legend.position = 'none') + ylim(c(-3, 6)) + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
vlnplot_CARTEx_200_extract_ident <- VlnPlot(expt.obj, c("CARTEx_200"), group.by = "extract.ident", pt.size = 0, cols=c("royalblue", "orchid")) + theme(legend.position = 'none') + ylim(c(-3, 6)) + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
vlnplot_CARTEx_84_extract_ident <- VlnPlot(expt.obj, c("CARTEx_84"), group.by = "extract.ident", pt.size = 0, cols=c("royalblue", "orchid")) + theme(legend.position = 'none') + ylim(c(-3, 6)) + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)

vlnplot_CARTEx_combined_extract_ident <- (vlnplot_CARTEx_630_extract_ident | vlnplot_CARTEx_200_extract_ident | vlnplot_CARTEx_84_extract_ident)
generate_figs(vlnplot_CARTEx_combined_extract_ident, paste('./plots/', experiment, '_extract_vlnplot_CARTEx_combined_extract_ident', sep = ''), c(14,5))


featplot_CARTEx_630_extract_ident <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', cols=c("orchid", "royalblue")) + theme(legend.position = 'none') + xlim(c(0, 11.5)) + ylim(c(-3, 6))
featplot_CARTEx_200_extract_ident <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', cols=c("orchid", "royalblue")) + theme(legend.position = 'none') + xlim(c(0, 11.5)) + ylim(c(-3, 6))
featplot_CARTEx_84_extract_ident <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', cols=c("orchid", "royalblue")) + theme(legend.position = 'none') + xlim(c(0, 11.5)) + ylim(c(-3, 6))

featplot_CARTEx_combined_extract_ident <- (featplot_CARTEx_630_extract_ident | featplot_CARTEx_200_extract_ident | featplot_CARTEx_84_extract_ident)
generate_figs(featplot_CARTEx_combined_extract_ident, paste('./plots/', experiment, '_extract_featplot_CARTEx_combined_extract_ident', sep = ''), c(14,5))


# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()






