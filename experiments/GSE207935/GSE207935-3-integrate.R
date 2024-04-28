#!/usr/bin/env Rscript

### Script name: GSE207935-3-integrate.R
### Description: integrate seurat object for GSE207935
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE207935'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_rpca.html

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
# aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_scored_cross.rds")
aging.obj <- readRDS(paste(PATH_EXPERIMENTS, "GSE136184/data/GSE136184_aging_extract.rds", sep = ''))
ref.obj <- readRDS(paste(PATH_EXPERIMENTS, "GSE164378/data/GSE164378_CD8pos_scored.rds", sep = ''))

# add experiment identifier
expt.obj@meta.data[['identifier']] <- experiment
aging.obj@meta.data[['identifier']] <- "GSE136184"
ref.obj@meta.data[['identifier']] <- "GSE164378"

expt.obj@meta.data[['identifier2']] <- expt.obj@meta.data$Affstat
# aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[['identifier3']] <- expt.obj@meta.data$Stim
aging.obj@meta.data[['identifier3']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier3']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[["split.ident"]] <- "Query"
aging.obj@meta.data[["split.ident"]] <- "Query"
ref.obj@meta.data[["split.ident"]] <- "Reference"

integration.list <- c(experiment = expt.obj, aging = aging.obj, reference = ref.obj)

# normalize and identify variable features for each dataset independently
integration.list <- lapply(X = integration.list, FUN = function(x) {
  x <- UpdateSeuratObject(x)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

integration.features <- SelectIntegrationFeatures(object.list = integration.list)

integration.list <- lapply(X = integration.list, FUN = function(x) {
  x <- ScaleData(x, features = integration.features)
  x <- RunPCA(x, features = integration.features)
})

# integration.list <- PrepSCTIntegration(object.list = integration.list, anchor.features = integration.features)
integration.anchors <- FindIntegrationAnchors(object.list = integration.list, anchor.features = integration.features, reduction = "rpca", reference = c(3), k.anchor = 20) # k.anchor = 100, k.filter = NA, dims = 1:10
integration.obj <- IntegrateData(anchorset = integration.anchors)
DefaultAssay(integration.obj) <- "integrated"

# Run the standard workflow for visualization and clustering
integration.obj <- ScaleData(integration.obj)
integration.obj <- RunPCA(integration.obj, npcs = 30)
integration.obj <- RunUMAP(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindNeighbors(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindClusters(integration.obj, resolution = 0.5)

head(integration.obj)

integration.obj$identifier <- factor(integration.obj$identifier, levels = c(experiment, "GSE136184", "GSE164378"))
unique(integration.obj$identifier)

integration.obj$identifier2 <- factor(integration.obj$identifier2, levels = c(names(table(expt.obj$identifier2)), names(table(aging.obj$identifier2)), names(table(ref.obj$identifier2))))
table(integration.obj$identifier2)

integration.obj$identifier3 <- factor(integration.obj$identifier3, levels = c(names(table(expt.obj$identifier3)), names(table(aging.obj$identifier3)), names(table(ref.obj$identifier3))))
table(integration.obj$identifier3)

unique(integration.obj$azimuth)
integration.obj$azimuth <- factor(integration.obj$azimuth, levels = c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"))

unique(integration.obj$monaco)
integration.obj$monaco <- factor(integration.obj$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))

unique(integration.obj$dice)
integration.obj$dice <- factor(integration.obj$dice, levels = c("T cells, CD8+, naive", "T cells, CD8+, naive, stimulated"))

unique(integration.obj$split.ident)
integration.obj$split.ident <- factor(integration.obj$split.ident, levels = c("Query", "Reference"))


saveRDS(integration.obj, file = paste('./data/', experiment, '_integrated.rds', sep = ''))

# integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

integrated_umap_identifier <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier2", split.by = "identifier", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier, paste('./plots/', experiment, '_integrated_umap_identifier', sep = ''), c(12, 5))

integrated_umap_identifier_highlight <- DimPlotHighlightIdents(integration.obj, identifier, 'umap', 'blue', 0.1, 3)
generate_figs(integrated_umap_identifier_highlight, paste('./plots/', experiment, '_integrated_umap_identifier_highlight', sep = ''), c(12, 5))

integrated_umap_identifier_expts <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier", split.by = "split.ident", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier_expts, paste('./plots/', experiment, '_integrated_umap_identifier_expts', sep = ''), c(12, 5))

integrated_umap_seurat_clusters <- DimPlot(integration.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_seurat_clusters, paste('./plots/', experiment, '_integrated_umap_seurat_clusters', sep = ''))

integrated_umap_seurat_clusters_highlight <- DimPlotHighlightIdents(integration.obj, seurat_clusters, 'umap', 'blue', 0.1, 5)
generate_figs(integrated_umap_seurat_clusters_highlight, paste('./plots/', experiment, '_integrated_umap_seurat_clusters_highlight', sep = ''), c(22, 20))

integrated_umap_azimuth <- DimPlot(integration.obj, reduction = "umap", group.by = "azimuth", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_azimuth, paste('./plots/', experiment, '_integrated_umap_azimuth', sep = ''))

integrated_umap_azimuth_highlight <- DimPlotHighlightIdents(integration.obj, azimuth, 'umap', 'blue', 0.1, 2)
generate_figs(integrated_umap_azimuth_highlight, paste('./plots/', experiment, '_integrated_umap_azimuth_highlight', sep = ''), c(10, 6))

integrated_umap_monaco <- DimPlot(integration.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_monaco, paste('./plots/', experiment, '_integrated_umap_monaco', sep = ''))

integrated_umap_monaco_highlight <- DimPlotHighlightIdents(integration.obj, monaco, 'umap', 'blue', 0.1, 2)
generate_figs(integrated_umap_monaco_highlight, paste('./plots/', experiment, '_integrated_umap_monaco_highlight', sep = ''), c(10, 6))

integrated_umap_dice <- DimPlot(integration.obj, reduction = "umap", group.by = "dice", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_dice, paste('./plots/', experiment, '_integrated_umap_dice', sep = ''))

integrated_umap_dice_highlight <- DimPlotHighlightIdents(integration.obj, dice, 'umap', 'blue', 0.1, 2)
generate_figs(integrated_umap_dice_highlight, paste('./plots/', experiment, '_integrated_umap_dice_highlight', sep = ''), c(10, 6))



featureplot_Tcell_markers <- FeaturePlot(integration.obj, features = c("CD4", "CD8A", "CD8B", "PDCD1"), raster = FALSE)
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_integrated_featureplot_Tcell_markers', sep = ''))



# BAR CHARTS of CELL TYPES
md <- integration.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_dice_seurat_clusters', sep = ''))

####################################################################################################
####################################### Examine query datasets #####################################
####################################################################################################

# integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

query.list <- SplitObject(integration.obj, split.by = "split.ident")
query.obj <- query.list$Query
rm(aging.obj, expt.obj, ref.obj, integration.obj)

md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

head(query.obj)


####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_630 <- Z(scores)
query.obj@meta.data$CARTEx_630i <- integerize(query.obj@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_200 <- Z(scores)
query.obj@meta.data$CARTEx_200i <- integerize(query.obj@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_84 <- Z(scores)
query.obj@meta.data$CARTEx_84i <- integerize(query.obj@meta.data$CARTEx_84)

saveRDS(query.obj, file = paste('./data/', experiment, '_query_scored.rds', sep = ''))

# query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

head(query.obj)







