#!/usr/bin/env Rscript

### Script name: GSE151511-3-integrate.R
### Description: integrate seurat object for GSE151511
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE151511'
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

expt.obj@meta.data[['identifier2']] <- expt.obj@meta.data$Responder
# aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[['identifier2_CRS']] <- expt.obj@meta.data$CRS
aging.obj@meta.data[['identifier2_CRS']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier2_CRS']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[['identifier2_ICANS']] <- expt.obj@meta.data$ICANS
aging.obj@meta.data[['identifier2_ICANS']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier2_ICANS']] <- ref.obj@meta.data$monaco

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
integration.anchors <- FindIntegrationAnchors(object.list = integration.list, anchor.features = integration.features, reduction = "rpca", reference = c(3)) # k.anchor = 100, k.filter = NA, dims = 1:10
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

integration.obj$identifier2_CRS <- factor(integration.obj$identifier2_CRS, levels = c(names(table(expt.obj$identifier2_CRS)), names(table(aging.obj$identifier2_CRS)), names(table(ref.obj$identifier2_CRS))))
table(integration.obj$identifier2_CRS)

integration.obj$identifier2_ICANS <- factor(integration.obj$identifier2_ICANS, levels = c(names(table(expt.obj$identifier2_ICANS)), names(table(aging.obj$identifier2_ICANS)), names(table(ref.obj$identifier2_ICANS))))
table(integration.obj$identifier2_ICANS)

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

integrated_umap_identifier_expts <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier", split.by = "split.ident", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier_expts, paste('./plots/', experiment, '_integrated_umap_identifier_expts', sep = ''), c(12, 5))

integrated_umap_seurat_clusters <- DimPlot(integration.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_seurat_clusters, paste('./plots/', experiment, '_integrated_umap_seurat_clusters', sep = ''))

integrated_umap_azimuth <- DimPlot(integration.obj, reduction = "umap", group.by = "azimuth", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_azimuth, paste('./plots/', experiment, '_integrated_umap_azimuth', sep = ''))

integrated_umap_monaco <- DimPlot(integration.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_monaco, paste('./plots/', experiment, '_integrated_umap_monaco', sep = ''))

integrated_umap_dice <- DimPlot(integration.obj, reduction = "umap", group.by = "dice", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_dice, paste('./plots/', experiment, '_integrated_umap_dice', sep = ''))

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
unique(query.obj$identifier2)

# CARTEx violin plot
query.obj <- SetIdent(query.obj, value = "identifier2")
split.ident.order = c("Responder", "Non-responder", "NE", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly")
Idents(query.obj) <- factor(Idents(query.obj), levels = split.ident.order)

query.obj$identifier2 <- factor(query.obj$identifier2, levels = split.ident.order)
vlnplot_CARTEx_84 <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Responder','Non-responder')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Responder','Newborn')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Non-responder','Newborn')), label = "p.signif", label.y = 4.5)
generate_figs(vlnplot_CARTEx_84, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84', sep = ''))


unique(query.obj$identifier2_CRS)
query.obj$identifier2_CRS <- factor(query.obj$identifier2_CRS, levels = c("1", "2", "3", "4", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
vlnplot_CARTEx_84_CRS <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2_CRS', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('1','4')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('1','Newborn')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('4','Newborn')), label = "p.signif", label.y = 4.5)
generate_figs(vlnplot_CARTEx_84_CRS, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_CRS', sep = ''))

unique(query.obj$identifier2_ICANS)
query.obj$identifier2_ICANS <- factor(query.obj$identifier2_ICANS, levels = c("0", "1", "2", "3", "4", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
vlnplot_CARTEx_84_ICANS <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2_ICANS', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('0','4')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('0','Newborn')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('4','Newborn')), label = "p.signif", label.y = 4.5)
generate_figs(vlnplot_CARTEx_84_ICANS, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_ICANS', sep = ''))



# BAR CHARTS of CELL TYPES
md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_identifier2, paste('./plots/', experiment, '_query_barplot_azimuth_identifier2', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_identifier2, paste('./plots/', experiment, '_query_barplot_monaco_identifier2', sep = ''))

md_temp <- md[, .N, by = c('dice', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_identifier2, paste('./plots/', experiment, '_query_barplot_dice_identifier2', sep = ''))



vlnplot_CARTEx_84_azimuth <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'azimuth', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_azimuth, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_azimuth', sep = ''))

vlnplot_CARTEx_84_monaco <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'monaco', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_monaco, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_monaco', sep = ''))

vlnplot_CARTEx_84_dice <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'dice', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_dice, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_dice', sep = ''))


# other metadata

md_temp <- md[, .N, by = c('azimuth', 'identifier2_CRS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_identifier2_CRS <- ggplot(md_temp, aes(x = identifier2_CRS, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_identifier2_CRS, paste('./plots/', experiment, '_query_barplot_azimuth_identifier2_CRS', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'identifier2_CRS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_identifier2_CRS <- ggplot(md_temp, aes(x = identifier2_CRS, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_identifier2_CRS, paste('./plots/', experiment, '_query_barplot_monaco_identifier2_CRS', sep = ''))

md_temp <- md[, .N, by = c('dice', 'identifier2_CRS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_identifier2_CRS <- ggplot(md_temp, aes(x = identifier2_CRS, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_identifier2_CRS, paste('./plots/', experiment, '_query_barplot_dice_identifier2_CRS', sep = ''))

md_temp <- md[, .N, by = c('azimuth', 'identifier2_ICANS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_identifier2_ICANS <- ggplot(md_temp, aes(x = identifier2_ICANS, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_identifier2_ICANS, paste('./plots/', experiment, '_query_barplot_azimuth_identifier2_ICANS', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'identifier2_ICANS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_identifier2_ICANS <- ggplot(md_temp, aes(x = identifier2_ICANS, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_identifier2_ICANS, paste('./plots/', experiment, '_query_barplot_monaco_identifier2_ICANS', sep = ''))

md_temp <- md[, .N, by = c('dice', 'identifier2_ICANS')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_identifier2_ICANS <- ggplot(md_temp, aes(x = identifier2_ICANS, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_identifier2_ICANS, paste('./plots/', experiment, '_query_barplot_dice_identifier2_ICANS', sep = ''))





