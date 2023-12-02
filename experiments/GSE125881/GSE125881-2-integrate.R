# integrate seurat object
# Vandon Duong

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(patchwork)
# https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################


# https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_rpca.html

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
# aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_scored_cross.rds")
aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_aging_extract.rds")
ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_CD8pos_scored.rds")

# add experiment identifier
expt.obj@meta.data[['identifier']] <- experiment
aging.obj@meta.data[['identifier']] <- "GSE136184"
ref.obj@meta.data[['identifier']] <- "GSE164378"

expt.obj@meta.data[['identifier2']] <- expt.obj@meta.data$Group
# aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$extract.ident
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[['identifier3']] <- expt.obj@meta.data$TimePoint
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


saveRDS(integration.obj, file = paste('./data/', experiment, '_integrated.rds', sep = ''))

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

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

integrated_umap_monaco <- DimPlot(integration.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_monaco, paste('./plots/', experiment, '_integrated_umap_monaco', sep = ''))

integrated_umap_dice <- DimPlot(integration.obj, reduction = "umap", group.by = "dice", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_dice, paste('./plots/', experiment, '_integrated_umap_dice', sep = ''))

unique(integration.obj$dice)
integration.obj$dice <- factor(integration.obj$dice, levels = c("T cells, CD8+, naive", "T cells, CD8+, naive, stimulated"))
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

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

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


####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

query.obj <- AddModuleScore(query.obj, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
query.obj@meta.data$Activation <- scale(query.obj@meta.data$State1)
query.obj@meta.data$Anergy <- scale(query.obj@meta.data$State2)
query.obj@meta.data$Stemness <- scale(query.obj@meta.data$State3)
query.obj@meta.data$Senescence <- scale(query.obj@meta.data$State4)

query.obj@meta.data$Activationi <- integerize(query.obj@meta.data$Activation)
query.obj@meta.data$Anergyi <- integerize(query.obj@meta.data$Anergy)
query.obj@meta.data$Stemnessi <- integerize(query.obj@meta.data$Stemness)
query.obj@meta.data$Senescencei <- integerize(query.obj@meta.data$Senescence)

query.obj@meta.data$State1 <- NULL
query.obj@meta.data$State2 <- NULL
query.obj@meta.data$State3 <- NULL
query.obj@meta.data$State4 <- NULL

saveRDS(query.obj, file = paste('./data/', experiment, '_query_scored.rds', sep = ''))

query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

head(query.obj)


# CARTEx violin plot
query.obj <- SetIdent(query.obj, value = "identifier2")
split.ident.order = c("IP", "ExpansionPeak", "Contraction", "Late", "Young_Naive", "Elderly_Terminal")
Idents(query.obj) <- factor(Idents(query.obj), levels = split.ident.order)

query.obj$identifier2 <- factor(query.obj$identifier2, levels = split.ident.order)
vlnplot_CARTEx_84 <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','ExpansionPeak')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Late')), label = "p.signif", label.y = 4.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Young_Naive')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('ExpansionPeak','Young_Naive')), label = "p.signif", label.y = 5.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Late','Young_Naive')), label = "p.signif", label.y = 3.5)
generate_figs(vlnplot_CARTEx_84, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84', sep = ''))


query.obj$identifier3 <- factor(query.obj$identifier3, levels = c("IP", "d12", "d21", "d28", "d29", "d38", "d83", "d89", "d102", "d112", "Young_Naive", "Elderly_Terminal"))
vlnplot_CARTEx_84_identifier3 <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier3', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_identifier3, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_identifier3', sep = ''))





# BAR CHARTS of CELL TYPES
md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

barplot_azimuth_identifier2 <- BarPlotStackSplit(query.obj, 'azimuth', 'identifier2')
generate_figs(barplot_azimuth_identifier2, paste('./plots/', experiment, '_query_barplot_azimuth_identifier2', sep = ''), c(8,4))

barplot_monaco_identifier2 <- BarPlotStackSplit(query.obj, 'monaco', 'identifier2')
generate_figs(barplot_monaco_identifier2, paste('./plots/', experiment, '_query_barplot_monaco_identifier2', sep = ''), c(8,4))

barplot_dice_identifier2 <- BarPlotStackSplit(query.obj, 'dice', 'identifier2')
generate_figs(barplot_dice_identifier2, paste('./plots/', experiment, '_query_barplot_dice_identifier2', sep = ''), c(8,4))



vlnplot_CARTEx_84_azimuth <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'azimuth', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_azimuth, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_azimuth', sep = ''))

vlnplot_CARTEx_84_monaco <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'monaco', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_monaco, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_monaco', sep = ''))

vlnplot_CARTEx_84_dice <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'dice', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_dice, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84_dice', sep = ''))


