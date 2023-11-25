# integrate seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE125881'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)

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

integerize = function(score){
  score_mod = round(score)
  score_mod[score_mod < -4] <- -5
  score_mod[score_mod > 4] <- 5
  return (score_mod)
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################


# https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_rpca.html

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_scored_cross.rds")
ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_CD8pos_scored.rds")

# add experiment identifier
expt.obj@meta.data[['identifier']] <- experiment
aging.obj@meta.data[['identifier']] <- "GSE136184"
ref.obj@meta.data[['identifier']] <- "GSE164378"

expt.obj@meta.data[['identifier2']] <- expt.obj@meta.data$Group
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

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

saveRDS(integration.obj, file = paste('./data/', experiment, '_integrated.rds', sep = ''))

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

integrated_umap_identifier <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier2", split.by = "identifier", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier, paste('./plots/', experiment, '_integrated_umap_identifier', sep = ''), c(12, 5))

# library(rlang)
# library(stringr)
# https://github.com/satijalab/seurat/issues/1396
plot_umap_highlight = function(atlas, identity){
  plot.list <- list()
  for (i in sapply(unique(x = atlas[[deparse(substitute(identity))]]), levels)) {
    plot.list[[i]] <- DimPlot(
      object = atlas, cells.highlight = Cells(integration.obj[, integration.obj[[deparse(substitute(identity))]] == i])
    ) + NoLegend() + ggtitle(i)
  }
  combined_plots <- CombinePlots(plots = plot.list, ncol = 3)
  return(combined_plots)
}

# two ways to capture cells which match
# WhichCells(object = integration.obj, expression = identifier == experiment)
# Cells(integration.obj[, integration.obj[['identifier']] == experiment])

integrated_umap_identifier_highlight <- plot_umap_highlight(integration.obj, identifier)
generate_figs(integrated_umap_identifier_highlight, paste('./plots/', experiment, '_integrated_umap_identifier_highlight', sep = ''), c(12, 5))

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
split.ident.order = c("IP", "ExpansionPeak", "Contraction", "Late", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly")
Idents(query.obj) <- factor(Idents(query.obj), levels = split.ident.order)

query.obj$identifier2 <- factor(query.obj$identifier2, levels = split.ident.order)
vlnplot_CARTEx_84 <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','ExpansionPeak')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Late')), label = "p.signif", label.y = 4.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Newborn')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('ExpansionPeak','Newborn')), label = "p.signif", label.y = 5.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Late','Newborn')), label = "p.signif", label.y = 3.5)
generate_figs(vlnplot_CARTEx_84, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84', sep = ''))

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


