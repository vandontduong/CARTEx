# aggregate seurat object
# Vandon Duong

set.seed(123)

experiment = 'infusionproduct'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(SingleR)
library(scuttle)


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

DimPlotHighlightIdents <- function(atlas, identity, reduction_map, highlight_color, pt_size, ncols){
  plot.list <- list()
  for (i in sapply(unique(x = atlas[[deparse(substitute(identity))]]), levels)) {
    plot.list[[i]] <- DimPlot(
      object = atlas, reduction = reduction_map, raster = FALSE, cols.highlight = highlight_color, pt.size = pt_size, sizes.highlight = pt_size,
      cells.highlight = Cells(atlas[, atlas[[deparse(substitute(identity))]] == i])
    ) + NoLegend() + ggtitle(i)
  }
  # combined_plots <- CombinePlots(plots = plot.list, ncol = ncols)
  combined_plots <- Reduce(`+`, plot.list) + patchwork::plot_layout(ncol = ncols)
  return(combined_plots)
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_rpca.html

GSE151511.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE151511/data/GSE151511_annotated.rds")
GSE125881.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE125881/data/GSE125881_annotated.rds")
# GSE120575.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE120575/data/GSE120575_annotated.rds")
aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_annotated_cross.rds")
ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_CD8pos_annotated.rds")

# add experiment identifier
GSE151511.obj@meta.data[['identifier']] <- "GSE151511"
GSE125881.obj@meta.data[['identifier']] <- "GSE125881"
aging.obj@meta.data[['identifier']] <- "GSE136184"
ref.obj@meta.data[['identifier']] <- "GSE164378"

GSE151511.obj@meta.data[['identifier2']] <- GSE151511.obj@meta.data$Responder
GSE125881.obj@meta.data[['identifier2']] <- GSE125881.obj@meta.data$Group
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

GSE151511.obj@meta.data[['identifier_IP']] <- "IP"
GSE125881.obj@meta.data[['identifier_IP']] <- plyr::mapvalues(x = GSE125881.obj@meta.data[['Group']],
                                          from = c('IP', 'ExpansionPeak', 'Contraction', 'Late'),
                                          to = c('IP', 'Other', 'Other', 'Other'))
GSE125881.obj@meta.data[['identifier_IP']] <- factor(GSE125881.obj@meta.data[['identifier_IP']], levels = c('IP', 'Other'))
aging.obj@meta.data[['identifier_IP']] <- "Other"
ref.obj@meta.data[['identifier_IP']] <- "Other"

GSE151511.obj@meta.data[['identifier_IP2']] <- "IP_GSE151511"
GSE125881.obj@meta.data[['identifier_IP2']] <- plyr::mapvalues(x = GSE125881.obj@meta.data[['Group']],
                                                              from = c('IP', 'ExpansionPeak', 'Contraction', 'Late'),
                                                              to = c('IP_GSE125881', 'Other', 'Other', 'Other'))
GSE125881.obj@meta.data[['identifier_IP2']] <- factor(GSE125881.obj@meta.data[['identifier_IP2']], levels = c('IP_GSE125881', 'Other'))
aging.obj@meta.data[['identifier_IP2']] <- "Other"
ref.obj@meta.data[['identifier_IP2']] <- "Other"

GSE151511.obj@meta.data[["split.ident"]] <- "Query"
GSE125881.obj@meta.data[["split.ident"]] <- "Query"
aging.obj@meta.data[["split.ident"]] <- "Query"
ref.obj@meta.data[["split.ident"]] <- "Reference"

integration.list <- c(GSE151511 = GSE151511.obj, GSE125881 = GSE125881.obj, aging = aging.obj, reference = ref.obj)

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
integration.anchors <- FindIntegrationAnchors(object.list = integration.list, anchor.features = integration.features, reduction = "rpca", reference = c(4), k.anchor = 20) # k.anchor = 100, k.filter = NA, dims = 1:10
integration.obj <- IntegrateData(anchorset = integration.anchors)
DefaultAssay(integration.obj) <- "integrated"

# Run the standard workflow for visualization and clustering
integration.obj <- ScaleData(integration.obj)
integration.obj <- RunPCA(integration.obj, npcs = 30)
integration.obj <- RunUMAP(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindNeighbors(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindClusters(integration.obj, resolution = 0.5)

head(integration.obj)

# add levels
unique(integration.obj$identifier)
integration.obj$identifier <- factor(integration.obj$identifier, levels = c("GSE151511", "GSE125881", "GSE136184", "GSE164378"))

unique(integration.obj$identifier_IP)
integration.obj$identifier_IP <- factor(integration.obj$identifier_IP, levels = c("IP", "Other"))

unique(integration.obj$identifier_IP2)
integration.obj$identifier_IP2 <- factor(integration.obj$identifier_IP2, levels = c("IP_GSE151511", "IP_GSE125881", "Other"))


saveRDS(integration.obj, file = paste('./data/', experiment, '_integrated.rds', sep = ''))

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

integrated_umap_identifier <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier2", split.by = "identifier", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier, paste('./plots/', experiment, '_integrated_umap_identifier', sep = ''), c(12, 5))

integrated_umap_identifier_highlight <- DimPlotHighlightIdents(integration.obj, identifier, 'umap', 'blue', 0.1, 2)
generate_figs(integrated_umap_identifier_highlight, paste('./plots/', experiment, '_integrated_umap_identifier_highlight', sep = ''), c(22, 20))

integrated_umap_identifier_expts <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier", split.by = "split.ident", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier_expts, paste('./plots/', experiment, '_integrated_umap_identifier_expts', sep = ''), c(12, 5))

integrated_umap_identifier_IP_highlight <- DimPlotHighlightIdents(integration.obj, identifier_IP, 'umap', 'blue', 0.1, 2)
generate_figs(integrated_umap_identifier_IP_highlight, paste('./plots/', experiment, '_integrated_umap_identifier_IP_highlight', sep = ''), c(10, 6))

integrated_umap_identifier_IP2_highlight <- DimPlotHighlightIdents(integration.obj, identifier_IP2, 'umap', 'blue', 0.1, 3)
generate_figs(integrated_umap_identifier_IP2_highlight, paste('./plots/', experiment, '_integrated_umap_identifier_IP2_highlight', sep = ''), c(15, 6))

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






