# explore seurat object
# Vandon Duong

experiment = 'GSE196606'

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
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "group")

head(CART.combined.CD8pos)

####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_group <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "group", shuffle = TRUE, seed = 123)
generate_figs(umap_group, paste('./plots/', experiment, '_umap_group', sep = ''))

umap_cell_type <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "cellType_corrected", shuffle = TRUE, seed = 123)
generate_figs(umap_cell_type, paste('./plots/', experiment, '_umap_cell_type', sep = ''))

umap_orig_ident <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_ident, paste('./plots/', experiment, '_umap_orig_ident', sep = ''))


# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- CART.combined.CD8pos@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "Phase")]
pie_phase <- md[, .N, by = c("Phase")]
setorder(pie_phase, cols = "Phase")
pie_phase$percent <- round(100*pie_phase$N / sum(pie_phase$N), digits = 1)
pie_phase$label <- paste(pie_phase$Phase, " (", pie_phase$percent,"%)", sep = "")
pie_phase$cols <- hcl.colors(n = nrow(pie_phase), palette = "Temps")
barplot_phase <- ggplot(data= pie_phase, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = pie_phase$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_barplot_phase', sep = ''))
umap_phase <- DimPlot(CART.combined.CD8pos, group.by = "Phase", cols = pie_phase$cols)
generate_figs(umap_phase, paste('./plots/', experiment, '_umap_phase', sep = ''))


# UMAP of scores
umap_CARTEx_630i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_630i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_630i, paste('./plots/', experiment, '_umap_CARTEx_630i', sep = ''))

umap_CARTEx_200i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_200i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_200i, paste('./plots/', experiment, '_umap_CARTEx_200i', sep = ''))

umap_CARTEx_84i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_84i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_84i, paste('./plots/', experiment, '_umap_CARTEx_84i', sep = ''))

umap_sig_activationi <- DimPlot(CART.combined.CD8pos, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_umap_sig_activationi', sep = ''))

umap_sig_anergyi <- DimPlot(CART.combined.CD8pos, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_umap_sig_anergyi', sep = ''))

umap_sig_stemnessi <- DimPlot(CART.combined.CD8pos, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_umap_sig_stemnessi', sep = ''))

umap_sig_senescencei <- DimPlot(CART.combined.CD8pos, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_umap_sig_senescencei', sep = ''))

# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(CART.combined.CD8pos, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_CARTEx_200 <- FeaturePlot(CART.combined.CD8pos, features = c("CARTEx_200"), order = TRUE) + fix.sc
umap_CARTEx_630 <- FeaturePlot(CART.combined.CD8pos, features = c("CARTEx_630"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(CART.combined.CD8pos, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(CART.combined.CD8pos, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(CART.combined.CD8pos, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(CART.combined.CD8pos, features = c("Senescence"), order = TRUE) + fix.sc


generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
generate_figs(umap_CARTEx_200, paste('./plots/', experiment, '_umap_CARTEx_200', sep = ''))
generate_figs(umap_CARTEx_630, paste('./plots/', experiment, '_umap_CARTEx_630', sep = ''))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_umap_sig_activation', sep = ''))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_umap_sig_anergy', sep = ''))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_umap_sig_stemness', sep = ''))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_umap_sig_senescence', sep = ''))


####################################################################################################
############################################ More plots ###########################################
####################################################################################################

scatterplot_CARTEx_activation <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Activation', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_anergy <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Anergy', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_stemness <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Stemness', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_senescence <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Senescence', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_allstates <- CombinePlots(plots = list(scatterplot_CARTEx_activation, scatterplot_CARTEx_anergy, scatterplot_CARTEx_stemness, scatterplot_CARTEx_senescence))

generate_figs(scatterplot_CARTEx_activation, paste('./plots/', experiment, '_scatterplot_CARTEx_activation', sep = ''))
generate_figs(scatterplot_CARTEx_anergy, paste('./plots/', experiment, '_scatterplot_CARTEx_anergy', sep = ''))
generate_figs(scatterplot_CARTEx_stemness, paste('./plots/', experiment, '_scatterplot_CARTEx_stemness', sep = ''))
generate_figs(scatterplot_CARTEx_senescence, paste('./plots/', experiment, '_scatterplot_CARTEx_senescence', sep = ''))
generate_figs(scatterplot_CARTEx_allstates, paste('./plots/', experiment, '_scatterplot_CARTEx_allstates', sep = ''), c(7,5))


CARTEx_630 <- rownames(read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1))
CARTEx_84 <- rownames(read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1))
TSG <- rownames(read.csv("../../signatures/tumor-suppressor-genes.csv", row.names = 1, header = TRUE))





# Correlations of canonical exhaustion markers

canonical_exhaustion <- c('PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'TCF7', 'TOX', 'GZMB', 'IFNG', 'IL10', 'PRDM1', 'TNFRSF9')


dotplot_canonical_exhaustion <- DotPlot(CART.combined.CD8pos, features = canonical_exhaustion) + RotatedAxis()
generate_figs(dotplot_canonical_exhaustion, paste('./plots/', experiment, '_dotplot_canonical_exhaustion', sep = ''), c(10,4))

scatterplot_CARTEx_PDCD1 <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'PDCD1', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_PDCD1, paste('./plots/', experiment, '_scatterplot_CARTEx_PDCD1', sep = ''))

featureplot_canonical_exhaustion <- FeaturePlot(CART.combined.CD8pos, features = c('CARTEx_84', canonical_exhaustion), order = TRUE, ncol = 3)
generate_figs(featureplot_canonical_exhaustion, paste('./plots/', experiment, '_featureplot_canonical_exhaustion', sep = ''))


# Tumor suppressor genes

dotplot_TSG <- DotPlot(CART.combined.CD8pos, features = intersect(TSG, CARTEx_630)) + RotatedAxis()
generate_figs(dotplot_TSG, paste('./plots/', experiment, '_dotplot_TSG', sep = ''), c(18,4))






