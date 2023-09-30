# explore seurat object
# Vandon Duong

experiment = 'GSE164378'

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

pbmcref <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
pbmcref <- SetIdent(pbmcref, value = "celltype.l2")

head(pbmcref)


####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(pbmcref, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_celltype_l1 <- DimPlot(pbmcref, reduction = "umap", group.by = "celltype.l1", shuffle = TRUE, seed = 123)
generate_figs(umap_celltype_l1, paste('./plots/', experiment, '_umap_celltype_l1', sep = ''))

umap_celltype_l2 <- DimPlot(pbmcref, reduction = "umap", group.by = "celltype.l2", shuffle = TRUE, seed = 123)
generate_figs(umap_celltype_l2, paste('./plots/', experiment, '_umap_celltype_l2', sep = ''))

umap_celltype_l3 <- DimPlot(pbmcref, reduction = "umap", group.by = "celltype.l3", shuffle = TRUE, seed = 123)
generate_figs(umap_celltype_l3, paste('./plots/', experiment, '_umap_celltype_l3', sep = ''))

umap_donor <- DimPlot(pbmcref, reduction = "umap", group.by = "donor", shuffle = TRUE, seed = 123)
generate_figs(umap_donor, paste('./plots/', experiment, '_umap_donor', sep = ''))

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- pbmcref@meta.data %>% as.data.table
md[, .N, by = c("Phase", "Phase_ref")]
pie_phase <- md[, .N, by = c("Phase")]
setorder(pie_phase, cols = "Phase")
pie_phase$percent <- round(100*pie_phase$N / sum(pie_phase$N), digits = 1)
pie_phase$label <- paste(pie_phase$Phase, " (", pie_phase$percent,"%)", sep = "")
pie_phase$cols <- hcl.colors(n = nrow(pie_phase), palette = "Temps")
barplot_phase <- ggplot(data= pie_phase, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = pie_phase$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_barplot_phase', sep = ''))
umap_phase <- DimPlot(pbmcref, group.by = "Phase", cols = pie_phase$cols)
generate_figs(umap_phase, paste('./plots/', experiment, '_umap_phase', sep = ''))

pie_phase_ref <- md[, .N, by = c("Phase_ref")]
setorder(pie_phase_ref, cols = "Phase_ref")
pie_phase_ref$percent <- round(100*pie_phase_ref$N / sum(pie_phase_ref$N), digits = 1)
pie_phase_ref$label <- paste(pie_phase_ref$Phase_ref, " (", pie_phase_ref$percent,"%)", sep = "")
pie_phase_ref$cols <- hcl.colors(n = nrow(pie_phase_ref), palette = "Temps")
barplot_phase_ref <- ggplot(data= pie_phase_ref, aes(x=Phase_ref, y=percent)) + geom_bar(stat="identity", fill = pie_phase_ref$cols)
generate_figs(barplot_phase_ref, paste('./plots/', experiment, '_barplot_phase_ref', sep = ''))
umap_phase_ref <- DimPlot(pbmcref, group.by = "Phase_ref", cols = pie_phase_ref$cols)
generate_figs(umap_phase_ref, paste('./plots/', experiment, '_umap_phase_ref', sep = ''))

# UMAP of scores
umap_CARTEx_630i <- DimPlot(pbmcref, group.by = "CARTEx_630i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_630i, paste('./plots/', experiment, '_umap_CARTEx_630i', sep = ''))

umap_CARTEx_200i <- DimPlot(pbmcref, group.by = "CARTEx_200i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_200i, paste('./plots/', experiment, '_umap_CARTEx_200i', sep = ''))

umap_CARTEx_84i <- DimPlot(pbmcref, group.by = "CARTEx_84i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_84i, paste('./plots/', experiment, '_umap_CARTEx_84i', sep = ''))

umap_sig_activationi <- DimPlot(pbmcref, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_umap_sig_activationi', sep = ''))

umap_sig_anergyi <- DimPlot(pbmcref, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_umap_sig_anergyi', sep = ''))

umap_sig_stemnessi <- DimPlot(pbmcref, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_umap_sig_stemnessi', sep = ''))

umap_sig_senescencei <- DimPlot(pbmcref, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_umap_sig_senescencei', sep = ''))

# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(pbmcref, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_CARTEx_200 <- FeaturePlot(pbmcref, features = c("CARTEx_200"), order = TRUE) + fix.sc
umap_CARTEx_630 <- FeaturePlot(pbmcref, features = c("CARTEx_630"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(pbmcref, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(pbmcref, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(pbmcref, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(pbmcref, features = c("Senescence"), order = TRUE) + fix.sc


generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
generate_figs(umap_CARTEx_200, paste('./plots/', experiment, '_umap_CARTEx_200', sep = ''))
generate_figs(umap_CARTEx_630, paste('./plots/', experiment, '_umap_CARTEx_630', sep = ''))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_umap_sig_activation', sep = ''))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_umap_sig_anergy', sep = ''))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_umap_sig_stemness', sep = ''))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_umap_sig_senescence', sep = ''))









# Violin plots

pbmcref_CD8T <- subset(x = pbmcref, idents = c("CD8 TEM", "CD8 TCM", "CD8 Naive", "CD8 Proliferating"))

vlnplot_CARTEx_84_celltype_l2 <- VlnPlot(pbmcref_CD8T, features = c("CARTEx_84"), group.by = 'celltype.l2', y.max = 4, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_celltype_l2, paste('./plots/', experiment, '_vlnplot_CARTEx_84_celltype_l2', sep = ''))





