# explore seurat object
# Vandon Duong

experiment = 'GSE153931'

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
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "orig.donor")

head(CART.combined.CD8pos)


####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_patient <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.donor", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_umap_patient', sep = ''))

umap_orig_virus  <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.virus", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_virus, paste('./plots/', experiment, '_umap_orig_virus', sep = ''))

umap_orig_virus2  <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.virus2", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_virus2, paste('./plots/', experiment, '_umap_orig_virus2', sep = ''))

umap_orig_severity  <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.severity", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_severity, paste('./plots/', experiment, '_umap_orig_severity', sep = ''))

umap_orig_severity_score  <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.severity_score", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_severity_score, paste('./plots/', experiment, '_umap_orig_severity_score', sep = ''))

umap_orig_peptide  <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.peptide", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_peptide, paste('./plots/', experiment, '_umap_orig_peptide', sep = ''))

umap_orig_hospital <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.hospital", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_hospital, paste('./plots/', experiment, '_umap_orig_hospital', sep = ''))

umap_orig_unit <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.unit", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_unit, paste('./plots/', experiment, '_umap_orig_unit', sep = ''))

umap_orig_diabetes <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.diabetes", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_diabetes, paste('./plots/', experiment, '_umap_orig_diabetes', sep = ''))

umap_orig_asthma <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.asthma", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_asthma, paste('./plots/', experiment, '_umap_orig_asthma', sep = ''))

umap_cardiovascular <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "cardiovascular", shuffle = TRUE, seed = 123)
generate_figs(umap_cardiovascular, paste('./plots/', experiment, '_umap_cardiovascular', sep = ''))

umap_tumour <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "tumour", shuffle = TRUE, seed = 123)
generate_figs(umap_tumour, paste('./plots/', experiment, '_umap_tumour', sep = ''))



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
umap_sig_activation <- FeaturePlot(CART.combined.CD8pos, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(CART.combined.CD8pos, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(CART.combined.CD8pos, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(CART.combined.CD8pos, features = c("Senescence"), order = TRUE) + fix.sc


generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
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


heatmap_allgenes <- DoHeatmap(CART.combined.CD8pos, group.by = "orig.severity_score", group.bar = TRUE, angle = 0)
generate_figs(heatmap_allgenes, paste('./plots/', experiment, '_heatmap_allgenes', sep = ''), c(7,10))

heatmap_CARTEx <- DoHeatmap(CART.combined.CD8pos, features = CARTEx_84, group.by = "orig.severity_score", group.bar = TRUE, angle = 0)
generate_figs(heatmap_CARTEx, paste('./plots/', experiment, '_heatmap_CARTEx', sep = ''), c(12,10))

heatmap_TSG <- DoHeatmap(CART.combined.CD8pos, features = intersect(TSG, CARTEx_630), group.by = "orig.severity_score", group.bar = TRUE, angle = 0)
generate_figs(heatmap_TSG, paste('./plots/', experiment, '_heatmap_TSG', sep = ''), c(12,10))


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

# violin plots

vlnplot_CARTEx_84_orig_severity <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'orig.severity', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_orig_severity, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_severity', sep = ''))

CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "orig.severity_score")
CART.combined.CD8pos@meta.data$orig.severity_score <- factor(CART.combined.CD8pos@meta.data$orig.severity_score, levels = c("Mild", "Moderate", "Severe", "Critical", "void"))
levels(CART.combined.CD8pos)
vlnplot_CARTEx_84_orig_severity_score <- VlnPlot(subset(CART.combined.CD8pos,cells=colnames(CART.combined.CD8pos)[Idents(CART.combined.CD8pos)!="void"]), features = c("CARTEx_84"), group.by = 'orig.severity_score', y.max = 5, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Critical','Mild')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Critical','Moderate')), label = "p.signif", label.y = 3.0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Critical','Severe')), label = "p.signif", label.y = 2.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Moderate')), label = "p.signif", label.y = 2.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Moderate','Severe')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 3.5)
generate_figs(vlnplot_CARTEx_84_orig_severity_score, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_severity_score', sep = ''))


CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "cardiovascular")
CART.combined.CD8pos@meta.data$cardiovascular <- factor(CART.combined.CD8pos@meta.data$cardiovascular, levels = c("None", "DVT", "Heart", "HTN", "void"))
levels(CART.combined.CD8pos)
vlnplot_CARTEx_84_cardiovascular <- VlnPlot(subset(CART.combined.CD8pos,cells=colnames(CART.combined.CD8pos)[Idents(CART.combined.CD8pos)!="void"]), features = c("CARTEx_84"), group.by = 'cardiovascular', y.max = 5, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','DVT')), label = "p.signif", label.y = 3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','Heart')), label = "p.signif", label.y = 3.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','HTN')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('DVT','HTN')), label = "p.signif", label.y = 2.5)
generate_figs(vlnplot_CARTEx_84_cardiovascular, paste('./plots/', experiment, '_vlnplot_CARTEx_84_cardiovascular', sep = ''))


CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "orig.diabetes")
CART.combined.CD8pos@meta.data$orig.diabetes <- factor(CART.combined.CD8pos@meta.data$orig.diabetes, levels = c("None", "Type1", "Type2", "void"))
levels(CART.combined.CD8pos)
vlnplot_CARTEx_84_orig_diabetes <- VlnPlot(subset(CART.combined.CD8pos,cells=colnames(CART.combined.CD8pos)[Idents(CART.combined.CD8pos)!="void"]), features = c("CARTEx_84"), group.by = 'orig.diabetes', y.max = 4, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','Type1')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','Type2')), label = "p.signif", label.y = 3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Type1','Type2')), label = "p.signif", label.y = 2.5)
generate_figs(vlnplot_CARTEx_84_orig_diabetes, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_diabetes', sep = ''))

CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "orig.asthma")
CART.combined.CD8pos@meta.data$orig.asthma <- factor(CART.combined.CD8pos@meta.data$orig.asthma, levels = c("None", "Asthma", "void"))
levels(CART.combined.CD8pos)
vlnplot_CARTEx_84_orig_asthma <- VlnPlot(subset(CART.combined.CD8pos,cells=colnames(CART.combined.CD8pos)[Idents(CART.combined.CD8pos)!="void"]), features = c("CARTEx_84"), group.by = 'orig.asthma', y.max = 4, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('None','Asthma')), label = "p.signif", label.y = 3)
generate_figs(vlnplot_CARTEx_84_orig_asthma, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_asthma', sep = ''))




CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "orig.unit")
CART.combined.CD8pos@meta.data$orig.unit <- factor(CART.combined.CD8pos@meta.data$orig.unit, levels = c("Not_admitted", "Ward", "ITU", "void"))
levels(CART.combined.CD8pos)
vlnplot_CARTEx_84_orig_unit <- VlnPlot(subset(CART.combined.CD8pos,cells=colnames(CART.combined.CD8pos)[Idents(CART.combined.CD8pos)!="void"]), features = c("CARTEx_84"), group.by = 'orig.unit', y.max = 3.5, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Not_admitted','Ward')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Not_admitted','ITU')), label = "p.signif", label.y = 2.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ward','ITU')), label = "p.signif", label.y = 2)
generate_figs(vlnplot_CARTEx_84_orig_unit, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_unit', sep = ''))




vlnplot_CARTEx_84_orig_hospital <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'orig.hospital', y.max = 3, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Yes','No')), label = "p.signif", label.y = 2.5)
generate_figs(vlnplot_CARTEx_84_orig_hospital, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_hospital', sep = ''))

vlnplot_CARTEx_84_orig_virus2 <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'orig.virus2', y.max = 3, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_orig_virus2, paste('./plots/', experiment, '_vlnplot_CARTEx_84_orig_virus2', sep = ''))


# crosstab with dplyr
# https://www.statology.org/dplyr-crosstab/

library(dplyr)
library(tidyr)
df_metadata <- as.data.frame(CART.combined.CD8pos@meta.data)
df_metadata <- subset(df_metadata, !duplicated(subset(df_metadata, select = c(orig.donor))))
df_metadata %>% group_by(orig.severity_score, orig.unit) %>% tally() %>% spread(orig.severity_score, n)

df_metadata <- as.data.frame(CART.combined.CD8pos@meta.data)
df_metadata <- subset(df_metadata, !duplicated(subset(df_metadata, select = c(orig.donor))))
df_metadata %>% group_by(orig.diabetes, cardiovascular) %>% tally() %>% spread(orig.diabetes, n)












