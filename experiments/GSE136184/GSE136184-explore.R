# explore seurat object
# Vandon Duong

experiment = 'GSE136184'

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
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "AgeGroup")

head(CART.combined.CD8pos)


####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_age_group <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "AgeGroup", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group, paste('./plots/', experiment, '_umap_age_group', sep = ''))

umap_age_group_2 <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "AgeGroup2", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group_2, paste('./plots/', experiment, '_umap_age_group_2', sep = ''))

umap_sex <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Sex", shuffle = TRUE, seed = 123)
generate_figs(umap_sex, paste('./plots/', experiment, '_umap_sex', sep = ''))

umap_visit <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "visit", shuffle = TRUE, seed = 123)
generate_figs(umap_visit, paste('./plots/', experiment, '_umap_visit', sep = ''))

umap_code <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Code", shuffle = TRUE, seed = 123)
generate_figs(umap_code, paste('./plots/', experiment, '_umap_code', sep = ''))

umap_cohort <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Cohort", shuffle = TRUE, seed = 123)
generate_figs(umap_cohort, paste('./plots/', experiment, '_umap_cohort', sep = ''))

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- CART.combined.CD8pos@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup", "Phase")]
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

CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "Cohort")

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
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "AgeGroup")
vlnplot_CARTEx_84_age_group <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'AgeGroup', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_age_group, paste('./plots/', experiment, '_vlnplot_CARTEx_84_age_group', sep = ''))

vlnplot_CARTEx_84_age_group_2 <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'AgeGroup2', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_age_group_2, paste('./plots/', experiment, '_vlnplot_CARTEx_84_age_group_2', sep = ''))

# young vs old
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "AgeGroup")
vlnplot_CARTEx_84_young_old <- VlnPlot(subset(CART.combined.CD8pos, idents = c("Young", "Old")), features = c("CARTEx_84"), group.by = 'AgeGroup', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Young','Old')), label = "p.signif", label.y = 4)
generate_figs(vlnplot_CARTEx_84_young_old, paste('./plots/', experiment, '_vlnplot_CARTEx_84_young_old', sep = ''))



# analyze only cross-sectional / longitudinal
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "Cohort")
CART.combined.CD8pos.cross <- subset(CART.combined.CD8pos, idents = c("Cross-sectional"))
CART.combined.CD8pos.long <- subset(CART.combined.CD8pos, idents = c("Longitudinal"))
rm(CART.combined.CD8pos)

CART.combined.CD8pos.cross <- SetIdent(CART.combined.CD8pos.cross, value = "AgeGroup")
CART.combined.CD8pos.long <- SetIdent(CART.combined.CD8pos.long, value = "AgeGroup")

vlnplot_CARTEx_84_young_old_cross <- VlnPlot(subset(CART.combined.CD8pos.cross, idents = c("Young", "Old")), features = c("CARTEx_84"), group.by = 'AgeGroup', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Young','Old')), label = "p.signif", label.y = 4)
generate_figs(vlnplot_CARTEx_84_young_old_cross, paste('./plots/', experiment, '_vlnplot_CARTEx_84_young_old_cross', sep = ''))

vlnplot_CARTEx_84_young_old_long <- VlnPlot(subset(CART.combined.CD8pos.long, idents = c("Young", "Old")), features = c("CARTEx_84"), group.by = 'AgeGroup', y.max = 5, pt.size = 0) + theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Young','Old')), label = "p.signif", label.y = 4)
generate_figs(vlnplot_CARTEx_84_young_old_long, paste('./plots/', experiment, '_vlnplot_CARTEx_84_young_old_long', sep = ''))

CART.combined.CD8pos.long <- SetIdent(CART.combined.CD8pos.long, value = "Code")

vlnplot_CARTEx_84_donor_long <- VlnPlot(CART.combined.CD8pos.long, features = c("CARTEx_84"), group.by = 'Code', split.by = 'visit', y.max = 5, pt.size = 0)
generate_figs(vlnplot_CARTEx_84_donor_long, paste('./plots/', experiment, '_vlnplot_CARTEx_84_donor_long', sep = ''))



umap_age_group_cross <- DimPlot(CART.combined.CD8pos.cross, reduction = "umap", group.by = "AgeGroup", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group_cross, paste('./plots/', experiment, '_umap_age_group_cross', sep = ''))

umap_age_group_2_cross <- DimPlot(CART.combined.CD8pos.cross, reduction = "umap", group.by = "AgeGroup2", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group_2_cross, paste('./plots/', experiment, '_umap_age_group_2_cross', sep = ''))

umap_sex_cross <- DimPlot(CART.combined.CD8pos.cross, reduction = "umap", group.by = "Sex", shuffle = TRUE, seed = 123)
generate_figs(umap_sex_cross, paste('./plots/', experiment, '_umap_sex_cross', sep = ''))

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- CART.combined.CD8pos.cross@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup", "Phase")]
pie_phase <- md[, .N, by = c("Phase")]
setorder(pie_phase, cols = "Phase")
pie_phase$percent <- round(100*pie_phase$N / sum(pie_phase$N), digits = 1)
pie_phase$label <- paste(pie_phase$Phase, " (", pie_phase$percent,"%)", sep = "")
pie_phase$cols <- hcl.colors(n = nrow(pie_phase), palette = "Temps")
barplot_phase_cross <- ggplot(data= pie_phase, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = pie_phase$cols)
generate_figs(barplot_phase_cross, paste('./plots/', experiment, '_barplot_phase_cross', sep = ''))
umap_phase_cross <- DimPlot(CART.combined.CD8pos.cross, group.by = "Phase", cols = pie_phase$cols)
generate_figs(umap_phase_cross, paste('./plots/', experiment, '_umap_phase_cross', sep = ''))

umap_CARTEx_84_cross <- FeaturePlot(CART.combined.CD8pos.cross, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_CARTEx_200_cross <- FeaturePlot(CART.combined.CD8pos.cross, features = c("CARTEx_200"), order = TRUE) + fix.sc
umap_CARTEx_630_cross <- FeaturePlot(CART.combined.CD8pos.cross, features = c("CARTEx_630"), order = TRUE) + fix.sc

generate_figs(umap_CARTEx_84_cross, paste('./plots/', experiment, '_umap_CARTEx_84_cross', sep = ''))
generate_figs(umap_CARTEx_200_cross, paste('./plots/', experiment, '_umap_CARTEx_200_cross', sep = ''))
generate_figs(umap_CARTEx_630_cross, paste('./plots/', experiment, '_umap_CARTEx_630_cross', sep = ''))









