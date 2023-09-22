# explore seurat object
# Vandon Duong

experiment = 'GSE125881'

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
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "Patient")

head(CART.combined.CD8pos)


####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_patient <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Patient", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_umap_patient', sep = ''))

umap_timepoint <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "TimePoint", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint, paste('./plots/', experiment, '_umap_timepoint', sep = ''))

umap_group <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Group", shuffle = TRUE, seed = 123)
generate_figs(umap_group, paste('./plots/', experiment, '_umap_group', sep = ''))

umap_disease <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Disease", shuffle = TRUE, seed = 123)
generate_figs(umap_disease, paste('./plots/', experiment, '_umap_disease', sep = ''))


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


heatmap_allgenes <- DoHeatmap(CART.combined.CD8pos, group.by = "TimePoint", group.bar = TRUE, angle = 0)
generate_figs(heatmap_allgenes, paste('./plots/', experiment, '_heatmap_allgenes', sep = ''), c(7,10))

heatmap_CARTEx <- DoHeatmap(CART.combined.CD8pos, features = CARTEx_84, group.by = "TimePoint", group.bar = TRUE, angle = 0)
generate_figs(heatmap_CARTEx, paste('./plots/', experiment, '_heatmap_CARTEx', sep = ''), c(7,10))

heatmap_TSG <- DoHeatmap(CART.combined.CD8pos, features = intersect(TSG, CARTEx_630), group.by = "TimePoint", group.bar = TRUE, angle = 0)
generate_figs(heatmap_TSG, paste('./plots/', experiment, '_heatmap_TSG', sep = ''), c(7,10))



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



# kinetics plots

vlnplot_CARTEx_84_Group <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'Group', y.max = 10, pt.size = 0) + theme(legend.position = 'none') +
  geom_boxplot(width=0.2, color="black", alpha=0.75) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','ExpansionPeak')), label = "p.signif", label.y = 9) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('ExpansionPeak','Contraction')), label = "p.signif", label.y = 8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Contraction','Late')), label = "p.signif", label.y = 7)
generate_figs(vlnplot_CARTEx_84_Group, paste('./plots/', experiment, '_vlnplot_CARTEx_84_Group', sep = ''))


vlnplot_CARTEx_84_timepoint <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'TimePoint', y.max = 10, pt.size = 0) + theme(legend.position = 'none') +
  geom_boxplot(width=0.2, color="black", alpha=0.75) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','d12')), label = "p.signif", label.y = 9.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d12','d21')), label = "p.signif", label.y = 8.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d21','d28')), label = "p.signif", label.y = 7.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d28','d29')), label = "p.signif", label.y = 6.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d29','d38')), label = "p.signif", label.y = 5.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d38','d83')), label = "p.signif", label.y = 4.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d83','d89')), label = "p.signif", label.y = 3.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d89','d102')), label = "p.signif", label.y = 4.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('d102','d112')), label = "p.signif", label.y = 5.2)
generate_figs(vlnplot_CARTEx_84_timepoint, paste('./plots/', experiment, '_vlnplot_CARTEx_84_timepoint', sep = ''))


vlnplot_CARTEx_84_patients <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'Patient', y.max = 10, pt.size = 0) + theme(legend.position = 'none') +
  geom_boxplot(width=0.2, color="black", alpha=0.75)
generate_figs(vlnplot_CARTEx_84_patients, paste('./plots/', experiment, '_vlnplot_CARTEx_84_patients', sep = ''))


# SPLIT ANALYSIS
levels(CART.combined.CD8pos$Group)

pie_phase_IP <- md[md$Group == "IP"][, .N, by = c("Phase")]
pie_phase_IP$percent <- round(100*pie_phase_IP$N / sum(pie_phase_IP$N), digits = 1)
pie_phase_IP$label <- paste(pie_phase_IP$Phase, " (", pie_phase_IP$percent,"%)", sep = "")
pie_phase_IP <- pie_phase_IP[c(1,3,2),]
pie_phase_IP$cols <- hcl.colors(n = nrow(pie_phase_IP), palette = "Temps")
pie(pie_phase_IP$N, labels = pie_phase_IP$label, col = pie_phase_IP$cols)

pie_phase_ExpansionPeak <- md[md$Group == "ExpansionPeak"][, .N, by = c("Phase")]
pie_phase_ExpansionPeak$percent <- round(100*pie_phase_ExpansionPeak$N / sum(pie_phase_ExpansionPeak$N), digits = 1)
pie_phase_ExpansionPeak$label <- paste(pie_phase_ExpansionPeak$Phase, " (", pie_phase_ExpansionPeak$percent,"%)", sep = "")
pie_phase_ExpansionPeak <- pie_phase_ExpansionPeak[c(1,3,2),]
pie_phase_ExpansionPeak$cols <- hcl.colors(n = nrow(pie_phase_ExpansionPeak), palette = "Temps")
pie(pie_phase_ExpansionPeak$N, labels = pie_phase_ExpansionPeak$label, col = pie_phase_ExpansionPeak$cols)

pie_phase_Contraction <- md[md$Group == "Contraction"][, .N, by = c("Phase")]
pie_phase_Contraction$percent <- round(100*pie_phase_Contraction$N / sum(pie_phase_Contraction$N), digits = 1)
pie_phase_Contraction$label <- paste(pie_phase_Contraction$Phase, " (", pie_phase_Contraction$percent,"%)", sep = "")
pie_phase_Contraction <- pie_phase_Contraction[c(1,3,2),]
pie_phase_Contraction$cols <- hcl.colors(n = nrow(pie_phase_Contraction), palette = "Temps")
pie(pie_phase_Contraction$N, labels = pie_phase_Contraction$label, col = pie_phase_Contraction$cols)

pie_phase_Late <- md[md$Group == "Late"][, .N, by = c("Phase")]
pie_phase_Late$percent <- round(100*pie_phase_Late$N / sum(pie_phase_Late$N), digits = 1)
pie_phase_Late$label <- paste(pie_phase_Late$Phase, " (", pie_phase_Late$percent,"%)", sep = "")
pie_phase_Late <- pie_phase_Late[c(1,3,2),]
pie_phase_Late$cols <- hcl.colors(n = nrow(pie_phase_Late), palette = "Temps")
pie(pie_phase_Late$N, labels = pie_phase_Late$label, col = pie_phase_Late$cols)

pie_phase_IP$Group <- c(rep("IP",3))
pie_phase_ExpansionPeak$Group <- c(rep("ExpansionPeak",3))
pie_phase_Contraction$Group <- c(rep("Contraction",3))
pie_phase_Late$Group <- c(rep("Late",3))

pie_agg <- rbind(pie_phase_IP, pie_phase_ExpansionPeak, pie_phase_Contraction, pie_phase_Late)
pie_agg$Group <- factor(pie_agg$Group, levels = c("IP","ExpansionPeak", "Contraction", "Late"))


barplot_phase_group <- ggplot(data = pie_agg, aes(fill = Phase, y = percent, x = Group)) +
  geom_bar(position="fill", stat = "identity") + scale_fill_manual(values = c("#089392","#EAE29C","#CF597E"))
generate_figs(barplot_phase_group, paste('./plots/', experiment, '_barplot_phase_group', sep = ''))



# examine infusion product
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "Group")

obj.group <- SplitObject(CART.combined.CD8pos, split.by = "Group")

vlnplot_CARTEx_84_patients_IP <- VlnPlot(obj.group$IP, features = c("CARTEx_84"), group.by = 'Patient', y.max = 4, pt.size = 0) + theme(legend.position = 'none') +
  geom_boxplot(width=0.2, color="black", alpha=0.75)
generate_figs(vlnplot_CARTEx_84_patients_IP, paste('./plots/', experiment, '_vlnplot_CARTEx_84_patients_IP', sep = ''))













