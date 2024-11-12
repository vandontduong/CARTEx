#!/usr/bin/env Rscript

### Script name: mdandersonTCM-3-explore.R
### Description: explore seurat object for mdandersonTCM
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'mdandersonTCM'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# H, normal tissues from healthy donors; U, tumor-adjacent uninvolved tissues; P, primary tumor tissues; M, metastatic tumor tissues

vlnplot_tissue_type_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'TissueType', cols = c("seagreen", "pink", "salmon", "skyblue"), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = c("H", "P", "M", "U")))

generate_figs(vlnplot_tissue_type_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_tissue_type_exhaustion_markers', sep = ''), c(8,6))



# percentage of CARTEx detected

featplot_CARTEx_630_tissue_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'TissueType', cols=c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))
featplot_CARTEx_200_tissue_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'TissueType', cols=c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))
featplot_CARTEx_84_tissue_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'TissueType', cols=c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))

featplot_CARTEx_combined_tissue_type <- (featplot_CARTEx_630_tissue_type | featplot_CARTEx_200_tissue_type | featplot_CARTEx_84_tissue_type)
generate_figs(featplot_CARTEx_combined_tissue_type, paste('./plots/', experiment, '_featplot_CARTEx_combined_tissue_type', sep = ''), c(10,4))

generate_figs(featplot_CARTEx_200_tissue_type, paste('./plots/', experiment, '_featplot_CARTEx_200_tissue_type', sep = ''), c(1.5,2))



barplot_tissue_type_cancer_type <- BarPlotStackSplit(expt.obj, 'TissueType', 'CancerType', color_set = c("seagreen", "pink", "salmon", "skyblue")) + theme(legend.position = 'none', plot.title = element_blank())
generate_figs(barplot_tissue_type_cancer_type, paste('./plots/', experiment, '_prepare_barplot_tissue_type_cancer_type', sep = ''), c(10,2))

barplot_tissue_type_seurat_clusters <- BarPlotStackSplit(expt.obj, 'TissueType', 'seurat_clusters', color_set = c("seagreen", "pink", "salmon", "skyblue")) + theme(legend.position = 'none', plot.title = element_blank())
generate_figs(barplot_tissue_type_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_tissue_type_seurat_clusters', sep = ''), c(6,1.75))

# scale_x_discrete(labels = c('H', 'P', 'M', 'U'))







# examine differentiation


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('TissueType', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$TissueType <- factor(expt.obj.agg$TissueType, levels = c("Healthy donor", "Primary tumor tissue", "Metastatic tumor tissue", "Uninvolved normal tissue"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(expt.obj.agg$TissueType)
table(md$monaco)

# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, TissueType, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "TissueType", "pblabels"))


aggplot_CARTEx_200_TissueType_monaco_split_countsized <- md %>% ggplot(aes(x = TissueType, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab('CARTEx') + 
  scale_x_discrete(labels = c('H', 'P', 'M', 'U')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_TissueType_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_TissueType_monaco_split_countsized', sep = ''), c(4.5,3.5)) 










####################################################################################################
#################################### UMAP score split by metadata ##################################
####################################################################################################


# expt.obj@meta.data$monaco <- plyr::mapvalues(x = expt.obj@meta.data$monaco,
#                                          from = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'),
#                                          to = c('N', 'CM', 'EM', 'TE'))
# expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('N', 'CM', 'EM', 'TE'))


fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_monaco, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_monaco', sep = ''), c(8,2))

umap_activation_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_activation_monaco, paste('./plots/', experiment, '_prepare_umap_activation_monaco', sep = ''), c(8,2))

umap_anergy_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_anergy_monaco, paste('./plots/', experiment, '_prepare_umap_anergy_monaco', sep = ''), c(8,2))




umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))



table(expt.obj$TissueType)
table(expt.obj$CancerType)

# H, normal tissues from healthy donors; U, tumor-adjacent uninvolved tissues; P, primary tumor tissues; M, metastatic tumor tissues

umap_CARTEx_200_tissue_type <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'TissueType', split_ids = c('H', 'P', 'M', 'U'), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_tissue_type, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_tissue_type', sep = ''), c(8,2))



umap_CARTEx_200_seurat_clusters <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'seurat_clusters', split_ids = seq(0,10,1), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_seurat_clusters', sep = ''), c(22,2))



umap_TSR_tissue_type <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'TissueType', split_ids = c('H', 'P', 'M', 'U'), color_scale = fix.sc)
generate_figs(umap_TSR_tissue_type, paste('./plots/', experiment, '_prepare_umap_TSR_tissue_type', sep = ''), c(8,2))











# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 






####################################################################################################
######################################### Score Correlations #######################################
####################################################################################################



table_monaco_phase <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Phase)
corr_monaco_phase <- CorrespondenceAnalysisPlot(table_monaco_phase, "Monaco", "Phase")
generate_figs(corr_monaco_phase, paste('./plots/', experiment, '_prepare_corr_monaco_phase', sep = ''), c(4,4))

table_CARTEx_200_NKlike <- table(expt.obj@meta.data$CARTEx_200i, expt.obj@meta.data$NKlike_Texi)
corr_CARTEx_200_NKlike <- CorrespondenceAnalysisPlot(table_CARTEx_200_NKlike, "CARTEx_200i", "NKlike_Texi")
generate_figs(corr_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_corr_CARTEx_200_NKlike', sep = ''), c(4,4))






scatter_CARTEx_200_LCMV <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'LCMV_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'TissueType', cols =  c("seagreen", "pink", "salmon", "skyblue")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_LCMV, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'NKlike_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'TissueType', cols =  c("seagreen", "pink", "salmon", "skyblue")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'BBD_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'TissueType', cols =  c("seagreen", "pink", "salmon", "skyblue")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_BBD, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'PD1_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'TissueType', cols =  c("seagreen", "pink", "salmon", "skyblue")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_PD1, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))


# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "TissueType", size = 0.1, alpha = 1, palette = c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE,
                                         margin.params = list(fill = "TissueType", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "TissueType", size = 0.1, alpha = 1, palette = c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE,
                                           margin.params = list(fill = "TissueType", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "TissueType", size = 0.1, alpha = 1, palette = c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE,
                                        margin.params = list(fill = "TissueType", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "TissueType", size = 0.1, alpha = 1, palette = c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE,
                                        margin.params = list(fill = "TissueType", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))



scatter_CARTEx_200_TSR <- ggscatterhist(md, x = "CARTEx_200", y = "TSR", color = "TissueType", size = 0.1, alpha = 1, palette = c("seagreen", "pink", "salmon", "skyblue"), shuffle = TRUE,
                                        margin.params = list(fill = "TissueType", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_TSR), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_TSR', sep = ''), c(2,2))





####################################################################################################
############################################ Venn diagram ##########################################
####################################################################################################


library(ggvenn)


compare_sigs <- list(
  CARTEx = rownames(read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)),
  TSR = rownames(read.csv(paste(PATH_SIGNATURES, "Chu_2023_Nat_Med_T_stress_response.csv", sep = ''), header = FALSE, row.names = 1))
)

par(mar=c(7,5,1,1))
compare_sigs_CARTEx_TSR <- ggvenn(compare_sigs) 
generate_figs(compare_sigs_CARTEx_TSR, paste('./plots/', experiment, '_prepare_compare_sigs_CARTEx_TSR', sep = ''), c(3,3))









####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_TissueType <- DoHeatmap(expt.obj.downsample, group.by = "TissueType")
heatmap_allgenes_exCARTEx630_TissueType <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "TissueType")
heatmap_CARTEx630_TissueType <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "TissueType")
heatmap_CARTEx84_TissueType <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "TissueType")
heatmap_TSG_TissueType <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "TissueType")

generate_figs(heatmap_allgenes_TissueType, paste0('./plots/', experiment, '_explore_heatmap_allgenes_TissueType'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_TissueType, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_TissueType'), c(15, 12))
generate_figs(heatmap_CARTEx630_TissueType, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_TissueType'), c(15, 12))
generate_figs(heatmap_CARTEx84_TissueType, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_TissueType'), c(15, 12))
generate_figs(heatmap_TSG_TissueType, paste0('./plots/', experiment, '_explore_heatmap_TSG_TissueType'), c(15, 12))









