#!/usr/bin/env Rscript

### Script name: GSE136874-3-explore.R
### Description: explore seurat object for GSE136874
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39


vlnplot_CAR_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'CAR', ncol = 3, cols = c('dodgerblue', 'indianred'), y.max = 3)

vlnplot_CAR_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('LAG3'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))


generate_figs(vlnplot_CAR_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_CAR_exhaustion_markers', sep = ''), c(5,3.5))


VlnPlot(expt.obj, features = c('SLAMF6', 'CD69'), group.by = 'monaco', ncol = 2, y.max = 3, cols = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3'))


# percentage of CARTEx detected

featplot_CARTEx_630_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 35)) + ylim(c(-3, 5))
featplot_CARTEx_200_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 35)) + ylim(c(-3, 5))
featplot_CARTEx_84_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 35)) + ylim(c(-3, 5))

featplot_CARTEx_combined_CAR <- (featplot_CARTEx_630_CAR | featplot_CARTEx_200_CAR | featplot_CARTEx_84_CAR)
generate_figs(featplot_CARTEx_combined_CAR, paste('./plots/', experiment, '_featplot_CARTEx_combined_CAR', sep = ''), c(10,4))


featplot_CARTEx_200_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected') + xlim(c(0, 35)) + ylim(c(-3, 5))
generate_figs(featplot_CARTEx_200_CAR, paste('./plots/', experiment, '_featplot_CARTEx_200_CAR', sep = ''), c(1.5,2))


# examine differentiation

vlnplot_CARTEx_CAR_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'CAR', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_CAR_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_CAR_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_CAR_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'CAR', pt.size = 0, cols = c('dodgerblue', 'indianred')) +
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4),
                                         VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_CAR_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_CAR_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_CAR_monaco <- ggplot(md, aes(x = CAR, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_CAR_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_CAR_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('CAR', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$CAR <- factor(expt.obj.agg$CAR, levels = c("CD19", "GD2"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)

aggplot_CARTEx_200_CAR_monaco_split <- md %>% ggplot(aes(x = CAR, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_CAR_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_CAR_monaco_split', sep = ''), c(6,5)) 



# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, CAR, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "CAR", "pblabels"))


aggplot_CARTEx_200_CAR_monaco_split_countsized <- md %>% ggplot(aes(x = CAR, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c('CD19', 'GD2')) + ylab('CARTEx') +
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_CAR_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_CAR_monaco_split_countsized', sep = ''), c(3,3)) 




# code for volcano plots located in GSE136874-6-transition.R



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

umap_CARTEx_200_CAR <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'CAR', split_ids = c("CD19", "GD2"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_CAR, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_CAR', sep = ''), c(4,2))

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






scatter_CARTEx_200_LCMV <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'LCMV_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'CAR', cols =  c('dodgerblue', 'indianred')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_LCMV, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'NKlike_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'CAR', cols =  c('dodgerblue', 'indianred')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'BBD_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'CAR', cols =  c('dodgerblue', 'indianred')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_BBD, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'PD1_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'CAR', cols =  c('dodgerblue', 'indianred')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_PD1, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))


# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "CAR", size = 0.1, alpha = 1, palette = c('dodgerblue', 'indianred'), shuffle = TRUE,
                                         margin.params = list(fill = "CAR", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "CAR", size = 0.1, alpha = 1, palette = c('dodgerblue', 'indianred'), shuffle = TRUE,
                                           margin.params = list(fill = "CAR", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "CAR", size = 0.1, alpha = 1, palette = c('dodgerblue', 'indianred'), shuffle = TRUE,
                                        margin.params = list(fill = "CAR", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "CAR", size = 0.1, alpha = 1, palette = c('dodgerblue', 'indianred'), shuffle = TRUE,
                                        margin.params = list(fill = "CAR", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))














####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

# expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_CAR <- DoHeatmap(expt.obj, group.by = "CAR")
heatmap_allgenes_exCARTEx630_CAR <- DoHeatmap(expt.obj, features = allgenes_exCARTEx, group.by = "CAR")
heatmap_CARTEx630_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_630_weights), group.by = "CAR")
heatmap_CARTEx200_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_200_weights), group.by = "CAR")
heatmap_CARTEx84_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_84_weights), group.by = "CAR")
heatmap_TSG_CAR <- DoHeatmap(expt.obj, features = TSG, group.by = "CAR")

generate_figs(heatmap_allgenes_CAR, paste0('./plots/', experiment, '_explore_heatmap_allgenes_CAR'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_CAR, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx630_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx200_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx84_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_CAR'), c(15, 12))
generate_figs(heatmap_TSG_CAR, paste0('./plots/', experiment, '_explore_heatmap_TSG_CAR'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'CAR')
# FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))
VlnPlot(expt.obj, features = c('HNF1B'), group.by = 'CAR')



