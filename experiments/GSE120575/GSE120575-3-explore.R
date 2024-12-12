#!/usr/bin/env Rscript

### Script name: GSE120575-3-explore.R
### Description: explore seurat object for GSE120575
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE120575'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39


vlnplot_response_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'NT5E', 'ENTPD1'), group.by = 'characteristics_response', ncol = 3, cols = c("firebrick", "seagreen"), y.max = 3)


vlnplot_response_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_response_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_response_exhaustion_markers', sep = ''), c(8,6))


# triple checkpoint expression
triple_checkpoint_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")[c("PDCD1", "HAVCR2", "LAG3"),]
expt.obj@meta.data$triple_checkpoint_expression <- diff(triple_checkpoint_expression@p)

scatterplot_triple_checkpoint_expression_cols <- colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$triple_checkpoint_expression)))
scatterplot_CARTEx_200_activation_triple_checkpoint_expression <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'Activation', group.by = "triple_checkpoint_expression", cols = scatterplot_triple_checkpoint_expression_cols, shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_200_activation_triple_checkpoint_expression, paste('./plots/', experiment, '_prepare_scatterplot_CARTEx_200_activation_triple_checkpoint_expression', sep = ''), c(8,6))


# correspondence analysis from contingency tables

table_triple_checkpoint_expression_phase <- table(expt.obj@meta.data$triple_checkpoint_expression, expt.obj@meta.data$Phase)
corr_triple_checkpoint_expression_phase <- CorrespondenceAnalysisPlot(table_triple_checkpoint_expression_phase, "Triple Checkpoint Expression", "Phase")
generate_figs(corr_triple_checkpoint_expression_phase, paste('./plots/', experiment, '_prepare_corr_triple_checkpoint_expression_phase', sep = ''), c(4,4))

table_timepoint_response_phase <- table(expt.obj@meta.data$Timepoint_Response, expt.obj@meta.data$Phase)
corr_timepoint_response_phase <- CorrespondenceAnalysisPlot(table_timepoint_response_phase, "Timepoint-Response", "Phase")
generate_figs(corr_timepoint_response_phase, paste('./plots/', experiment, '_prepare_corr_timepoint_response_phase', sep = ''), c(4,4))

table_monaco_phase <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Phase)
corr_monaco_phase <- CorrespondenceAnalysisPlot(table_monaco_phase, "Monaco", "Phase")
generate_figs(corr_monaco_phase, paste('./plots/', experiment, '_prepare_corr_monaco_phase', sep = ''), c(4,4))

table_timepoint_response_monaco <- table(expt.obj@meta.data$Timepoint_Response, expt.obj@meta.data$monaco)
corr_timepoint_response_monaco <- CorrespondenceAnalysisPlot(table_timepoint_response_monaco, "Timepoint-Response", "Monaco")
generate_figs(corr_timepoint_response_monaco, paste('./plots/', experiment, '_prepare_corr_timepoint_response_monaco', sep = ''), c(4,4))



table_triple_checkpoint_expression_CARTEx_200i <- table(expt.obj@meta.data$triple_checkpoint_expression, expt.obj@meta.data$CARTEx_200i)
corr_triple_checkpoint_expression_CARTEx_200i <- CorrespondenceAnalysisPlot(table_triple_checkpoint_expression_CARTEx_200i, "Triple Checkpoint Expression", "CARTEx 200")
generate_figs(corr_triple_checkpoint_expression_CARTEx_200i, paste('./plots/', experiment, '_prepare_corr_triple_checkpoint_expression_CARTEx_200i', sep = ''), c(4,4))

table_NKlike_Texi_CARTEx_200i <- table(expt.obj@meta.data$NKlike_Texi, expt.obj@meta.data$CARTEx_200i)
corr_NKlike_Texi_CARTEx_200i <- CorrespondenceAnalysisPlot(table_NKlike_Texi_CARTEx_200i, "NK-like", "CARTEx 200")
generate_figs(corr_NKlike_Texi_CARTEx_200i, paste('./plots/', experiment, '_prepare_corr_NKlike_Texi_CARTEx_200i', sep = ''), c(4,4))

table_activationi_CARTEx_200i <- table(expt.obj@meta.data$Activationi, expt.obj@meta.data$CARTEx_200i)
table_activationi_CARTEx_200i <- CleanTable(table_activationi_CARTEx_200i)
corr_activationi_CARTEx_200i <- CorrespondenceAnalysisPlot(table_activationi_CARTEx_200i, "Activation", "CARTEx 200")
generate_figs(corr_activationi_CARTEx_200i, paste('./plots/', experiment, '_prepare_corr_activationi_CARTEx_200i', sep = ''), c(4,4))

table_anergyi_CARTEx_200i <- table(expt.obj@meta.data$Anergyi, expt.obj@meta.data$CARTEx_200i)
table_anergyi_CARTEx_200i <- CleanTable(table_anergyi_CARTEx_200i)
corr_anergyi_CARTEx_200i <- CorrespondenceAnalysisPlot(table_anergyi_CARTEx_200i, "Anergy", "CARTEx 200")
generate_figs(corr_anergyi_CARTEx_200i, paste('./plots/', experiment, '_prepare_corr_anergyi_CARTEx_200i', sep = ''), c(4,4))

table_anergyi_NKlike_Texi <- table(expt.obj@meta.data$Anergyi, expt.obj@meta.data$NKlike_Texi)
table_anergyi_NKlike_Texi <- CleanTable(table_anergyi_NKlike_Texi)
corr_anergyi_NKlike_Texi <- CorrespondenceAnalysisPlot(table_anergyi_NKlike_Texi, "Anergy", "NK-like")
generate_figs(corr_anergyi_NKlike_Texi, paste('./plots/', experiment, '_prepare_corr_anergyi_NKlike_Texi', sep = ''), c(4,4))




# COLORFUL gradient

meta_trunc <- data.frame(expt.obj@meta.data$CARTEx_200, expt.obj@meta.data$characteristics_response, expt.obj@meta.data$characteristics_therapy, expt.obj@meta.data$Timepoint)
colnames(meta_trunc) <- c('CARTEx_score', 'Response', 'Therapy', 'Timepoint')

gradient_timepoint_response <- ggplot(meta_trunc, aes(x=factor(Timepoint, level = c("Pre","Post")),y=CARTEx_score)) + 
  geom_boxplot() + geom_jitter(aes(color = Response), alpha = 0.5) + xlab("Therapy") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre','Post')), label = "p.signif") +
  scale_color_manual(values = c('NR' = 'firebrick', 'R' = 'seagreen')) +
  ylab('CARTEx 200')

generate_figs(gradient_timepoint_response, paste('./plots/', experiment, '_gradient_timepoint_response', sep = ''), c(4,6))


# percentage of CARTEx detected

featplot_CARTEx_630_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_200_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_84_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))

featplot_CARTEx_combined_response <- (featplot_CARTEx_630_response | featplot_CARTEx_200_response | featplot_CARTEx_84_response)
generate_figs(featplot_CARTEx_combined_response, paste('./plots/', experiment, '_featplot_CARTEx_combined_response', sep = ''), c(10,4))

featplot_CARTEx_200_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))
generate_figs(featplot_CARTEx_200_response, paste('./plots/', experiment, '_featplot_CARTEx_200_response', sep = ''), c(1.5,2))




### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare NR vs R
de_genes <- FindMarkers(expt.obj, ident.1 = "NR", ident.2 = "R", group.by = "characteristics_response", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_response <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                  pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                  selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                  shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                  xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_response, paste('./plots/', experiment, '_plot_volcano_response', sep = ''), c(6, 5))





signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), "CARTEx")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_response_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_response_200, paste('./plots/', experiment, '_plot_volcano_response_200', sep = ''), c(6, 5))





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

umap_CARTEx_200_response <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'characteristics_response', split_ids = c("NR", "R"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_response, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_response', sep = ''), c(4,2))

# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 




### examining baseline (exclude post)

# expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
Idents(expt.obj) <- "Timepoint_Response"
expt.obj <- subset(expt.obj, idents = c("Post-NR", "Post-R"), invert = TRUE)



umap_response_baseline <- DimPlot(expt.obj, reduction = "umap", group.by = "characteristics_response", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c("firebrick", "seagreen")) + theme(plot.title = element_blank())
generate_figs(umap_response_baseline, paste('./plots/', experiment, '_prepare_umap_response_baseline', sep = ''), c(3,2))


fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_200_baseline <- FeaturePlot(expt.obj, features = c("CARTEx_200"), order = TRUE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
generate_figs(umap_CARTEx_200_baseline, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_baseline', sep = ''), c(2.8,2))


barplot_monaco_timepoint_response_baseline <- BarPlotStackSplit(expt.obj, 'monaco', 'Timepoint_Response', color_set = c('deepskyblue','seagreen','darkgoldenrod','plum3'))
generate_figs(barplot_monaco_timepoint_response_baseline, paste('./plots/', experiment, '_prepare_barplot_monaco_timepoint_response_baseline', sep = ''), c(8,4))

# https://github.com/satijalab/seurat/issues/3366#issuecomment-674262907

vlnplot_response_exhaustion_markers_baseline <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'characteristics_response', ncol = 3, cols = c("firebrick", "seagreen"), y.max = 3)

vlnplot_response_exhaustion_markers_baseline <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                          VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                          VlnPlot(expt.obj, features=c('LAG3'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                          VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                          VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                          VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'characteristics_response', cols = c("firebrick", "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_response_exhaustion_markers_baseline, paste('./plots/', experiment, '_prepare_vlnplot_response_exhaustion_markers_baseline', sep = ''), c(8,6))



### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare Pre-NR vs Pre-R
de_genes <- FindMarkers(expt.obj, ident.1 = "Pre-NR", ident.2 = "Pre-R", group.by = "Timepoint_Response", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_baseline_response <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                  pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                  selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                  shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                  xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_baseline_response, paste('./plots/', experiment, '_plot_volcano_baseline_response', sep = ''), c(6, 5))



signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), "CARTEx")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_baseline_response_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                  pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                  selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                  shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                  xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_baseline_response_200, paste('./plots/', experiment, '_plot_volcano_baseline_response_200', sep = ''), c(6, 5))








# examine CARTEx C2 genes
cartex_C2 <- rownames(read.csv(paste(PATH_WEIGHTS, "cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))


# Compare Pre-NR vs Pre-R
de_genes <- FindMarkers(expt.obj, ident.1 = "Pre-NR", ident.2 = "Pre-R", group.by = "Timepoint_Response", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, "C2")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_baseline_response_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                  pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                  selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                  shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                  xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_baseline_response_C2genes, paste('./plots/', experiment, '_plot_volcano_baseline_response_C2genes', sep = ''), c(6, 5))




# percentage of CARTEx detected

featplot_CARTEx_630_baseline_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_200_baseline_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_84_baseline_response <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'characteristics_response', cols=c("firebrick", "seagreen"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 20)) + ylim(c(-3, 5))

featplot_CARTEx_combined_baseline_response <- (featplot_CARTEx_630_baseline_response | featplot_CARTEx_200_baseline_response | featplot_CARTEx_84_baseline_response)
generate_figs(featplot_CARTEx_combined_baseline_response, paste('./plots/', experiment, '_featplot_CARTEx_combined_baseline_response', sep = ''), c(10,4))

generate_figs(featplot_CARTEx_200_baseline_response, paste('./plots/', experiment, '_featplot_CARTEx_200_baseline_response', sep = ''), c(2,4))



# examine differentiation

vlnplot_CARTEx_response_monaco_split_baseline_response <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'characteristics_response', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_response_monaco_split_baseline_response, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_response_monaco_split_baseline_response', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_response_monaco_baseline_response <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'characteristics_response', pt.size = 0, cols = c('firebrick', 'seagreen')) +
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4),
                                         VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_response_monaco_baseline_response, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_response_monaco_baseline_response', sep = ''), c(12,5)) 


md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_response_monaco_baseline_response <- ggplot(md, aes(x = characteristics_response, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_response_monaco_baseline_response, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_response_monaco_baseline_response', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('characteristics_response', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$characteristics_response <- factor(expt.obj.agg$characteristics_response, levels = c("NR", "R"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)

aggplot_CARTEx_200_response_monaco_split_baseline_response <- md %>% ggplot(aes(x = characteristics_response, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_response_monaco_split_baseline_response, paste('./plots/', experiment, '_aggplot_CARTEx_200_response_monaco_split_baseline_response', sep = ''), c(6,5)) 

# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, characteristics_response, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "characteristics_response", "pblabels"))

aggplot_CARTEx_200_response_monaco_split_baseline_response_countsized <- md %>% ggplot(aes(x = characteristics_response, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab('CARTEx') +
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_response_monaco_split_baseline_response_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_response_monaco_split_baseline_response_countsized', sep = ''), c(3,3)) 





####################################################################################################
######################################### Score Correlations #######################################
####################################################################################################



table_monaco_phase <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Phase)
corr_monaco_phase <- CorrespondenceAnalysisPlot(table_monaco_phase, "Monaco", "Phase")
generate_figs(corr_monaco_phase, paste('./plots/', experiment, '_prepare_corr_monaco_phase_baseline_response', sep = ''), c(4,4))

table_CARTEx_200_NKlike <- table(expt.obj@meta.data$CARTEx_200i, expt.obj@meta.data$NKlike_Texi)
corr_CARTEx_200_NKlike <- CorrespondenceAnalysisPlot(table_CARTEx_200_NKlike, "CARTEx_200i", "NKlike_Texi")
generate_figs(corr_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_corr_CARTEx_200_NKlike_baseline_response', sep = ''), c(4,4))






scatter_CARTEx_200_LCMV <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'LCMV_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'characteristics_response', cols =  c("firebrick", "seagreen")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_LCMV, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'NKlike_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'characteristics_response', cols =  c("firebrick", "seagreen")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'BBD_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'characteristics_response', cols =  c("firebrick", "seagreen")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_BBD, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'PD1_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'characteristics_response', cols =  c("firebrick", "seagreen")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_PD1, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1_baseline_response', sep = ''), c(2,2))


# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "characteristics_response", size = 0.1, alpha = 1, palette = c("firebrick", "seagreen"), shuffle = TRUE,
                                         margin.params = list(fill = "characteristics_response", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "characteristics_response", size = 0.1, alpha = 1, palette = c("firebrick", "seagreen"), shuffle = TRUE,
                                           margin.params = list(fill = "characteristics_response", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "characteristics_response", size = 0.1, alpha = 1, palette = c("firebrick", "seagreen"), shuffle = TRUE,
                                        margin.params = list(fill = "characteristics_response", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD_baseline_response', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "characteristics_response", size = 0.1, alpha = 1, palette = c("firebrick", "seagreen"), shuffle = TRUE,
                                        margin.params = list(fill = "characteristics_response", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1_baseline_response', sep = ''), c(2,2))












