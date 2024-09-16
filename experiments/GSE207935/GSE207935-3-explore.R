#!/usr/bin/env Rscript

### Script name: GSE207935-3-explore.R
### Description: explore seurat object for GSE207935
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE207935'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# NT5E corresponds to CD73
# ENTPD1 corresponds to CD39

VlnPlot(expt.obj, features=c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'Affstat')

vlnplot_affstat_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Affstat', cols = c("steelblue", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_affstat_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_affstat_exhaustion_markers', sep = ''), c(6,5))

vlnplot_stim_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Stim', cols = c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                             VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Stim', cols = c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                             VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Stim', cols = c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                             VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Stim', cols =c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                             VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Stim', cols = c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                             VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Stim', cols = c("lightslategrey", "darkslategrey"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_stim_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_stim_exhaustion_markers', sep = ''), c(6,5))


VlnPlot(expt.obj, features=c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)

custom_labels <- c('Control\nRested', 'Control\nStimulated', 'STAT3_GOF\nRested', 'STAT3_GOF\nStimulated')

vlnplot_affstatstim_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('LAG3'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels))



custom_labels <- c('Control (R)', 'Control (S)', 'STAT3_GOF (R)', 'STAT3_GOF (S)')

vlnplot_affstatstim_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('LAG3'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'AffstatStim', cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), y.max = 20)+theme(axis.title.x = element_blank()) + guides(fill=FALSE)+ scale_x_discrete(labels = custom_labels))


generate_figs(vlnplot_affstatstim_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_affstatstim_exhaustion_markers', sep = ''), c(6,5))



### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Stimulated", ident.2 = "Control_Stimulated", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_stim <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_stim, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_stim', sep = ''), c(6, 5))


# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Rested", ident.2 = "Control_Rested", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_rest <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_rest, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_rest', sep = ''), c(6, 5))



# Compare STAT3 GOF vs control
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF", ident.2 = "Control", group.by = "Affstat", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_control <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                                 selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                 shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_control, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_control', sep = ''), c(6, 5))



# examine CARTEx C2 genes
cartex_C2 <- rownames(read.csv(paste(PATH_CONSTRUCTION, "./data/cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))


# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Stimulated", ident.2 = "Control_Stimulated", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, 'C2')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_stim_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_stim_C2genes, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_stim_C2genes', sep = ''), c(6, 5))


# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Rested", ident.2 = "Control_Rested", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, 'C2')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_rest_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_rest_C2genes, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_rest_C2genes', sep = ''), c(6, 5))



# Compare STAT3 GOF vs control
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF", ident.2 = "Control", group.by = "Affstat", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, 'C2')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_control_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                                 selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                 shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_control_C2genes, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_control_C2genes', sep = ''), c(6, 5))








# percentage of CARTEx detected

featplot_CARTEx_630_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_200_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_84_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 30)) + ylim(c(-1, 8))

featplot_CARTEx_combined_affstatstim <- (featplot_CARTEx_630_affstatstim | featplot_CARTEx_200_affstatstim | featplot_CARTEx_84_affstatstim)
generate_figs(featplot_CARTEx_combined_affstatstim, paste('./plots/', experiment, '_featplot_CARTEx_combined_affstatstim', sep = ''), c(10,4))


generate_figs(featplot_CARTEx_200_affstatstim, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_affstatstim', sep = ''), c(2,4))





# examine differentiation

vlnplot_CARTEx_affstatstim_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AffstatStim', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_affstatstim_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_affstatstim_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_affstatstim_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AffstatStim', pt.size = 0, cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
                                              theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                              ylab("CARTEx 200") + ylim(-2,4),
                                            VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                              theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                              ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_affstatstim_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_affstatstim_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_affstatstim_monaco <- ggplot(md, aes(x = AffstatStim, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_affstatstim_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_affstatstim_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AffstatStim', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

# for some reason, "_" changed to "-"
expt.obj.agg$AffstatStim <- factor(expt.obj.agg$AffstatStim, levels = c("Control-Rested", "Control-Stimulated", "STAT3-GOF-Rested", "STAT3-GOF-Stimulated"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)
table(md$AffstatStim)

custom_labels <- c('Control (R)', 'Control (S)', 'STAT3_GOF (R)', 'STAT3_GOF (S)')

aggplot_CARTEx_200_affstatstim_monaco_split <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank()) + scale_x_discrete(labels = custom_labels)
generate_figs(aggplot_CARTEx_200_affstatstim_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_monaco_split', sep = ''), c(6,5)) 


# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, AffstatStim, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)

md <- md %>% left_join(md_count, by = c("monaco", "AffstatStim", "pblabels"))
md$count <- md_count$count # some issue where the fix is to re-import

aggplot_CARTEx_200_affstatstim_monaco_split_countsized <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank()) + scale_x_discrete(labels = custom_labels)
generate_figs(aggplot_CARTEx_200_affstatstim_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_monaco_split_countsized', sep = ''), c(6,5)) 




####################################################################################################
####################################### Examine gene expression ####################################
####################################################################################################


# heatmap

C5_genes <- rownames(read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1))
C2_genes <- rownames(read.csv(paste(PATH_CONSTRUCTION, "./data/cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))
fix.sc <- scale_fill_gradientn(colours = c("blue","lightgrey","red"), limits = c(-2,2))


expt.obj@meta.data$AffstatStim2 <- plyr::mapvalues(x = expt.obj@meta.data$AffstatStim,
                                                          from = c('Control_Rested', 'Control_Stimulated', 'STAT3_GOF_Rested', 'STAT3_GOF_Stimulated'),
                                                          to = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"))
expt.obj@meta.data$AffstatStim2 <- factor(expt.obj@meta.data$AffstatStim2, levels = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"))


expt.obj.subset <- subset(expt.obj, downsample = 100)
Idents(expt.obj.subset) <- "AffstatStim"

heatmap_all <- DoHeatmap(expt.obj.subset, group.by = 'AffstatStim2', group.colors = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), angle = 0, vjust = 1, hjust = 0.5) + fix.sc
heatmap_C5 <- DoHeatmap(expt.obj.subset, group.by = 'AffstatStim2', group.colors = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), angle = 0, vjust = 1, hjust = 0.5, features = C5_genes) + fix.sc + ylab("C5 genes") + theme(axis.text.y = element_blank())
heatmap_C2 <- DoHeatmap(expt.obj.subset, group.by = 'AffstatStim2', group.colors = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), angle = 0, vjust = 1, hjust = 0.5, features = C2_genes) + fix.sc + ylab("C2 genes") + theme(axis.text.y = element_blank())

heatmap_C5 + heatmap_C2




expt.obj.agg <- AggregateExpression(expt.obj, group.by = 'AffstatStim2', return.seurat = TRUE)
# expt.obj.agg@assays$RNA$scale.data

# exhaustion genes
pb_heatmap_C5 <- DoHeatmap(expt.obj.agg, group.by = 'AffstatStim2', group.colors = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), draw.lines = FALSE, features = C5_genes, angle = 0, vjust = 10, hjust = 0.5) + fix.sc + labs(title = 'C5') + ylab("C5 genes") + theme(axis.text.y = element_blank(), legend.position="none", text = element_text(size = 16))

# fitness genes
pb_heatmap_C2 <- DoHeatmap(expt.obj.agg, group.by = 'AffstatStim2', group.colors = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), draw.lines = FALSE, features = C2_genes, angle = 0, vjust = 10, hjust = 0.5) + fix.sc + labs(title = 'C2') + ylab("C2 genes") + theme(axis.text.y = element_blank(), legend.position="none", text = element_text(size = 16))

# pb_heatmap_C5 + pb_heatmap_C2
pb_heatmap_combined <- ggpubr::ggarrange(pb_heatmap_C5, pb_heatmap_C2,
                                 common.legend = TRUE, legend = "right", align = "hv", ncol = 2)
generate_figs(pb_heatmap_combined, paste('./plots/', experiment, '_aggplot_heatmap_C5_C2', sep = ''), c(8,5)) 




gene_expression_long <- melt(expt.obj.agg@assays$RNA$scale.data)

pb_vlnplot_all <- ggplot(gene_expression_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_violin(trim = TRUE) + scale_fill_manual(values = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black") +
  # geom_boxplot(width=0.1) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','Ctrl (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOF (R)','GOF (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','GOF (R)')), label = "p.signif", label.y = 2.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (S)','GOF (S)')), label = "p.signif", label.y = 3.3) + 
  theme_classic() + labs(y = "Mean Expression Level") + ylim(-2,3.7) + theme(legend.position = "none", axis.title.x = element_blank()) + labs(title="All")

pb_vlnplot_C5 <- ggplot(subset(gene_expression_long, Var1 %in% C5_genes), aes(x = Var2, y = value, fill = Var2)) + 
  geom_violin(trim = TRUE) + scale_fill_manual(values = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black") +
  # geom_boxplot(width=0.1) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','Ctrl (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOF (R)','GOF (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','GOF (R)')), label = "p.signif", label.y = 2.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (S)','GOF (S)')), label = "p.signif", label.y = 3.3) + 
  theme_classic() + labs(y = "Mean Expression Level") + ylim(-2,3.7) + theme(legend.position = "none", axis.title.x = element_blank()) + labs(title="C5")

pb_vlnplot_C2 <- ggplot(subset(gene_expression_long, Var1 %in% C2_genes), aes(x = Var2, y = value, fill = Var2)) + 
  geom_violin(trim = TRUE) + scale_fill_manual(values = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black") +
  # geom_boxplot(width=0.1) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','Ctrl (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOF (R)','GOF (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','GOF (R)')), label = "p.signif", label.y = 2.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (S)','GOF (S)')), label = "p.signif", label.y = 3.3) + 
  theme_classic() + labs(y = "Mean Expression Level") + ylim(-2,3.7) + theme(legend.position = "none", axis.title.x = element_blank()) + labs(title="C2")

pb_vlnplot_C5 + pb_vlnplot_C2
generate_figs(pb_vlnplot_C5 + pb_vlnplot_C2, paste('./plots/', experiment, '_aggplot_vlnplot_C5_C2', sep = ''), c(5,2.5)) 


# geom_quasirandom(groupOnX = FALSE, size = 0.1) +


cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

pb_vlnplot_CARTEx_200 <- ggplot(subset(gene_expression_long, Var1 %in% rownames(cartex_200_weights)), aes(x = Var2, y = value, fill = Var2)) + 
  geom_violin(trim = TRUE) + scale_fill_manual(values = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black") +
  # geom_boxplot(width=0.1) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','Ctrl (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOF (R)','GOF (S)')), label = "p.signif", label.y = 1.7) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (R)','GOF (R)')), label = "p.signif", label.y = 2.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Ctrl (S)','GOF (S)')), label = "p.signif", label.y = 3.3) + 
  theme_classic() + labs(y = "Mean Expression Level") + ylim(-2,3.7) + theme(legend.position = "none", axis.title.x = element_blank()) + labs(title="CARTEx 200")

pb_vlnplot_all + pb_vlnplot_CARTEx_200
generate_figs(pb_vlnplot_all + pb_vlnplot_CARTEx_200, paste('./plots/', experiment, '_aggplot_vlnplot_all_CARTEx_200', sep = ''), c(5,2.5)) 


generate_figs(pb_vlnplot_all + pb_vlnplot_C5 + pb_vlnplot_C2, paste('./plots/', experiment, '_aggplot_vlnplot_all_C5_C2', sep = ''), c(8,2.5)) 





fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))

umap_CARTEx_200_splitby_affstatstim <- FeaturePlot(expt.obj, features = c("CARTEx_200"), order = FALSE, pt.size = 0.1, split.by = 'AffstatStim2') + fix.sc
umap_activation_splitby_affstatstim <- FeaturePlot(expt.obj, features = c("Activation"), order = FALSE, pt.size = 0.1, split.by = 'AffstatStim2') + fix.sc
umap_anergy_splitby_affstatstim <- FeaturePlot(expt.obj, features = c("Anergy"), order = FALSE, pt.size = 0.1, split.by = 'AffstatStim2') + fix.sc
umap_senescence_splitby_affstatstim <- FeaturePlot(expt.obj, features = c("Senescence"), order = FALSE, pt.size = 0.1, split.by = 'AffstatStim2') + fix.sc
umap_stemness_splitby_affstatstim <- FeaturePlot(expt.obj, features = c("Stemness"), order = FALSE, pt.size = 0.1, split.by = 'AffstatStim2') + fix.sc

generate_figs(umap_CARTEx_200_splitby_affstatstim, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_splitby_affstatstim', sep = ''), c(9.5, 2.5))
generate_figs(umap_activation_splitby_affstatstim, paste('./plots/', experiment, '_prepare_umap_activation_splitby_affstatstim', sep = ''), c(9.5, 2.5))
generate_figs(umap_anergy_splitby_affstatstim, paste('./plots/', experiment, '_prepare_umap_anergy_splitby_affstatstim', sep = ''), c(9.5, 2.5))
generate_figs(umap_senescence_splitby_affstatstim, paste('./plots/', experiment, '_prepare_umap_senescence_splitby_affstatstim', sep = ''), c(9.5, 2.5))
generate_figs(umap_stemness_splitby_affstatstim, paste('./plots/', experiment, '_prepare_umap_stemness_splitby_affstatstim', sep = ''), c(9.5, 2.5))















