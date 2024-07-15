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
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_stim <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 1, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_stim, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_stim', sep = ''), c(10, 8))


# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Rested", ident.2 = "Control_Rested", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_rest <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 1, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_rest, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_rest', sep = ''), c(10, 8))



# Compare STAT3 GOF vs control
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF", ident.2 = "Control", group.by = "Affstat", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_control <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 1, 
                                                 selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_control, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_control', sep = ''), c(10, 8))







# percentage of CARTEx detected

featplot_CARTEx_630_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_200_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_84_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 30)) + ylim(c(-1, 8))

featplot_CARTEx_combined_affstatstim <- (featplot_CARTEx_630_affstatstim | featplot_CARTEx_200_affstatstim | featplot_CARTEx_84_affstatstim)
generate_figs(featplot_CARTEx_combined_affstatstim, paste('./plots/', experiment, '_featplot_CARTEx_combined_affstatstim', sep = ''), c(10,5))





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
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
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
  theme_bw() + theme(axis.title.x = element_blank()) + scale_x_discrete(labels = custom_labels)
generate_figs(aggplot_CARTEx_200_affstatstim_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_monaco_split', sep = ''), c(6,5)) 


# incorporate size
# md_count <- expt.obj@meta.data %>% group_by(monaco, AffstatStim, pblabels) %>% summarize(count = n(), .groups = 'drop')
# md_count$pblabels <- as.character(md_count$pblabels)

# replace affstatestim
# expt.obj@meta.data$AffstatStim

# md <- md %>% left_join(md_count, by = c("monaco", "AffstatStim", "pblabels"))

# aggplot_CARTEx_200_affstatstim_monaco_split_countsized <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = monaco, size = count)) +
  # geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  # scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  # theme_bw() + theme(axis.title.x = element_blank()) + scale_x_discrete(labels = custom_labels)
# generate_figs(aggplot_CARTEx_200_affstatstim_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_monaco_split_countsized', sep = ''), c(6,5)) 






