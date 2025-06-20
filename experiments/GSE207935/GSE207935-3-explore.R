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
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare STAT3 GOF (stim)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Stimulated", ident.2 = "Control_Stimulated", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

write.csv(rownames(subset(de_genes, p_val < 10e-6 & avg_log2FC > 0.5)), "./data/GSE207935_STAT3_GOF_stim.csv", row.names=FALSE)

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_stim <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_stim, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_stim', sep = ''), c(6, 5))



signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), 'CARTEx')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_stim_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_stim_200, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_stim_200', sep = ''), c(6, 5))









# Compare STAT3 GOF (rest)
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF_Rested", ident.2 = "Control_Rested", group.by = "AffstatStim", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

write.csv(rownames(subset(de_genes, p_val < 10e-6 & avg_log2FC > 0.5)), "./data/GSE207935_STAT3_GOF_rest.csv", row.names=FALSE)

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_rest <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_rest, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_rest', sep = ''), c(6, 5))


signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), 'CARTEx')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_rest_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_rest_200, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_rest_200', sep = ''), c(6, 5))








# Compare STAT3 GOF vs control
de_genes <- FindMarkers(expt.obj, ident.1 = "STAT3_GOF", ident.2 = "Control", group.by = "Affstat", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

write.csv(rownames(subset(de_genes, p_val < 10e-6 & avg_log2FC > 0.5)), "./data/GSE207935_STAT3_GOF_combined.csv", row.names=FALSE)

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_STAT3GOF_control <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                                 selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                 shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_STAT3GOF_control, paste('./plots/', experiment, '_plot_volcano_STAT3GOF_control', sep = ''), c(6, 5))



















# analyze STAT3 GOF signatures

STAT3_GOF_sig_stim <- rownames(read.csv("./data/GSE207935_STAT3_GOF_stim.csv", header = TRUE, row.names = 1))
STAT3_GOF_sig_rest <- rownames(read.csv("./data/GSE207935_STAT3_GOF_rest.csv", header = TRUE, row.names = 1))
STAT3_GOF_sig_combined <- rownames(read.csv("./data/GSE207935_STAT3_GOF_combined.csv", header = TRUE, row.names = 1))

STAT3_GOF_sig_intersect <- intersect(STAT3_GOF_sig_stim, STAT3_GOF_sig_rest)
write.csv(STAT3_GOF_sig_intersect, "./data/GSE207935_STAT3_GOF_intersect.csv", row.names=FALSE)


STAT3_GOF_sig_venn <- list(
  Stim = STAT3_GOF_sig_stim, 
  Rest = STAT3_GOF_sig_rest, 
  Combined = STAT3_GOF_sig_combined
)

library(ggvenn)
plot_STAT3_GOF_sig_venn <- ggvenn(STAT3_GOF_sig_venn, fill_color = c("violetred", "steelblue", "orange"), stroke_size = 0.5, set_name_size = 3, text_size = 3)
generate_figs(plot_STAT3_GOF_sig_venn, paste('./plots/', experiment, '_plot_STAT3_GOF_sig_venn', sep = ''), c(2,2))

####







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

featplot_CARTEx_630_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_200_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-1, 8))
featplot_CARTEx_84_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-1, 8))

featplot_CARTEx_combined_affstatstim <- (featplot_CARTEx_630_affstatstim | featplot_CARTEx_200_affstatstim | featplot_CARTEx_84_affstatstim)
generate_figs(featplot_CARTEx_combined_affstatstim, paste('./plots/', experiment, '_featplot_CARTEx_combined_affstatstim', sep = ''), c(10,4))


featplot_CARTEx_200_affstatstim <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'AffstatStim', cols=c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-1, 8))
generate_figs(featplot_CARTEx_200_affstatstim, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_affstatstim', sep = ''), c(1.5,2))




# aggregate for CARTEx scores (BAR CHART)


expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AffstatStim', 'pblabels'), return.seurat = TRUE)
expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg$AffstatStim <- factor(expt.obj.agg$AffstatStim, levels = c("Control-Rested", "Control-Stimulated", "STAT3-GOF-Rested", "STAT3-GOF-Stimulated"))

md <- expt.obj.agg@meta.data %>% as.data.table
glimpse(md)

custom_labels <- c('CTRL\n(R)', 'CTRL\n(S)', 'GOF\n(R)', 'GOF\n(S)')

agg_barplot_CARTEx_200 <- md %>% ggplot(aes(AffstatStim, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AffstatStim), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control-Rested','Control-Stimulated')), label = "p.signif", label.y = 0.6) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control-Rested','STAT3-GOF-Rested')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control-Stimulated','STAT3-GOF-Rested')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control-Stimulated','STAT3-GOF-Stimulated')), label = "p.signif", label.y = 1.8) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) + scale_x_discrete(labels = custom_labels)
generate_figs(agg_barplot_CARTEx_200, paste('./plots/', experiment, '_agg_barplot_CARTEx_200', sep = ''), c(3.2,3))










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

custom_labels <- c('CTRL\n(R)', 'CTRL\n(S)', 'GOF\n(R)', 'GOF\n(S)')

aggplot_CARTEx_200_affstatstim_monaco_split_countsized <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab('CARTEx') + 
  scale_x_discrete(labels = custom_labels) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_affstatstim_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_monaco_split_countsized', sep = ''), c(4.5,3)) 


aggplot_TSR_affstatstim_monaco_split_countsized <- md %>% ggplot(aes(x = AffstatStim, y = TSR, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab('Stress Response') + 
  scale_x_discrete(labels = custom_labels) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_TSR_affstatstim_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_TSR_affstatstim_monaco_split_countsized', sep = ''), c(4.5,3)) 










# aggregate without controls || exhaustion as linear differentiation

table(expt.obj$Tex_linear_diff)

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

table(expt.obj$AffstatStim, expt.obj$Tex_linear_diff)
table(expt.obj$AffstatStim, expt.obj$Tex_linear_diff, expt.obj$pblabels)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AffstatStim', 'Tex_linear_diff', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

# for some reason, "_" changed to "-"
expt.obj.agg$AffstatStim <- factor(expt.obj.agg$AffstatStim, levels = c("Control-Rested", "Control-Stimulated", "STAT3-GOF-Rested", "STAT3-GOF-Stimulated"))
expt.obj.agg$Tex_linear_diff <- factor(expt.obj.agg$Tex_linear_diff, levels = c("Tex-prog", "Tex-term", "Unclassified"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$Tex_linear_diff)
table(md$AffstatStim)
table(md$AffstatStim, md$Tex_linear_diff)


aggplot_CARTEx_200_affstatstim_Tex_linear_diff_split <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = Tex_linear_diff)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Tex-prog' = 'royalblue', 'Tex-term' = 'maroon', 'Unclassified' = 'lightgrey')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_affstatstim_Tex_linear_diff_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_affstatstim_Tex_linear_diff_split', sep = ''), c(6,5)) 



# incorporate size
md_count <- expt.obj@meta.data %>% group_by(Tex_linear_diff, AffstatStim, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md_count$AffstatStim <- factor(md_count$AffstatStim, levels = c("Control-Rested", "Control-Stimulated", "STAT3-GOF-Rested", "STAT3-GOF-Stimulated"))
md_count$Tex_linear_diff <- factor(md_count$Tex_linear_diff, levels = c("Tex-prog", "Tex-term", "Unclassified"))

# fix element naming scheme to match
md_count <- md_count %>% mutate(Tex_linear_diff = recode(Tex_linear_diff, "Tex_prog" = "Tex-prog", "Tex_term" = "Tex-term"))
md_count <- md_count %>% mutate(AffstatStim = recode(AffstatStim, "Control_Rested" = "Control-Rested", "Control_Stimulated" = "STAT3_GOF-Stimulated", "Control_Rested" = "STAT3-GOF-Rested", "STAT3_GOF_Stimulated" = "STAT3-GOF-Stimulated",))



md <- md %>% left_join(md_count, by = c("Tex_linear_diff", "AffstatStim", "pblabels"))
md$count <- md_count$count # some issue where the fix is to re-import
head(md)

table(md$Tex_linear_diff,md$AffstatStim)

custom_labels <- c('CTRL\n(R)', 'CTRL\n(S)', 'GOF\n(R)', 'GOF\n(S)')

aggplot_CARTEx_200_Tex_linear_diff_affstatstim_split_countsized <- md %>% ggplot(aes(x = AffstatStim, y = CARTEx_200, color = Tex_linear_diff, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2.5) +
  scale_color_manual(values = c('Tex-prog' = 'royalblue', 'Tex-term' = 'maroon', 'Unclassified' = 'lightgrey')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = custom_labels) + ylab('CARTEx') +
  scale_color_manual(labels=c("P", "T", "U"), values = c('royalblue', 'maroon', 'lightgrey'), name = 'Tex')
generate_figs(aggplot_CARTEx_200_Tex_linear_diff_affstatstim_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_Tex_linear_diff_affstatstim_split_countsized', sep = ''), c(4.5,3)) 


### 















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

umap_senescence_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Senescence"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_senescence_monaco, paste('./plots/', experiment, '_prepare_umap_senescence_monaco', sep = ''), c(8,2))

umap_TSR_monaco <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_TSR_monaco, paste('./plots/', experiment, '_prepare_umap_TSR_monaco', sep = ''), c(8,2))

umap_Tex_Term_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Tex_Term"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_Tex_Term_monaco, paste('./plots/', experiment, '_prepare_umap_Tex_Term_monaco', sep = ''), c(8,2))

umap_Tex_KLR_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Tex_KLR"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_Tex_KLR_monaco, paste('./plots/', experiment, '_prepare_umap_Tex_KLR_monaco', sep = ''), c(8,2))





umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))

umap_CARTEx_200_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_affstatstim, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_affstatstim', sep = ''), c(8,2))

umap_activation_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_activation_affstatstim, paste('./plots/', experiment, '_prepare_umap_activation_affstatstim', sep = ''), c(8,2))

umap_anergy_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_anergy_affstatstim, paste('./plots/', experiment, '_prepare_umap_anergy_affstatstim', sep = ''), c(8,2))

umap_senescence_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Senescence"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_senescence_affstatstim, paste('./plots/', experiment, '_prepare_umap_senescence_affstatstim', sep = ''), c(8,2))

umap_stemness_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Stemness"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_stemness_affstatstim, paste('./plots/', experiment, '_prepare_umap_stemness_affstatstim', sep = ''), c(8,2))

umap_TSR_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_TSR_affstatstim, paste('./plots/', experiment, '_prepare_umap_TSR_affstatstim', sep = ''), c(8,2))

umap_Tex_Term_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Tex_Term"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_Tex_Term_affstatstim, paste('./plots/', experiment, '_prepare_umap_Tex_Term_affstatstim', sep = ''), c(8,2))

umap_Tex_KLR_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Tex_KLR"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc)
generate_figs(umap_Tex_KLR_affstatstim, paste('./plots/', experiment, '_prepare_umap_Tex_KLR_affstatstim', sep = ''), c(8,2))




# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 

# dmaps

dmap_CARTEx_200_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_CARTEx_200_affstatstim[[i]] <- dmap_CARTEx_200_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_CARTEx_200_affstatstim, paste('./plots/', experiment, '_prepare_dmap_CARTEx_200_affstatstim', sep = ''), c(8,2))

dmap_activation_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_activation_affstatstim[[i]] <- dmap_activation_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_activation_affstatstim, paste('./plots/', experiment, '_prepare_dmap_activation_affstatstim', sep = ''), c(8,2))

dmap_anergy_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_anergy_affstatstim[[i]] <- dmap_anergy_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_anergy_affstatstim, paste('./plots/', experiment, '_prepare_dmap_anergy_affstatstim', sep = ''), c(8,2))

dmap_senescence_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Senescence"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_senescence_affstatstim[[i]] <- dmap_senescence_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_senescence_affstatstim, paste('./plots/', experiment, '_prepare_dmap_senescence_affstatstim', sep = ''), c(8,2))

dmap_TSR_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_TSR_affstatstim[[i]] <- dmap_TSR_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_TSR_affstatstim, paste('./plots/', experiment, '_prepare_dmap_TSR_affstatstim', sep = ''), c(8,2))

dmap_Tex_Term_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Tex_Term"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_Tex_Term_affstatstim[[i]] <- dmap_Tex_Term_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_Tex_Term_affstatstim, paste('./plots/', experiment, '_prepare_dmap_Tex_Term_affstatstim', sep = ''), c(8,2))

dmap_Tex_KLR_affstatstim <- FeaturePlotSplitBy(expt.obj, features = c("Tex_KLR"), split_identity = 'AffstatStim', split_ids = c("Ctrl (R)", "Ctrl (S)", "GOF (R)", "GOF (S)"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_Tex_KLR_affstatstim[[i]] <- dmap_Tex_KLR_affstatstim[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_Tex_KLR_affstatstim, paste('./plots/', experiment, '_prepare_dmap_Tex_KLR_affstatstim', sep = ''), c(8,2))

dmap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_CARTEx_200_monaco[[i]] <- dmap_CARTEx_200_monaco[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_CARTEx_200_monaco, paste('./plots/', experiment, '_prepare_dmap_CARTEx_200_monaco', sep = ''), c(8,2))

dmap_TSR_monaco <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_TSR_monaco[[i]] <- dmap_TSR_monaco[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_TSR_monaco, paste('./plots/', experiment, '_prepare_dmap_TSR_monaco', sep = ''), c(8,2))

dmap_Tex_Term_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Tex_Term"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_Tex_Term_monaco[[i]] <- dmap_Tex_Term_monaco[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_Tex_Term_monaco, paste('./plots/', experiment, '_prepare_dmap_Tex_Term_monaco', sep = ''), c(8,2))

dmap_Tex_KLR_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Tex_KLR"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc, reduction = 'dm')
for(i in 1:4){dmap_Tex_KLR_monaco[[i]] <- dmap_Tex_KLR_monaco[[i]] + xlim(c(-0.025, 0.025)) + ylim(c(-0.025, 0.025))}
generate_figs(dmap_Tex_KLR_monaco, paste('./plots/', experiment, '_prepare_dmap_Tex_KLR_monaco', sep = ''), c(8,2))








###
###




FeaturePlot(expt.obj, features = c("CARTEx_200", "TSR"), pt.size = 0.1, blend = TRUE)

FeaturePlot(expt.obj, features = c("CARTEx_200", "Activation"), pt.size = 0.1, blend = TRUE)

FeaturePlot(expt.obj, features = c("Activation", "TSR"), pt.size = 0.1, blend = TRUE)






####################################################################################################
######################################### Score Correlations #######################################
####################################################################################################



table_monaco_phase <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Phase)
corr_monaco_phase <- CorrespondenceAnalysisPlot(table_monaco_phase, "Monaco", "Phase")
generate_figs(corr_monaco_phase, paste('./plots/', experiment, '_prepare_corr_monaco_phase', sep = ''), c(4,4))

table_CARTEx_200_NKlike <- table(expt.obj@meta.data$CARTEx_200i, expt.obj@meta.data$NKlike_Texi)
corr_CARTEx_200_NKlike <- CorrespondenceAnalysisPlot(table_CARTEx_200_NKlike, "CARTEx_200i", "NKlike_Texi")
generate_figs(corr_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_corr_CARTEx_200_NKlike', sep = ''), c(4,4))






scatter_CARTEx_200_LCMV <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'LCMV_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'AffstatStim', cols =  c("lightsteelblue", "steelblue", "palevioletred", "violetred")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_LCMV, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'NKlike_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'AffstatStim', cols =  c("lightsteelblue", "steelblue", "palevioletred", "violetred")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'BBD_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'AffstatStim', cols =  c("lightsteelblue", "steelblue", "palevioletred", "violetred")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_BBD, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'PD1_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'AffstatStim', cols =  c("lightsteelblue", "steelblue", "palevioletred", "violetred")) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_PD1, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))


# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "AffstatStim", size = 0.1, alpha = 1, palette = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE,
                                         margin.params = list(fill = "AffstatStim", color = "black", size = 0.2), xlim = c(-2, 10), ylim = c(-2, 10), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "AffstatStim", size = 0.1, alpha = 1, palette = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE,
                                           margin.params = list(fill = "AffstatStim", color = "black", size = 0.2), xlim = c(-2, 10), ylim = c(-2, 10), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "AffstatStim", size = 0.1, alpha = 1, palette = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE,
                                        margin.params = list(fill = "AffstatStim", color = "black", size = 0.2), xlim = c(-2, 10), ylim = c(-2, 10), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "AffstatStim", size = 0.1, alpha = 1, palette = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE,
                                        margin.params = list(fill = "AffstatStim", color = "black", size = 0.2), xlim = c(-2, 10), ylim = c(-2, 10), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))

scatter_CARTEx_200_TSR <- ggscatterhist(md, x = "CARTEx_200", y = "TSR", color = "AffstatStim", size = 0.1, alpha = 1, palette = c("lightsteelblue", "steelblue", "palevioletred", "violetred"), shuffle = TRUE,
                                        margin.params = list(fill = "AffstatStim", color = "black", size = 0.2), xlim = c(-2, 10), ylim = c(-2, 10), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_TSR), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_TSR', sep = ''), c(2,2))

# table of cells with scores binned
table(expt.obj@meta.data$CARTEx_200i, expt.obj@meta.data$TSRi)

expt.obj@meta.data$custom_anno <- case_when(
  # (expt.obj@meta.data$Tex_Term > 1 | expt.obj@meta.data$Tex_KLR > 1) & expt.obj@meta.data$CARTEx_200 <= 1 & expt.obj@meta.data$TSR <= 1 ~ "exhausted (differentiation)",
  expt.obj@meta.data$Tex_Term > 1 & expt.obj@meta.data$CARTEx_200 <= 1 & expt.obj@meta.data$TSR <= 1 ~ "exhausted (differentiation)",
  expt.obj@meta.data$CARTEx_200 > 1 & expt.obj@meta.data$TSR <= 1 ~ "exhausted (CARTEx)",
  expt.obj@meta.data$TSR > 1 & expt.obj@meta.data$CARTEx_200 <= 1 ~ "stressed",
  TRUE ~ "normal"
  # TRUE ~ NA_character_
)


expt.obj@meta.data$custom_anno <- factor(expt.obj@meta.data$custom_anno, levels = c("normal", "exhausted (CARTEx)", "exhausted (differentiation)", "stressed"))


umap_custom_anno <- DimPlot(expt.obj, reduction = "umap", group.by = "custom_anno", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c("lightgrey", "indianred", "orchid", "cadetblue")) + 
  theme(plot.title = element_blank()) + scale_color_manual(labels=c("normal", "exhausted (CARTEx)", "exhausted (differentiation)", "stressed"), values = c("lightgrey", "indianred", "orchid", "cadetblue"))
generate_figs(umap_custom_anno, paste('./plots/', experiment, '_prepare_umap_custom_anno', sep = ''), c(4.3, 2))
# generate_figs(umap_custom_anno, paste('./plots/', experiment, '_prepare_umap_custom_anno', sep = ''), c(8, 4))

# check how the custom annotations bin against AffstatStim
table(expt.obj@meta.data$AffstatStim, expt.obj@meta.data$custom_anno)


custom_labels <- c('Ctrl\n(R)', 'Ctrl\n(S)', 'GOF\n(R)', 'GOF\n(S)')
barplot_affstatstim_custom_anno <- BarPlotStackSplit(expt.obj, 'custom_anno', 'AffstatStim', color_set = c("lightgrey", "indianred", "orchid", "cadetblue")) +
  scale_x_discrete(labels = custom_labels) + theme(legend.position = "none")
generate_figs(barplot_affstatstim_custom_anno, paste('./plots/', experiment, '_prepare_barplot_affstatstim_custom_anno', sep = ''), c(1.75,1.75))



###
###











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


















