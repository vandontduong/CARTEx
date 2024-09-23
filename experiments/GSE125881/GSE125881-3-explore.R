#!/usr/bin/env Rscript

### Script name: GSE125881-3-explore.R
### Description: explore seurat object for GSE125881
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
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

vlnplot_group_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'Group', ncol = 3, cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)

vlnplot_group_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_group_exhaustion_markers, paste('./plots/', experiment, '_explore_vlnplot_group_exhaustion_markers', sep = ''), c(8,6))

vlnplot_group2_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Group2', cols = colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_group2_exhaustion_markers, paste('./plots/', experiment, '_explore_vlnplot_group2_exhaustion_markers', sep = ''), c(8,6))


# percentage of CARTEx detected

featplot_CARTEx_630_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 30)) + ylim(c(-4, 6))
featplot_CARTEx_200_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 30)) + ylim(c(-4, 6))
featplot_CARTEx_84_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 30)) + ylim(c(-4, 6))

featplot_CARTEx_combined_group <- (featplot_CARTEx_630_group | featplot_CARTEx_200_group | featplot_CARTEx_84_group)
generate_figs(featplot_CARTEx_combined_group, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_group', sep = ''), c(10,4))

featplot_CARTEx_630_group2 <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'Group2', cols=colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-4, 6))
featplot_CARTEx_200_group2 <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Group2', cols=colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-4, 6))
featplot_CARTEx_84_group2 <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'Group2', cols=colorRampPalette(c("cadetblue", "violet", "darkorchid"))(length(unique(expt.obj@meta.data$Group2))), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-4, 6))

featplot_CARTEx_combined_group2 <- (featplot_CARTEx_630_group2 | featplot_CARTEx_200_group2 | featplot_CARTEx_84_group2)
generate_figs(featplot_CARTEx_combined_group2, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_group2', sep = ''), c(10,4))

generate_figs(featplot_CARTEx_200_group2, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_group2', sep = ''), c(1.5,2))





# examine differentiation

vlnplot_CARTEx_group_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_group_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_group_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group', pt.size = 0, cols = c('red', 'violetred', 'violet', 'purple')) +
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4),
                                             VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco', sep = ''), c(12,5)) 



vlnplot_CARTEx_group2_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group2', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_group2_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group2_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_group2_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group2', pt.size = 0, cols = c('red', 'violetred', 'violet', 'purple')) +
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4),
                                         VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_group2_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group2_monaco', sep = ''), c(12,5)) 





md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_group_monaco <- ggplot(md, aes(x = Group, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_group_monaco', sep = ''), c(6,5)) 

swarmplot_CARTEx_group2_monaco <- ggplot(md, aes(x = Group2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group2_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_group2_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('Group2', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$Group2 <- factor(expt.obj.agg$Group2, levels = c("IP", "Early", "Late", "Very Late"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_group2_monaco_split <- md %>% ggplot(aes(x = Group2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_group2_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_group2_monaco_split', sep = ''), c(6,5)) 



# incorporate size
table(md$pblabels, md$monaco)
agg.ref <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Group2, expt.obj@meta.data$pblabels)
prop.table(agg.ref)

md_count <- expt.obj@meta.data %>% group_by(monaco, Group2, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "Group2", "pblabels"))



aggplot_CARTEx_200_group2_monaco_split_countsized <- md %>% ggplot(aes(x = Group2, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c('IP', 'E', 'L', 'VL')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_group2_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_group2_monaco_split_countsized', sep = ''), c(4,3)) 




####################################################################################################
#################################### Differential gene expression ##################################
####################################################################################################

### label CARTEx in exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))


# Compare ExpansionPeak vs Late
de_genes <- FindMarkers(expt.obj, ident.1 = "ExpansionPeak", ident.2 = "Late", group.by = "Group", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_ExpansionPeak_Late <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                   pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                   selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                   shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                   xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_ExpansionPeak_Late, paste('./plots/', experiment, '_explore_plot_volcano_ExpansionPeak_Late', sep = ''), c(6, 5))


# Compare ExpansionPeak vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = "ExpansionPeak", ident.2 = "IP", group.by = "Group", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_ExpansionPeak_IP <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                 selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                 shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_ExpansionPeak_IP, paste('./plots/', experiment, '_explore_plot_volcano_ExpansionPeak_IP', sep = ''), c(6, 5))

VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Patient', split.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple')) + theme(legend.position = 'none', axis.text.x = element_text(angle = 0, hjust = 0.5)) + ylab('CARTEx 200') + xlab(NULL)



table(expt.obj@meta.data$Group, expt.obj@meta.data$TimePoint)


# repeat but with Group2 labels


# Compare Early vs Very Late
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "Very Late", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_VeryLate <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                               pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                               selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                               shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                               xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Early_VeryLate, paste('./plots/', experiment, '_explore_plot_volcano_Early_VeryLate', sep = ''), c(6, 5))



# Compare Early vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "IP", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_IP <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Early_IP, paste('./plots/', experiment, '_explore_plot_volcano_Early_IP', sep = ''), c(6, 5))


# Compare Early, Late, Very Late vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = c("Early", "Late", "Very Late"), ident.2 = "IP", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_InVivo_IP <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_InVivo_IP, paste('./plots/', experiment, '_explore_plot_volcano_InVivo_IP', sep = ''), c(6, 5))




# examine CARTEx C2 genes
cartex_C2 <- rownames(read.csv(paste(PATH_CONSTRUCTION, "./data/cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))


# Compare Early vs Very Late
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "Very Late", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, "C2")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_VeryLate_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                               pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                               selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                               shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                               xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Early_VeryLate_C2genes, paste('./plots/', experiment, '_explore_plot_volcano_Early_VeryLate_C2genes', sep = ''), c(6, 5))

# Compare Early vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "IP", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, "C2")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_IP_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Early_IP_C2genes, paste('./plots/', experiment, '_explore_plot_volcano_Early_IP_C2genes', sep = ''), c(6, 5))





####################################################################################################
##################################### Heatmap of gene expression ###################################
####################################################################################################


allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)
top100genes <- head(VariableFeatures(expt.obj.downsample), 100)


heatmap_allgenes_group <- DoHeatmap(expt.obj.downsample, group.by = "Group")
heatmap_allgenes_exCARTEx630_group <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "Group")
heatmap_top100genes_group <- DoHeatmap(expt.obj.downsample, features = top100genes, group.by = "Group")
heatmap_top100genes_seurat_clusters <- DoHeatmap(expt.obj.downsample, features = top100genes, group.by = "seurat_clusters")
heatmap_CARTEx630_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "Group")
heatmap_CARTEx630_seurat_clusters <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "seurat_clusters")
heatmap_CARTEx200_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "Group")
heatmap_CARTEx84_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "Group")
heatmap_TSG_group <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "Group")
heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")

generate_figs(heatmap_allgenes_group, paste0('./plots/', experiment, '_explore_heatmap_allgenes_group'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_group, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_group'), c(15, 12))
generate_figs(heatmap_CARTEx630_group, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_group'), c(15, 12))
generate_figs(heatmap_CARTEx200_group, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_group'), c(15, 12))
generate_figs(heatmap_CARTEx84_group, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_group'), c(15, 12))
generate_figs(heatmap_TSG_group, paste0('./plots/', experiment, '_explore_heatmap_TSG_group'), c(15, 12))
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'Group')
FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))

# https://www.mcponline.org/article/S1535-9476(20)30681-2/fulltext
SYK_substrates <- FeaturePlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'))
generate_figs(SYK_substrates, paste0('./plots/', experiment, '_explore_SYK_substrates'), c(15, 15))

vlnplot_SYK_substrates <- VlnPlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'), group.by = 'Group')
generate_figs(vlnplot_SYK_substrates, paste0('./plots/', experiment, '_explore_vlnplot_SYK_substrates'), c(15, 15))


VlnPlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2'), group.by = 'Group')
FeaturePlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2'))












