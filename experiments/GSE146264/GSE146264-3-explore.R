#!/usr/bin/env Rscript

### Script name: GSE146264-3-explore.R
### Description: explore seurat object for GSE146264
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE146264'
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


vlnplot_severity_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + guides(fill=FALSE) + scale_x_discrete(labels = c('H', 'M', 'E')))

generate_figs(vlnplot_severity_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_severity_exhaustion_markers', sep = ''), c(5,4))



VlnPlot(expt.obj, features = "CARTEx_200", group.by = 'Severity')


featplot_CARTEx_630_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))
featplot_CARTEx_200_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))
featplot_CARTEx_84_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))

featplot_CARTEx_combined_severity <- (featplot_CARTEx_630_severity | featplot_CARTEx_200_severity | featplot_CARTEx_84_severity)
generate_figs(featplot_CARTEx_combined_severity, paste('./plots/', experiment, '_prepare_featplot_CARTEx_combined_severity', sep = ''), c(10,5))

featplot_CARTEx_200_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected') + xlim(c(0, 40)) + ylim(c(-3, 6))
generate_figs(featplot_CARTEx_200_severity, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_severity', sep = ''), c(1.5,2))






### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare Severe vs healthy
de_genes <- FindMarkers(expt.obj, ident.1 = "Extensive", ident.2 = "Healthy", group.by = "Severity", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_extensive_healthy <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                           pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                           selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 20, 
                                                           shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                           xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_severity_extensive_healthy, paste('./plots/', experiment, '_plot_volcano_severity_extensive_healthy', sep = ''), c(10, 8))




# Compare Severe vs Mild
de_genes <- FindMarkers(expt.obj, ident.1 = "Extensive", ident.2 = "Moderate", group.by = "Severity", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_extensive_moderate <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                            pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                            selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 20, 
                                                            shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                            xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_severity_extensive_moderate, paste('./plots/', experiment, '_plot_volcano_severity_extensive_moderate', sep = ''), c(10, 8))



# Compare Severe vs healthy
de_genes <- FindMarkers(expt.obj, ident.1 = c("Extensive", "Moderate"), ident.2 = "Healthy", group.by = "Severity", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_extensive_moderate_healthy <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                           pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                           selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 20, 
                                                           shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                           xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_severity_extensive_moderate_healthy, paste('./plots/', experiment, '_plot_volcano_severity_extensive_moderate_healthy', sep = ''), c(10, 8))



###
###



# examine differentiation

vlnplot_CARTEx_severity_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Severity', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_severity_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_severity_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_severity_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Severity', pt.size = 0, cols = c('seagreen', 'steelblue', 'firebrick')) +
                                              theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                              ylab("CARTEx 200") + ylim(-2,4),
                                            VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                              theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                              ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_severity_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_severity_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_severity_monaco <- ggplot(md, aes(x = Severity, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_severity_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_severity_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('Severity', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$Severity <- factor(expt.obj.agg$Severity, levels = c("Healthy", "Moderate", "Extensive"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)


# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, Severity, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "Severity", "pblabels"))

# not enough naive cells in pseudo-bulk
aggplot_CARTEx_200_severity_monaco_split_countsized <- md %>% ggplot(aes(x = Severity, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab("CARTEx") +
  scale_x_discrete(labels = c('C', 'M', 'S')) + 
  scale_color_manual(labels=c("CM", "EM", "TE"), values = c('seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_severity_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_monaco_split_countsized', sep = ''), c(3.5,3)) 


# not enough naive cells in pseudo-bulk
aggplot_TSR_severity_monaco_split_countsized <- md %>% ggplot(aes(x = Severity, y = TSR, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab("Stress Response") +
  scale_x_discrete(labels = c('C', 'M', 'S')) + 
  scale_color_manual(labels=c("CM", "EM", "TE"), values = c('seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_TSR_severity_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_TSR_severity_monaco_split_countsized', sep = ''), c(3.5,3)) 








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

umap_CARTEx_200_Severity <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Severity', split_ids = c('C', 'M', 'S'), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_Severity, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_Severity', sep = ''), c(6,2))

# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 




umap_TSR_Severity <- FeaturePlotSplitBy(expt.obj, features = c("TSR"), split_identity = 'Severity', split_ids = c('C', 'M', 'S'), color_scale = fix.sc)







