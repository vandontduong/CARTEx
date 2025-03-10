#!/usr/bin/env Rscript

### Script name: GSE160160-3-explore.R
### Description: explore seurat object for GSE160160
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
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

vlnplot_exposure_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'exposure', ncol = 3, cols = c('skyblue', 'cadetblue'), y.max = 3)

vlnplot_exposure_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'exposure', cols = c('skyblue', 'cadetblue'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("D0", "D20")))

generate_figs(vlnplot_exposure_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_exposure_exhaustion_markers', sep = ''), c(5,3.5))

# examine BET signaling / transcription genes
VlnPlot(expt.obj, features = c('MYC', 'BRD4', 'BRD3', 'BRD2', 'BRDT', 'CDK9', 'CCNT1', 'HIF1A', 'HIF1B'), group.by = 'exposure', ncol = 3, cols = c('skyblue', 'cadetblue'), y.max = 3)

# percentage of CARTEx detected

featplot_CARTEx_630_exposure <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'exposure', cols=c('skyblue', 'cadetblue'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 50)) + ylim(c(-3, 5)) + scale_x_continuous(breaks = c(0,20,40))
featplot_CARTEx_200_exposure <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'exposure', cols=c('skyblue', 'cadetblue'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 50)) + ylim(c(-3, 5)) + scale_x_continuous(breaks = c(0,20,40))
featplot_CARTEx_84_exposure <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'exposure', cols=c('skyblue', 'cadetblue'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 50)) + ylim(c(-3, 5)) + scale_x_continuous(breaks = c(0,20,40))

featplot_CARTEx_combined_exposure <- (featplot_CARTEx_630_exposure | featplot_CARTEx_200_exposure | featplot_CARTEx_84_exposure)
generate_figs(featplot_CARTEx_combined_exposure, paste('./plots/', experiment, '_featplot_CARTEx_combined_exposure', sep = ''), c(10,4))

featplot_CARTEx_200_exposure <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'exposure', cols=c('skyblue', 'cadetblue'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected') + xlim(c(0, 50)) + ylim(c(-3, 5)) + scale_x_continuous(breaks = c(0,20,40))
generate_figs(featplot_CARTEx_200_exposure, paste('./plots/', experiment, '_featplot_CARTEx_200_exposure', sep = ''), c(1.5,2))





# aggregate for CARTEx scores (BAR CHART)


expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('exposure', 'pblabels'), return.seurat = TRUE)
expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg$exposure <- factor(expt.obj.agg$exposure, levels = c("day0", "day20"))

md <- expt.obj.agg@meta.data %>% as.data.table
glimpse(md)



agg_barplot_CARTEx_200 <- md %>% ggplot(aes(exposure, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-2,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(agg_barplot_CARTEx_200, paste('./plots/', experiment, '_agg_barplot_CARTEx_200', sep = ''), c(2.6,3))








# examine differentiation

vlnplot_CARTEx_exposure_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'exposure', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_exposure_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_exposure_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_exposure_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'exposure', pt.size = 0, cols = c('skyblue', 'cadetblue')) +
                                         theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                         ylab("CARTEx 200") + ylim(-2,4),
                                       VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                         theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                         ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_exposure_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_exposure_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_exposure_monaco <- ggplot(md, aes(x = exposure, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_exposure_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_exposure_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('exposure', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$exposure <- factor(expt.obj.agg$exposure, levels = c("day0", "day20"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)

aggplot_CARTEx_200_exposure_monaco_split <- md %>% ggplot(aes(x = exposure, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_exposure_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_exposure_monaco_split', sep = ''), c(6,5)) 



# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, exposure, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "exposure", "pblabels"))

aggplot_CARTEx_200_exposure_monaco_split_countsized <- md %>% ggplot(aes(x = exposure, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c('D0', 'D20')) + ylab('CARTEx') +
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_exposure_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_exposure_monaco_split_countsized', sep = ''), c(3,3)) 





umap_DC1rank <- FeaturePlot(expt.obj, feature = 'DC1rank', order = FALSE, pt.size = 0.1) + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c("yellow", "orange", "red", "violetred", "darkorchid"))
generate_figs(umap_DC1rank, paste('./plots/', experiment, '_prepare_umap_DC1rank', sep = ''), c(3, 2))









# differentially expressed genes

# Compare CAE to day0
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

de_genes <- FindMarkers(expt.obj, ident.1 = "day20", ident.2 = "day0", group.by = "exposure", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAEvday0 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                         selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_CAEvday0, paste('./plots/', experiment, '_explore_volcano_CAEvday0', sep = ''), c(6, 5))


FeaturePlot(expt.obj, features = c('CCL3', 'CCL4', 'PPARG', 'METRNL', 'IFNG', 'IL2RA', 'PMCH', 'EPAS1', 'RDH10'))

signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), 'CARTEx')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAEvday0_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 28, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_CAEvday0_200, paste('./plots/', experiment, '_explore_volcano_CAEvday0_200', sep = ''), c(6, 5))

### 




# examine CARTEx C2 genes
cartex_C2 <- rownames(read.csv(paste(PATH_CONSTRUCTION, "./data/cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))

de_genes <- FindMarkers(expt.obj, ident.1 = "day20", ident.2 = "day0", group.by = "exposure", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, 'C2')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAEvday0_CARTEx_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                         selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_CAEvday0_CARTEx_C2genes, paste('./plots/', experiment, '_explore_volcano_CAEvday0_CARTEx_C2genes', sep = ''), c(6, 5))












# Compare CAE exhausted to CAE not exhausted
# 3, 1, 8, 0 vs 6, 9, 7

de_genes <- FindMarkers(expt.obj, ident.1 = c(0,1,3,8), ident.2 = c(6,7,9), group.by = "seurat_clusters", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), 'C5')

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAE2groups <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                           pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL, 
                                           selectLab = rownames(signif), drawConnectors = FALSE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                           shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                           xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_CAE2groups, paste('./plots/', experiment, '_explore_volcano_CAE2groups', sep = ''), c(6, 5))

de_genes_sorted_CAE2groups <- arrange(filter(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5), avg_log2FC)
write.csv(de_genes_sorted_CAE2groups, "data/de_genes_sorted_CAE2groups.csv")


# FeaturePlot(expt.obj, features = rownames(filter(de_genes_sorted_CAE2groups, avg_log2FC < 0)))







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

umap_CARTEx_200_exposure <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'exposure', split_ids = c("D0", "D20"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_exposure, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_exposure', sep = ''), c(4,2))


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



scatter_CARTEx_200_LCMV <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'LCMV_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'exposure', cols =  c('skyblue', 'cadetblue')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_LCMV, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'NKlike_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'exposure', cols =  c('skyblue', 'cadetblue')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'BBD_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'exposure', cols =  c('skyblue', 'cadetblue')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_BBD, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'PD1_Tex', shuffle = TRUE, pt.size = 0.1, group.by = 'exposure', cols =  c('skyblue', 'cadetblue')) + 
  theme_classic() + theme(plot.title = element_blank(), legend.position="none") + xlim(c(-3, 6)) + ylim(c(-3, 6))
generate_figs(scatter_CARTEx_200_PD1, paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))




# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "exposure", size = 0.1, alpha = 1, palette = c('skyblue', 'cadetblue'), shuffle = TRUE,
                                         margin.params = list(fill = "exposure", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "exposure", size = 0.1, alpha = 1, palette = c('skyblue', 'cadetblue'), shuffle = TRUE,
                                           margin.params = list(fill = "exposure", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "exposure", size = 0.1, alpha = 1, palette = c('skyblue', 'cadetblue'), shuffle = TRUE,
                                        margin.params = list(fill = "exposure", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "exposure", size = 0.1, alpha = 1, palette = c('skyblue', 'cadetblue'), shuffle = TRUE,
                                        margin.params = list(fill = "exposure", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))










# gene expression heatmaps


cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_exposure <- DoHeatmap(expt.obj.downsample, group.by = "exposure")
heatmap_allgenes_exCARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "exposure")
heatmap_CARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "exposure")
heatmap_CARTEx200_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "exposure")
heatmap_CARTEx84_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "exposure")
heatmap_TSG_exposure <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "exposure")

generate_figs(heatmap_allgenes_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exposure'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx200_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx84_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_exposure'), c(15, 12))
generate_figs(heatmap_TSG_exposure, paste0('./plots/', experiment, '_explore_heatmap_TSG_exposure'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'exposure')
FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))

# https://www.mcponline.org/article/S1535-9476(20)30681-2/fulltext
SYK_substrates <- FeaturePlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'))
generate_figs(SYK_substrates, paste0('./plots/', experiment, '_explore_SYK_substrates'), c(15, 15))

vlnplot_SYK_substrates <- VlnPlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'), group.by = 'exposure')
generate_figs(vlnplot_SYK_substrates, paste0('./plots/', experiment, '_explore_vlnplot_SYK_substrates'), c(15, 15))





