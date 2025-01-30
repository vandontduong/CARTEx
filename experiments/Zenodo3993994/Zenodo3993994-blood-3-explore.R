#!/usr/bin/env Rscript

### Script name: Zenodo3993994-blood-3-explore.R
### Description: annotate seurat object for Zenodo3993994
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'Zenodo3993994'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_blood_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39

vlnplot_sample_type_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'sampleType', ncol = 3, cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)

# custom_labels <- c('Parkinson (B)', 'Control (B)', 'Alzheimer (C)', 'MCI (C)', 'Control (C)')
custom_labels <- c('Parkinson', 'Control')

vlnplot_sample_type_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('LAG3'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_sample_type_exhaustion_markers, paste('./plots/', experiment, '_blood_explore_vlnplot_sample_type_exhaustion_markers', sep = ''), c(8,6))



# percentage of CARTEx detected

featplot_CARTEx_630_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_200_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_84_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 20)) + ylim(c(-3, 5))

featplot_CARTEx_combined_sample_type <- (featplot_CARTEx_630_sample_type | featplot_CARTEx_200_sample_type | featplot_CARTEx_84_sample_type)
generate_figs(featplot_CARTEx_combined_sample_type, paste('./plots/', experiment, '_blood_featplot_CARTEx_combined_sample_type', sep = ''), c(10,4))

featplot_CARTEx_200_sample_type <- featplot_CARTEx_200_sample_type + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx') + xlab('% detected')
generate_figs(featplot_CARTEx_200_sample_type, paste('./plots/', experiment, '_blood_prepare_featplot_CARTEx_200_sample_type', sep = ''), c(1.5,2))




# examine differentiation

vlnplot_CARTEx_disease_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Disease', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_disease_monaco_split, paste('./plots/', experiment, '_blood_prepare_vlnplot_CARTEx_disease_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_disease_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Disease', pt.size = 0, cols = c('Parkinson (B)' = 'violet', 'Control (B)' = 'cadetblue')) +
                                             theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                             ylab("CARTEx 200") + ylim(-2,4),
                                           VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                             theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                             ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_disease_monaco, paste('./plots/', experiment, '_blood_prepare_vlnplot_CARTEx_disease_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_disease_monaco <- ggplot(md, aes(x = Disease, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_disease_monaco, paste('./plots/', experiment, '_blood_prepare_swarmplot_CARTEx_disease_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('Disease', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$Disease <- factor(expt.obj.agg$Disease, levels = c('Parkinson (B)', 'Control (B)', 'Alzheimer (C)', 'MCI (C)', 'Control (C)'))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_disease_monaco_split <- md %>% ggplot(aes(x = Disease, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_disease_monaco_split, paste('./plots/', experiment, '_blood_aggplot_CARTEx_200_disease_monaco_split', sep = ''), c(6,5)) 



# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, Disease, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "Disease", "pblabels"))

aggplot_CARTEx_200_disease_monaco_split_countsized <- md %>% ggplot(aes(x = Disease, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab("CARTEx") +
  scale_x_discrete(labels = c('P', 'C')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_disease_monaco_split_countsized, paste('./plots/', experiment, '_blood_aggplot_CARTEx_200_disease_monaco_split_countsized', sep = ''), c(4,3.5)) 






####################################################################################################
#################################### Differential gene expression ##################################
####################################################################################################

### label CARTEx in exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))


# Compare Parkinson vs Control in blood
de_genes <- FindMarkers(expt.obj, ident.1 = "Parkinson(Blood)", ident.2 = "Control(Blood)", group.by = "sampleType", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Parkinson_blood <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Parkinson_blood, paste('./plots/', experiment, '_blood_explore_plot_volcano_Parkinson_blood', sep = ''), c(6, 5))



signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), "CARTEx")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Parkinson_blood_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_Parkinson_blood_200, paste('./plots/', experiment, '_blood_explore_plot_volcano_Parkinson_blood_200', sep = ''), c(6, 5))












md_count <- expt.obj@meta.data %>% group_by(clonotype) %>% summarize(count = n(), .groups = 'drop')
rnames <- colnames(expt.obj)
md <- expt.obj@meta.data %>% as.data.table
md <- md %>% left_join(md_count, by = c("clonotype"))
colnames(expt.obj) <- rnames

expt.obj@meta.data$clonotype_count <- md$count
expt.obj@meta.data$clonotype_count_log2 <- log2(md$count)

# expt.obj@meta.data$count
# expt.obj@meta.data$clonotype
# sum(is.na(expt.obj@meta.data$clonotype))
# expt.obj@meta.data$count[is.na(expt.obj@meta.data$count)] <- 0

umap_clonotype_count <- FeaturePlot(expt.obj, features = c("clonotype_count"), order = TRUE)
umap_clonotype_count_log2 <- FeaturePlot(expt.obj, features = c("clonotype_count_log2"), order = TRUE)
# scale_fill_continuous(low="thistle2", high="darkred", guide="colorbar",na.value="white")

generate_figs(umap_clonotype_count, paste('./plots/', experiment, '_blood_prepare_umap_clonotype_count', sep = ''), c(5.5,5))
generate_figs(umap_clonotype_count_log2, paste('./plots/', experiment, '_blood_prepare_umap_clonotype_count_log2', sep = ''), c(5.5,5))


VlnPlot(expt.obj, features = "clonotype_count_log2", group.by = 'Disease')
VlnPlot(expt.obj, features = "CARTEx_200", group.by = 'Disease')








####################################################################################################
#################################### UMAP score split by metadata ##################################
####################################################################################################


table(expt.obj@meta.data$sampleType)
table(expt.obj@meta.data$Disease)
expt.obj@meta.data$Disease <- plyr::mapvalues(x = expt.obj@meta.data$Disease,
                                             from = c('Parkinson (B)', 'Control (B)', 'Alzheimer (C)', 'MCI (C)', 'Control (C)'),
                                             to = c('P', 'C', 'A', 'M', 'C_C'))
expt.obj@meta.data$Disease <- factor(expt.obj@meta.data$Disease, levels = c('P', 'C'))

# to = c('Parkinson', 'Control', 'Alzheimer', 'MCI', 'Control'))

# expt.obj@meta.data$monaco <- plyr::mapvalues(x = expt.obj@meta.data$monaco,
#                                          from = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'),
#                                          to = c('N', 'CM', 'EM', 'TE'))
# expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('N', 'CM', 'EM', 'TE'))


fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_monaco, paste('./plots/', experiment, '_blood_prepare_umap_CARTEx_200_monaco', sep = ''), c(8,2))

umap_activation_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_activation_monaco, paste('./plots/', experiment, '_blood_prepare_umap_activation_monaco', sep = ''), c(8,2))

umap_anergy_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_anergy_monaco, paste('./plots/', experiment, '_blood_prepare_umap_anergy_monaco', sep = ''), c(8,2))




umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_blood_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))

umap_CARTEx_200_disease <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Disease', split_ids = c('P', 'C'), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_disease, paste('./plots/', experiment, '_blood_prepare_umap_CARTEx_200_disease', sep = ''), c(4,2))

# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 











