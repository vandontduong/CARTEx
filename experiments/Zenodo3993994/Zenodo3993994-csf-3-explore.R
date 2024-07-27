#!/usr/bin/env Rscript

### Script name: Zenodo3993994-csf-3-explore.R
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

expt.obj <- readRDS(paste('./data/', experiment, '_csf_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39

vlnplot_sample_type_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'sampleType', ncol = 3, cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)

# custom_labels <- c('Parkinson (B)', 'Control (B)', 'Alzheimer (C)', 'MCI (C)', 'Control (C)')
custom_labels <- c('Alzheimer (C)', 'MCI (C)', 'Control (C)')

vlnplot_sample_type_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('LAG3'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels),
                                                    VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'sampleType', cols = c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_sample_type_exhaustion_markers, paste('./plots/', experiment, '_csf_explore_vlnplot_sample_type_exhaustion_markers', sep = ''), c(8,6))



# percentage of CARTEx detected

featplot_CARTEx_630_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 20)) + ylim(c(-1, 6))
featplot_CARTEx_200_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 20)) + ylim(c(-1, 6))
featplot_CARTEx_84_sample_type <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'sampleType', cols=c('Parkinson(Blood)' = 'violet', 'Control(Blood)' = 'cadetblue', 'Alzheimer(CSF)' = 'firebrick', 'MCI(CSF)' = 'indianred', 'Control(CSF)' = 'lightsalmon'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 20)) + ylim(c(-1, 6))

featplot_CARTEx_combined_sample_type <- (featplot_CARTEx_630_sample_type | featplot_CARTEx_200_sample_type | featplot_CARTEx_84_sample_type)
generate_figs(featplot_CARTEx_combined_sample_type, paste('./plots/', experiment, '_csf_featplot_CARTEx_combined_sample_type', sep = ''), c(10,5))



# examine differentiation

vlnplot_CARTEx_disease_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Disease', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_disease_monaco_split, paste('./plots/', experiment, '_csf_prepare_vlnplot_CARTEx_disease_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_disease_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Disease', pt.size = 0, cols = c("lightsteelblue", "steelblue", "palevioletred", "violetred")) +
                                             theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                             ylab("CARTEx 200") + ylim(-2,4),
                                           VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                             theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                             ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_disease_monaco, paste('./plots/', experiment, '_csf_prepare_vlnplot_CARTEx_disease_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_disease_monaco <- ggplot(md, aes(x = Disease, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_disease_monaco, paste('./plots/', experiment, '_csf_prepare_swarmplot_CARTEx_disease_monaco', sep = ''), c(6,5)) 


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
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_disease_monaco_split, paste('./plots/', experiment, '_csf_aggplot_CARTEx_200_disease_monaco_split', sep = ''), c(6,5)) 




####################################################################################################
#################################### Differential gene expression ##################################
####################################################################################################

### label CARTEx in exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

# Compare AD/MCI vs Control in blood
de_genes <- FindMarkers(expt.obj, ident.1 = c("Alzheimer(CSF)", "MCI(CSF)"), ident.2 = "Control(CSF)", group.by = "sampleType", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_AD_MCI_CSF <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                           pCutoff = 10e-6, FCcutoff = 0.5, 
                                           selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                           xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_AD_MCI_CSF, paste('./plots/', experiment, '_csf_explore_plot_volcano_AD_MCI_CSF', sep = ''), c(10, 8))





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

generate_figs(umap_clonotype_count, paste('./plots/', experiment, '_csf_prepare_umap_clonotype_count', sep = ''), c(5.5,5))
generate_figs(umap_clonotype_count_log2, paste('./plots/', experiment, '_csf_prepare_umap_clonotype_count_log2', sep = ''), c(5.5,5))


VlnPlot(expt.obj, features = "clonotype_count_log2", group.by = 'Disease')
VlnPlot(expt.obj, features = "CARTEx_200", group.by = 'Disease')

