#!/usr/bin/env Rscript

### Script name: GSE196606-3-explore.R
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
experiment = 'GSE196606'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################


# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))


vlnplot_group_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'group', cols = c("grey", 'orange', "seagreen"), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_group_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_group_exhaustion_markers', sep = ''), c(8,6))



featplot_CARTEx_630_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'group', cols=c("grey", 'orange', "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_200_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'group', cols=c("grey", 'orange', "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_84_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'group', cols=c("grey", 'orange', "seagreen"), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 20)) + ylim(c(-3, 5))

featplot_CARTEx_combined_group <- (featplot_CARTEx_630_group | featplot_CARTEx_200_group | featplot_CARTEx_84_group)
generate_figs(featplot_CARTEx_combined_group, paste('./plots/', experiment, '_prepare_featplot_CARTEx_combined_group', sep = ''), c(10,5))

generate_figs(featplot_CARTEx_200_group, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_group', sep = ''), c(1.5,2))












# examine differentiation

vlnplot_CARTEx_group_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'group', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_group_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_group_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'group', pt.size = 0, cols = c('grey', 'orange', 'seagreen')) +
                                         theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                         ylab("CARTEx 200") + ylim(-2,4),
                                       VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                         theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                         ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_group_monaco <- ggplot(md, aes(x = group, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_group_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('group', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$group <- factor(expt.obj.agg$group, levels = c("Control", "LRBA", "NBEAL2"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

table(md$monaco)

aggplot_CARTEx_200_group_monaco_split <- md %>% ggplot(aes(x = group, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_group_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_group_monaco_split', sep = ''), c(6,5)) 



# incorporate size
md_count <- expt.obj@meta.data %>% group_by(monaco, group, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "group", "pblabels"))


aggplot_CARTEx_200_group_monaco_split_countsized <- md %>% ggplot(aes(x = group, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c('C', 'L', 'N')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_group_monaco_split_countsized', sep = ''), c(3.5,3)) 





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

umap_CARTEx_200_group <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'group', split_ids = c("CTRL", "LRBA", "NBEAL2"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_group, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_group', sep = ''), c(6,2))

# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
### 























