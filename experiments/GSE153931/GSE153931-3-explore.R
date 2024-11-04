#!/usr/bin/env Rscript

### Script name: GSE153931-3-explore.R
### Description: integrate seurat object for GSE153931
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE153931'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

table(expt.obj@meta.data$orig.severity, expt.obj@meta.data$orig.severity_x)

table(expt.obj@meta.data$orig.severity, expt.obj@meta.data$orig.severity_score)

table(expt.obj@meta.data$orig.severity_x, expt.obj@meta.data$orig.severity_score)


# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# NT5E corresponds to CD73
# ENTPD1 corresponds to CD39

table(expt.obj@meta.data$orig.severity_x, expt.obj@meta.data$orig.virus2)

md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("orig.donor", "orig.severity_x", "orig.virus2")]





# subset the data
expt.obj.subset <- subset(expt.obj, subset = orig.severity_x != 'void')


vlnplot_severity_x_CV_exhaustion_markers <- VlnPlot(expt.obj.subset, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'orig.severity_x', ncol = 3, cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)

vlnplot_severity_x_CV_exhaustion_markers <- plot_grid(VlnPlot(expt.obj.subset, features=c('PDCD1'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                      VlnPlot(expt.obj.subset, features=c('HAVCR2'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                      VlnPlot(expt.obj.subset, features=c('LAG3'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                      VlnPlot(expt.obj.subset, features=c('CTLA4'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                      VlnPlot(expt.obj.subset, features=c('TIGIT'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                                      VlnPlot(expt.obj.subset, features=c('ENTPD1'), group.by = 'orig.severity_x', cols = c('cadetblue', 'violetred', 'indianred'), y.max = 5)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_severity_x_CV_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_severity_x_CV_exhaustion_markers', sep = ''), c(8,6))

vlnplot_severity_score_CV_exhaustion_markers <- VlnPlot(expt.obj.subset, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'orig.severity_score', ncol = 3, 
                                                        cols = c('cadetblue', 'violetred', 'indianred', 'firebrick'), y.max = 5)
generate_figs(vlnplot_severity_score_CV_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_severity_score_CV_exhaustion_markers', sep = ''), c(8,6))

vlnplot_severity_mod_CV_exhaustion_markers <- VlnPlot(expt.obj.subset, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'severity_mod', ncol = 3, 
                                                        cols = c('cadetblue', 'indianred'), y.max = 5)
generate_figs(vlnplot_severity_mod_CV_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_severity_mod_CV_exhaustion_markers', sep = ''), c(8,6))


# percentage of CARTEx detected

featplot_CARTEx_630_severity_x <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'orig.severity_x', cols=c('cadetblue', 'violetred', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))
featplot_CARTEx_200_severity_x <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'orig.severity_x', cols=c('cadetblue', 'violetred', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))
featplot_CARTEx_84_severity_x <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'orig.severity_x', cols=c('cadetblue', 'violetred', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))

featplot_CARTEx_combined_severity_x <- (featplot_CARTEx_630_severity_x | featplot_CARTEx_200_severity_x | featplot_CARTEx_84_severity_x)
generate_figs(featplot_CARTEx_combined_severity_x, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_severity_x', sep = ''), c(10,5))

featplot_CARTEx_630_severity_mod <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'severity_mod', cols=c('cadetblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 630') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))
featplot_CARTEx_200_severity_mod <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'severity_mod', cols=c('cadetblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank())+ ylab('CARTEx 200') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))
featplot_CARTEx_84_severity_mod <- FeatureScatter(expt.obj.subset, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'severity_mod', cols=c('cadetblue', 'indianred'), shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(legend.position = 'none', plot.title = element_blank()) + ylab('CARTEx 84') + xlab('% detected') + xlim(c(0, 30)) + ylim(c(-3, 5))

featplot_CARTEx_combined_severity_mod <- (featplot_CARTEx_630_severity_mod | featplot_CARTEx_200_severity_mod | featplot_CARTEx_84_severity_mod)
generate_figs(featplot_CARTEx_combined_severity_mod, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_severity_mod', sep = ''), c(10,5))

generate_figs(featplot_CARTEx_200_severity_mod, paste('./plots/', experiment, '_prepare_featplot_CARTEx_200_severity_mod', sep = ''), c(1.5,2))




# examine differentiation

vlnplot_CARTEx_severity_x_monaco_split <- VlnPlot(expt.obj.subset, features = 'CARTEx_200', group.by = 'orig.severity_x', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_severity_x_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_severity_x_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_severity_x_monaco <- plot_grid(VlnPlot(expt.obj.subset, features = 'CARTEx_200', group.by = 'orig.severity_x', pt.size = 0, cols = c('cadetblue', 'violetred', 'indianred')) +
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4),
                                         VlnPlot(expt.obj.subset, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_severity_x_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_severity_x_monaco', sep = ''), c(12,5)) 


md <- expt.obj.subset@meta.data %>% as.data.table

swarmplot_CARTEx_severity_x_monaco <- ggplot(md, aes(x = orig.severity_x, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_severity_x_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_severity_x_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj.subset@meta.data$pblabels <- PseudoBulkLabels(expt.obj.subset, 5)

expt.obj.agg <- AggregateExpression(expt.obj.subset, group.by = c('orig.severity_x', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$orig.severity_x <- factor(expt.obj.agg$orig.severity_x, levels = c("Mild", "Moderate", "Severe", "Critical"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_severity_x_monaco_split <- md %>% ggplot(aes(x = orig.severity_x, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_severity_x_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_x_monaco_split', sep = ''), c(6,5)) 


# incorporate size
md_count <- expt.obj.subset@meta.data %>% group_by(monaco, orig.severity_x, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "orig.severity_x", "pblabels"))

aggplot_CARTEx_200_severity_x_monaco_split_countsized <- md %>% ggplot(aes(x = orig.severity_x, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank()) + ylab('CARTEx')
generate_figs(aggplot_CARTEx_200_severity_x_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_x_monaco_split_countsized', sep = ''), c(6,5)) 




# REPEAT BUT WITH SEVERITY_MOD
# aggregate without controls

expt.obj.subset@meta.data$pblabels <- PseudoBulkLabels(expt.obj.subset, 5)

expt.obj.agg <- AggregateExpression(expt.obj.subset, group.by = c('severity_mod', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$severity_mod <- factor(expt.obj.agg$severity_mod, levels = c("Mild", "Severe"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_severity_mod_monaco_split <- md %>% ggplot(aes(x = severity_mod, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_severity_mod_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_mod_monaco_split', sep = ''), c(6,5)) 


# incorporate size
md_count <- expt.obj.subset@meta.data %>% group_by(monaco, severity_mod, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "severity_mod", "pblabels"))

aggplot_CARTEx_200_severity_mod_monaco_split_countsized <- md %>% ggplot(aes(x = severity_mod, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + ylab('CARTEx')+ 
  scale_x_discrete(labels = c('M', 'S')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_severity_mod_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_mod_monaco_split_countsized', sep = ''), c(3,3)) 






####################################################################################################
#################################### Differential gene expression ##################################
####################################################################################################


### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare severe vs mild
de_genes <- FindMarkers(expt.obj, ident.1 = "Severe", ident.2 = "Mild", group.by = "orig.severity_x", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_x_severe_mild <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                       pCutoff = 10e-6, FCcutoff = 0.5, 
                                                       selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                       xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_severity_x_severe_mild, paste('./plots/', experiment, '_plot_volcano_severity_x_severe_mild', sep = ''), c(10, 8))


# Compare severe/moderate vs mild
de_genes <- FindMarkers(expt.obj, ident.1 = c("Severe","Moderate"), ident.2 = "Mild", group.by = "orig.severity_x", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_x_severemoderate_mild <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                               pCutoff = 10e-6, FCcutoff = 0.5, 
                                                               selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                               xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_severity_x_severemoderate_mild, paste('./plots/', experiment, '_plot_volcano_severity_x_severemoderate_mild', sep = ''), c(10, 8))



# Compare severe vs void
de_genes <- FindMarkers(expt.obj, ident.1 = "Severe", ident.2 = "void", group.by = "orig.severity_x", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_x_severe_void <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                       pCutoff = 10e-6, FCcutoff = 0.5, 
                                                       selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                       xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_severity_x_severe_void, paste('./plots/', experiment, '_plot_volcano_severity_x_severe_void', sep = ''), c(10, 8))

# Compare mild vs void
de_genes <- FindMarkers(expt.obj, ident.1 = "Mild", ident.2 = "void", group.by = "orig.severity_x", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_x_mild_void <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                     pCutoff = 10e-6, FCcutoff = 0.5, 
                                                     selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                     xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_severity_x_mild_void, paste('./plots/', experiment, '_plot_volcano_severity_x_mild_void', sep = ''), c(10, 8))







####################################################################################################
#################################### UMAP score split by metadata ##################################
####################################################################################################

# use subset (executed earlier)
# table(expt.obj.subset$orig.severity_x)

expt.obj@meta.data$orig.severity_x <- factor(expt.obj@meta.data$orig.severity_x, levels = c('Mild', 'Moderate', 'Severe'))
expt.obj@meta.data$severity_mod <- factor(expt.obj@meta.data$severity_mod, levels = c('Mild', 'Severe'))



# expt.obj@meta.data$monaco <- plyr::mapvalues(x = expt.obj@meta.data$monaco,
#                                          from = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'),
#                                          to = c('N', 'CM', 'EM', 'TE'))
# expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('N', 'CM', 'EM', 'TE'))

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))


umap_CARTEx_200_orig_virus <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'orig.virus', split_ids = c("CV", "FLU", "RSV"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_orig_virus, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_orig_virus', sep = ''), c(6,2))



umap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_monaco, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_monaco', sep = ''), c(8,2))

umap_activation_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_activation_monaco, paste('./plots/', experiment, '_prepare_umap_activation_monaco', sep = ''), c(8,2))

umap_anergy_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_anergy_monaco, paste('./plots/', experiment, '_prepare_umap_anergy_monaco', sep = ''), c(8,2))




umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))

umap_CARTEx_200_severity_x <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'orig.severity_x', split_ids = c("Mild", "Moderate", "Severe"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_severity_x, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_severity_x', sep = ''), c(6,2))

umap_CARTEx_200_severity_mod <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'severity_mod', split_ids = c("Mild", "Severe"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_severity_mod, paste('./plots/', experiment, '_prepare_umap_CARTEx_200_severity_mod', sep = ''), c(4,2))

umap_activation_severity_mod <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'severity_mod', split_ids = c("Mild", "Severe"), color_scale = fix.sc)
generate_figs(umap_activation_severity_mod, paste('./plots/', experiment, '_prepare_umap_activation_severity_mod', sep = ''), c(4,2))


FeatureScatter(expt.obj, "CARTEx_200", "Activation", group.by = "severity_mod", split.by = "monaco", pt.size = 0.1, cols = c('cadetblue', 'indianred'))
FeatureScatter(expt.obj, "CARTEx_200", "Activation", group.by = "monaco", split.by = "severity_mod", pt.size = 0.1, cols = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))





### 





####################################################################################################
######################################### Score Correlations #######################################
####################################################################################################

expt.obj$severity_mod <- factor(expt.obj$severity_mod, levels = c("Mild", "Severe"))

table_monaco_phase <- table(expt.obj@meta.data$monaco, expt.obj@meta.data$Phase)
corr_monaco_phase <- CorrespondenceAnalysisPlot(table_monaco_phase, "Monaco", "Phase")
generate_figs(corr_monaco_phase, paste('./plots/', experiment, '_prepare_corr_monaco_phase', sep = ''), c(4,4))

table_CARTEx_200_NKlike <- table(expt.obj@meta.data$CARTEx_200i, expt.obj@meta.data$NKlike_Texi)
corr_CARTEx_200_NKlike <- CorrespondenceAnalysisPlot(table_CARTEx_200_NKlike, "CARTEx_200i", "NKlike_Texi")
generate_figs(corr_CARTEx_200_NKlike, paste('./plots/', experiment, '_prepare_corr_CARTEx_200_NKlike', sep = ''), c(4,4))



# alternative figure with marginal histograms
md <- expt.obj@meta.data %>% as.data.table

scatter_CARTEx_200_LCMV <- ggscatterhist(md, x = "CARTEx_200", y = "LCMV_Tex", color = "severity_mod", size = 0.1, alpha = 1, palette = c('cadetblue', 'indianred'), shuffle = TRUE,
                                         margin.params = list(fill = "severity_mod", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_LCMV), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_LCMV', sep = ''), c(2,2))

scatter_CARTEx_200_NKlike <- ggscatterhist(md, x = "CARTEx_200", y = "NKlike_Tex", color = "severity_mod", size = 0.1, alpha = 1, palette = c('cadetblue', 'indianred'), shuffle = TRUE,
                                           margin.params = list(fill = "severity_mod", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_NKlike), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_NKlike', sep = ''), c(2,2))

scatter_CARTEx_200_BBD <- ggscatterhist(md, x = "CARTEx_200", y = "BBD_Tex", color = "severity_mod", size = 0.1, alpha = 1, palette = c('cadetblue', 'indianred'), shuffle = TRUE,
                                        margin.params = list(fill = "severity_mod", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_BBD), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_BBD', sep = ''), c(2,2))

scatter_CARTEx_200_PD1 <- ggscatterhist(md, x = "CARTEx_200", y = "PD1_Tex", color = "severity_mod", size = 0.1, alpha = 1, palette = c('cadetblue', 'indianred'), shuffle = TRUE,
                                        margin.params = list(fill = "severity_mod", color = "black", size = 0.2), xlim = c(-3, 6), ylim = c(-3, 6), legend = "none", ggtheme = theme_classic())
generate_figs(print(scatter_CARTEx_200_PD1), paste('./plots/', experiment, '_prepare_scatter_CARTEx_200_PD1', sep = ''), c(2,2))















