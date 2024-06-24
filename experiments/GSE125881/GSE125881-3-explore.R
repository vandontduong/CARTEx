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

vlnplot_group2_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                              VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Group2', cols = c('red', 'violetred', 'violet', 'purple'), y.max = 4)+theme(axis.title.x = element_blank()) + guides(fill=FALSE))

generate_figs(vlnplot_group2_exhaustion_markers, paste('./plots/', experiment, '_explore_vlnplot_group2_exhaustion_markers', sep = ''), c(8,6))


# percentage of CARTEx detected

featplot_CARTEx_630_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_200_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 20)) + ylim(c(-3, 5))
featplot_CARTEx_84_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'Group', cols=c('red', 'violetred', 'violet', 'purple'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 20)) + ylim(c(-3, 5))

featplot_CARTEx_combined_group <- (featplot_CARTEx_630_group | featplot_CARTEx_200_group | featplot_CARTEx_84_group)
generate_figs(featplot_CARTEx_combined_group, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_group', sep = ''), c(10,5))


# examine differentiation

vlnplot_CARTEx_group_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_group_monaco_split, paste('./plots/', experiment, '_cs_prepare_vlnplot_CARTEx_group_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_group_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Group', pt.size = 0, cols = c('red', 'violetred', 'violet', 'purple')) +
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4),
                                             VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_group_monaco, paste('./plots/', experiment, '_cs_prepare_vlnplot_CARTEx_group_monaco', sep = ''), c(12,5)) 



md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_group_monaco <- ggplot(md, aes(x = Group, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group_monaco, paste('./plots/', experiment, '_cs_prepare_swarmplot_CARTEx_group_monaco', sep = ''), c(6,5)) 

swarmplot_CARTEx_group2_monaco <- ggplot(md, aes(x = Group2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group2_monaco, paste('./plots/', experiment, '_cs_prepare_swarmplot_CARTEx_group2_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('Group2', 'monaco', 'pblabels'), return.seurat = TRUE)


# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_630 <- Z(SignatureScore(expt.obj.agg, cartex_630_weights))
expt.obj.agg@meta.data$CARTEx_630i <- integerize(expt.obj.agg@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_200 <- Z(SignatureScore(expt.obj.agg, cartex_200_weights))
expt.obj.agg@meta.data$CARTEx_200i <- integerize(expt.obj.agg@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_84 <- Z(SignatureScore(expt.obj.agg, cartex_84_weights))
expt.obj.agg@meta.data$CARTEx_84i <- integerize(expt.obj.agg@meta.data$CARTEx_84)


activation.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = TRUE, row.names = 1))
senescence.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = TRUE, row.names = 1))

expt.obj.agg <- AddModuleScore(expt.obj.agg, features = list(activation.sig, anergy.sig, stemness.sig, senescence.sig), name="State", search = TRUE)

# z score normalization
expt.obj.agg@meta.data$Activation <- scale(expt.obj.agg@meta.data$State1)
expt.obj.agg@meta.data$Anergy <- scale(expt.obj.agg@meta.data$State2)
expt.obj.agg@meta.data$Stemness <- scale(expt.obj.agg@meta.data$State3)
expt.obj.agg@meta.data$Senescence <- scale(expt.obj.agg@meta.data$State4)

expt.obj.agg@meta.data$Activationi <- integerize(expt.obj.agg@meta.data$Activation)
expt.obj.agg@meta.data$Anergyi <- integerize(expt.obj.agg@meta.data$Anergy)
expt.obj.agg@meta.data$Stemnessi <- integerize(expt.obj.agg@meta.data$Stemness)
expt.obj.agg@meta.data$Senescencei <- integerize(expt.obj.agg@meta.data$Senescence)

expt.obj.agg@meta.data$Activationi <- factor(expt.obj.agg@meta.data$Activationi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Anergyi <- factor(expt.obj.agg@meta.data$Anergyi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Stemnessi <- factor(expt.obj.agg@meta.data$Stemnessi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Senescencei <- factor(expt.obj.agg@meta.data$Senescencei, levels = c(4,3,2,1,0,-1,-2,-3,-4))

expt.obj.agg@meta.data$State1 <- NULL
expt.obj.agg@meta.data$State2 <- NULL
expt.obj.agg@meta.data$State3 <- NULL
expt.obj.agg@meta.data$State4 <- NULL

# examine other signatures
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))

expt.obj.agg <- AddModuleScore(expt.obj.agg, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

expt.obj.agg@meta.data$NKlike_Tex <- scale(expt.obj.agg@meta.data$Signature1)
expt.obj.agg@meta.data$LCMV_Tex <- scale(expt.obj.agg@meta.data$Signature2)
expt.obj.agg@meta.data$BBD_Tex <- scale(expt.obj.agg@meta.data$Signature3)
expt.obj.agg@meta.data$PD1_Tex <- scale(expt.obj.agg@meta.data$Signature4)

expt.obj.agg@meta.data$NKlike_Texi <- integerize(expt.obj.agg@meta.data$NKlike_Tex)
expt.obj.agg@meta.data$LCMV_Texi <- integerize(expt.obj.agg@meta.data$LCMV_Tex)
expt.obj.agg@meta.data$BBD_Texi <- integerize(expt.obj.agg@meta.data$BBD_Tex)
expt.obj.agg@meta.data$PD1_Texi <- integerize(expt.obj.agg@meta.data$PD1_Tex)

expt.obj.agg@meta.data$Signature1 <- NULL
expt.obj.agg@meta.data$Signature2 <- NULL
expt.obj.agg@meta.data$Signature3 <- NULL
expt.obj.agg@meta.data$Signature4 <- NULL


expt.obj.agg$Group2 <- factor(expt.obj.agg$Group2, levels = c("IP", "Early", "Late", "Very Late"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_group2_monaco_split <- md %>% ggplot(aes(x = Group2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_group2_monaco_split, paste('./plots/', experiment, '_cs_aggplot_CARTEx_200_group2_monaco_split', sep = ''), c(6,5)) 





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
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_ExpansionPeak_Late <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                   pCutoff = 10e-6, FCcutoff = 1, 
                                                   selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                   xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_ExpansionPeak_Late, paste('./plots/', experiment, '_explore_plot_volcano_ExpansionPeak_Late', sep = ''), c(10, 8))


# Compare ExpansionPeak vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = "ExpansionPeak", ident.2 = "IP", group.by = "Group", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_ExpansionPeak_IP <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 1, 
                                                 selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_ExpansionPeak_IP, paste('./plots/', experiment, '_explore_plot_volcano_ExpansionPeak_IP', sep = ''), c(10, 8))

VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'Patient', split.by = 'Group', cols = c('red', 'violetred', 'violet', 'purple')) + theme(legend.position = 'none', axis.text.x = element_text(angle = 0, hjust = 0.5)) + ylab('CARTEx 200') + xlab(NULL)



table(expt.obj@meta.data$Group, expt.obj@meta.data$TimePoint)


# repeat but with Group2 labels


# Compare Early vs Very Late
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "Very Late", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_VeryLate <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                   pCutoff = 10e-6, FCcutoff = 1, 
                                                   selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                   xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_Early_VeryLate, paste('./plots/', experiment, '_explore_plot_volcano_Early_VeryLate', sep = ''), c(10, 8))



# Compare Early vs IP
de_genes <- FindMarkers(expt.obj, ident.1 = "Early", ident.2 = "IP", group.by = "Group2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_Early_IP <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                 pCutoff = 10e-6, FCcutoff = 1, 
                                                 selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                 xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_Early_IP, paste('./plots/', experiment, '_explore_plot_volcano_Early_IP', sep = ''), c(10, 8))




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












