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

vlnplot_severity_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'Severity', ncol = 3, cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)

vlnplot_severity_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('LAG3'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE),
                                                 VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'Severity', cols = c('seagreen', 'steelblue', 'firebrick'), y.max = 5)+theme(axis.title.x = element_blank()) + guides(fill=FALSE))


generate_figs(vlnplot_severity_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_severity_exhaustion_markers', sep = ''), c(8,6))




VlnPlot(expt.obj, features = "CARTEx_200", group.by = 'Severity')


### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)

# Compare Severity
de_genes <- FindMarkers(expt.obj, ident.1 = "Extensive", ident.2 = "Healthy", group.by = "Severity", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 1)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_severity_extensive_healthy <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                           pCutoff = 10e-6, FCcutoff = 1, 
                                                           selectLab = rownames(signif), drawConnectors = TRUE, title = NULL, subtitle = NULL, 
                                                           xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_severity_extensive_healthy, paste('./plots/', experiment, '_plot_volcano_severity_extensive_healthy', sep = ''), c(10, 8))



# percentage of CARTEx detected

featplot_CARTEx_630_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 50)) + ylim(c(-3, 5))
featplot_CARTEx_200_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 50)) + ylim(c(-3, 5))
featplot_CARTEx_84_severity <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'Severity', cols=c('seagreen', 'steelblue', 'firebrick'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 50)) + ylim(c(-3, 5))

featplot_CARTEx_combined_severity <- (featplot_CARTEx_630_severity | featplot_CARTEx_200_severity | featplot_CARTEx_84_severity)
generate_figs(featplot_CARTEx_combined_severity, paste('./plots/', experiment, '_featplot_CARTEx_combined_severity', sep = ''), c(10,5))



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

aggplot_CARTEx_200_severity_monaco_split <- md %>% ggplot(aes(x = Severity, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_severity_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_severity_monaco_split', sep = ''), c(6,5)) 








