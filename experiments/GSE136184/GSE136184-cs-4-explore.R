#!/usr/bin/env Rscript

### Script name: GSE136184-4-explore.R
### Description: explore seurat object for GSE136184
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136184'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_cs_scored.rds', sep = ''))


umap_gene_ADGRG1 <- FeaturePlot(expt.obj, features = "ADGRG1")
generate_figs(umap_gene_ADGRG1, paste('./plots/', experiment, '_cs_prepare_umap_gene_ADGRG1', sep = ''), c(6,5) )

umap_gene_MT2A <- FeaturePlot(expt.obj, features = "MT2A")
generate_figs(umap_gene_MT2A, paste('./plots/', experiment, '_cs_prepare_umap_gene_MT2A', sep = ''), c(6,5) )


# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39

vlnplot_age_group_2_cols <- colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$AgeGroup2)))
vlnplot_age_group_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'AgeGroup2', ncol = 3, cols = vlnplot_age_group_2_cols, y.max = 4)
# generate_figs(vlnplot_age_group_exhaustion_markers, paste('./plots/', experiment, '_cs_prepare_vlnplot_age_group_exhaustion_markers', sep = ''), c(8,6))

# ENTPD1 not detected
vlnplot_age_group_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'AgeGroup2', cols = vlnplot_age_group_2_cols, y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("N", "U30", "U50", "U70", "E")),
                                                  VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'AgeGroup2', cols = vlnplot_age_group_2_cols, y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("N", "U30", "U50", "U70", "E")),
                                                  VlnPlot(expt.obj, features=c('LAG3'), group.by = 'AgeGroup2', cols = vlnplot_age_group_2_cols, y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("N", "U30", "U50", "U70", "E")),
                                                  VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'AgeGroup2', cols = vlnplot_age_group_2_cols, y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("N", "U30", "U50", "U70", "E")),
                                                  VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'AgeGroup2', cols = vlnplot_age_group_2_cols, y.max = 4)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE) + scale_x_discrete(labels=c("N", "U30", "U50", "U70", "E")))

generate_figs(vlnplot_age_group_exhaustion_markers, paste('./plots/', experiment, '_cs_prepare_vlnplot_age_group_exhaustion_markers', sep = ''), c(7,4)) 


# percentage of CARTEx detected

featplot_CARTEx_630_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'AgeGroup2', cols=vlnplot_age_group_2_cols, shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 15)) + ylim(c(-3, 5))
featplot_CARTEx_200_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'AgeGroup2', cols=vlnplot_age_group_2_cols, shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 15)) + ylim(c(-3, 5))
featplot_CARTEx_84_group <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'AgeGroup2', cols=vlnplot_age_group_2_cols, shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 15)) + ylim(c(-3, 5))

featplot_CARTEx_combined_group <- (featplot_CARTEx_630_group | featplot_CARTEx_200_group | featplot_CARTEx_84_group)
generate_figs(featplot_CARTEx_combined_group, paste('./plots/', experiment, '_explore_featplot_CARTEx_combined_group', sep = ''), c(10,4))






# examine differentiation

vlnplot_CARTEx_group_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AgeGroup2', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_group_monaco_split, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco_split', sep = ''), c(12,5)) 

vlnplt_age_group_2_cols <- colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$AgeGroup2)))
custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_group_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AgeGroup2', pt.size = 0, cols = vlnplt_age_group_2_cols) +
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4),
                                         VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                           theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                           ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_vlnplot_CARTEx_group_monaco', sep = ''), c(12,5)) 


md <- expt.obj@meta.data %>% as.data.table

swarmplot_CARTEx_group_monaco <- ggplot(md, aes(x = AgeGroup2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_group_monaco, paste('./plots/', experiment, '_prepare_swarmplot_CARTEx_group_monaco', sep = ''), c(6,5)) 


# aggregate without controls

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AgeGroup2', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$AgeGroup2 <- factor(expt.obj.agg$AgeGroup2, levels = c("Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table

aggplot_CARTEx_200_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_group_monaco_split, paste('./plots/', experiment, '_aggplot_CARTEx_200_group_monaco_split', sep = ''), c(6,5)) 



# incorporate size

md_count <- expt.obj@meta.data %>% group_by(monaco, AgeGroup2, pblabels) %>% summarize(count = n(), .groups = 'drop')
md_count$pblabels <- as.character(md_count$pblabels)
md <- md %>% left_join(md_count, by = c("monaco", "AgeGroup2", "pblabels"))



aggplot_CARTEx_200_group_monaco_split_countsized <- md %>% ggplot(aes(x = AgeGroup2, y = CARTEx_200, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c('N', 'U30', 'U50', 'U70', 'E')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_CARTEx_200_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_CARTEx_200_group_monaco_split_countsized', sep = ''), c(5,3)) 

aggplot_activation_group_monaco_split_countsized <- md %>% ggplot(aes(x = AgeGroup2, y = Activation, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c('N', 'U30', 'U50', 'U70', 'E')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_activation_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_activation_group_monaco_split_countsized', sep = ''), c(5,3)) 

aggplot_anergy_group_monaco_split_countsized <- md %>% ggplot(aes(x = AgeGroup2, y = Anergy, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c('N', 'U30', 'U50', 'U70', 'E')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_anergy_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_anergy_group_monaco_split_countsized', sep = ''), c(5,3)) 

aggplot_senescence_group_monaco_split_countsized <- md %>% ggplot(aes(x = AgeGroup2, y = Senescence, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c('N', 'U30', 'U50', 'U70', 'E')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_senescence_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_senescence_group_monaco_split_countsized', sep = ''), c(5,3)) 

aggplot_stemness_group_monaco_split_countsized <- md %>% ggplot(aes(x = AgeGroup2, y = Stemness, color = monaco, size = count)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  theme_classic() + theme(text = element_text(size = 18), axis.title.x = element_blank()) +
  scale_x_discrete(labels = c('N', 'U30', 'U50', 'U70', 'E')) + 
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(aggplot_stemness_group_monaco_split_countsized, paste('./plots/', experiment, '_aggplot_stemness_group_monaco_split_countsized', sep = ''), c(5,3)) 




####################################################################################################
#################################### UMAP score split by metadata ##################################
####################################################################################################


# expt.obj@meta.data$monaco <- plyr::mapvalues(x = expt.obj@meta.data$monaco,
#                                          from = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'),
#                                          to = c('N', 'CM', 'EM', 'TE'))
# expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('N', 'CM', 'EM', 'TE'))

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_monaco, paste('./plots/', experiment, '_cs_prepare_umap_CARTEx_200_monaco', sep = ''), c(8,2))

umap_activation_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_activation_monaco, paste('./plots/', experiment, '_cs_prepare_umap_activation_monaco', sep = ''), c(8,2))

umap_anergy_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_anergy_monaco, paste('./plots/', experiment, '_cs_prepare_umap_anergy_monaco', sep = ''), c(8,2))




umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_cs_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))

umap_CARTEx_200_group <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'AgeGroup2', split_ids = c("N", "U30", "U50", "U70", "E"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_group, paste('./plots/', experiment, '_cs_prepare_umap_CARTEx_200_group', sep = ''), c(10,2))


### 






# heatmap

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_age_group <- DoHeatmap(expt.obj.downsample, group.by = "AgeGroup2")
heatmap_allgenes_exCARTEx630_age_group <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "AgeGroup2")
heatmap_CARTEx630_age_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "AgeGroup2")
heatmap_CARTEx200_age_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "AgeGroup2")
heatmap_CARTEx84_age_group <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "AgeGroup2")
heatmap_TSG_age_group <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "AgeGroup2")

generate_figs(heatmap_allgenes_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_allgenes_age_group'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_allgenes_exCARTEx630_age_group'), c(15, 12))
generate_figs(heatmap_CARTEx630_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_CARTEx630_age_group'), c(15, 12))
generate_figs(heatmap_CARTEx200_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_CARTEx200_age_group'), c(15, 12))
generate_figs(heatmap_CARTEx84_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_CARTEx84_age_group'), c(15, 12))
generate_figs(heatmap_TSG_age_group, paste0('./plots/', experiment, '_cs_explore_heatmap_TSG_age_group'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_cs_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_cs_explore_heatmap_CARTEx200_monaco'), c(15, 12))



annotation_colors <- list("AgeGroup2"= c("Newborn" = "red",
                                         "Under 30" = "orange",
                                         "Under 50" = "gray",
                                         "Under 70" = "green",
                                         "Elderly" = "purple"),
                          "monaco" = c("Naive CD8 T cells" = "deepskyblue",
                                       "Central memory CD8 T cells" = "seagreen",
                                       "Effector memory CD8 T cells" = "darkgoldenrod",
                                       "Terminal effector CD8 T cells" = "plum3"),
                          "Phase" = c("G1" = "#089392",
                                      "S" = "#EAE29C",
                                      "G2M" = "#CF597E"))

annotation_colors <- list("AgeGroup2"= c("Newborn" = "red",
                                         "Under 30" = "orange",
                                         "Under 50" = "gray",
                                         "Under 70" = "green",
                                         "Elderly" = "purple"),
                          "monaco" = c("Naive CD8 T cells" = "deepskyblue",
                                       "Central memory CD8 T cells" = "seagreen",
                                       "Effector memory CD8 T cells" = "darkgoldenrod",
                                       "Terminal effector CD8 T cells" = "plum3"))


annotation_colors <- list("AgeGroup2"= c("Newborn" = "red",
                                         "Under 30" = "orange",
                                         "Under 50" = "gray",
                                         "Under 70" = "green",
                                         "Elderly" = "purple"),
                          "Phase" = c("G1" = "#089392",
                                      "S" = "#EAE29C",
                                      "G2M" = "#CF597E"))


test <- HeatmapWithMultiGroups(expt.obj.downsample, rownames(cartex_200_weights), annotation_colors, 'AgeGroup2')
generate_figs(as.grob(test), paste0('./plots/', experiment, '_cs_explore_heatmap_test'), c(15, 12))


