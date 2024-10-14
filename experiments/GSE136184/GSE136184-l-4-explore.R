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
expt.obj <- readRDS(paste('./data/', experiment, '_l_scored.rds', sep = ''))












####################################################################################################
#################################### UMAP score split by metadata ##################################
####################################################################################################


# expt.obj@meta.data$monaco <- plyr::mapvalues(x = expt.obj@meta.data$monaco,
#                                          from = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'),
#                                          to = c('N', 'CM', 'EM', 'TE'))
# expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('N', 'CM', 'EM', 'TE'))

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_200_monaco <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_monaco, paste('./plots/', experiment, '_l_prepare_umap_CARTEx_200_monaco', sep = ''), c(8,2))

umap_activation_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Activation"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_activation_monaco, paste('./plots/', experiment, '_l_prepare_umap_activation_monaco', sep = ''), c(8,2))

umap_anergy_monaco <- FeaturePlotSplitBy(expt.obj, features = c("Anergy"), split_identity = 'monaco', split_ids = c("N", "CM", "EM", "TE"), color_scale = fix.sc)
generate_figs(umap_anergy_monaco, paste('./plots/', experiment, '_l_prepare_umap_anergy_monaco', sep = ''), c(8,2))




umap_CARTEx_200_phase <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'Phase', split_ids = c("G1", "S", "G2M"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_phase, paste('./plots/', experiment, '_l_prepare_umap_CARTEx_200_phase', sep = ''), c(6,2))

umap_CARTEx_200_visit <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'visit', split_ids = c("1", "2", "3"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_visit, paste('./plots/', experiment, '_l_prepare_umap_CARTEx_200_visit', sep = ''), c(6,2))

umap_CARTEx_200_group <- FeaturePlotSplitBy(expt.obj, features = c("CARTEx_200"), split_identity = 'AgeGroup2', split_ids = c("U50", "U70", "E"), color_scale = fix.sc)
generate_figs(umap_CARTEx_200_group, paste('./plots/', experiment, '_l_prepare_umap_CARTEx_200_group', sep = ''), c(6,2))


### 









