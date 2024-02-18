#!/usr/bin/env Rscript

### Script name: mdandersonTCM-3-explore.R
### Description: explore seurat object for mdandersonTCM
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'mdandersonTCM'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_CARTEx630_TissueType <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "TissueType")
heatmap_CARTEx84_TissueType <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "TissueType")
heatmap_TSG_TissueType <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "TissueType")

generate_figs(heatmap_CARTEx630_TissueType, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_TissueType'), c(15, 12))
generate_figs(heatmap_CARTEx84_TissueType, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_TissueType'), c(15, 12))
generate_figs(heatmap_TSG_TissueType, paste0('./plots/', experiment, '_explore_heatmap_TSG_TissueType'), c(15, 12))




