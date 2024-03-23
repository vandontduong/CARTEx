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
experiment = 'GSE151511'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_responder <- DoHeatmap(expt.obj.downsample, group.by = "Responder")
heatmap_allgenes_exCARTEx630_responder <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "Responder")
heatmap_CARTEx630_responder <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "Responder")
heatmap_CARTEx200_responder <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "Responder")
heatmap_CARTEx84_responder <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "Responder")
heatmap_TSG_responder <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "Responder")

generate_figs(heatmap_allgenes_responder, paste0('./plots/', experiment, '_explore_heatmap_allgenes_Responder'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_responder, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_Responder'), c(15, 12))
generate_figs(heatmap_CARTEx630_responder, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_Responder'), c(15, 12))
generate_figs(heatmap_CARTEx200_responder, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_Responder'), c(15, 12))
generate_figs(heatmap_CARTEx84_responder, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_Responder'), c(15, 12))
generate_figs(heatmap_TSG_responder, paste0('./plots/', experiment, '_explore_heatmap_TSG_Responder'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'Responder')
FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))
VlnPlot(expt.obj, features = c('HNF1A'), group.by = 'Responder')
VlnPlot(expt.obj, features = c('HNF1B'), group.by = 'Responder')



