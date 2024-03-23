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
experiment = 'GSE136874'
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

# expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_CAR <- DoHeatmap(expt.obj, group.by = "CAR")
heatmap_allgenes_exCARTEx630_CAR <- DoHeatmap(expt.obj, features = allgenes_exCARTEx, group.by = "CAR")
heatmap_CARTEx630_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_630_weights), group.by = "CAR")
heatmap_CARTEx200_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_200_weights), group.by = "CAR")
heatmap_CARTEx84_CAR <- DoHeatmap(expt.obj, features = rownames(cartex_84_weights), group.by = "CAR")
heatmap_TSG_CAR <- DoHeatmap(expt.obj, features = TSG, group.by = "CAR")

generate_figs(heatmap_allgenes_CAR, paste0('./plots/', experiment, '_explore_heatmap_allgenes_CAR'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_CAR, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx630_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx200_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_CAR'), c(15, 12))
generate_figs(heatmap_CARTEx84_CAR, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_CAR'), c(15, 12))
generate_figs(heatmap_TSG_CAR, paste0('./plots/', experiment, '_explore_heatmap_TSG_CAR'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'CAR')
# FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))
VlnPlot(expt.obj, features = c('HNF1B'), group.by = 'CAR')



