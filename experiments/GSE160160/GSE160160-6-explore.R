#!/usr/bin/env Rscript

### Script name: GSE160160-4-explore.R
### Description: explore seurat object for GSE160160
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
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

heatmap_allgenes_exposure <- DoHeatmap(expt.obj.downsample, group.by = "exposure")
heatmap_allgenes_exCARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "exposure")
heatmap_CARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "exposure")
heatmap_CARTEx200_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "exposure")
heatmap_CARTEx84_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "exposure")
heatmap_TSG_exposure <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "exposure")

generate_figs(heatmap_allgenes_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exposure'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx200_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx84_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_exposure'), c(15, 12))
generate_figs(heatmap_TSG_exposure, paste0('./plots/', experiment, '_explore_heatmap_TSG_exposure'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'exposure')
FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))

# https://www.mcponline.org/article/S1535-9476(20)30681-2/fulltext
SYK_substrates <- FeaturePlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'))
generate_figs(SYK_substrates, paste0('./plots/', experiment, '_explore_SYK_substrates'), c(15, 15))

vlnplot_SYK_substrates <- VlnPlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'), group.by = 'exposure')
generate_figs(vlnplot_SYK_substrates, paste0('./plots/', experiment, '_explore_vlnplot_SYK_substrates'), c(15, 15))




