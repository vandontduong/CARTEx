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
                                       "Terminal effector CD8 T cells" = "plum3"))


test <- HeatmapWithMultiGroups(expt.obj.downsample, rownames(cartex_200_weights), annotation_colors, 'AgeGroup2')
generate_figs(as.grob(test), paste0('./plots/', experiment, '_cs_explore_heatmap_test'), c(15, 12))






