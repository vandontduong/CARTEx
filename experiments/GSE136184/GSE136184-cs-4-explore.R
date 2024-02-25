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


# https://divingintogeneticsandgenomics.com/post/enhancement-of-scrnaseq-heatmap-using-complexheatmap/
# https://bioinformatics.stackexchange.com/questions/4851/how-i-can-reproduce-this-heat-map

meta_data <- expt.obj.downsample@meta.data
ordered_meta_data <- meta_data[order(meta_data$monaco), ]

annotation_colors <- list("AgeGroup2"= c("Newborn" = "red",
                                         "Under 30" = "orange",
                                         "Under 50" = "gray",
                                         "Under 70" = "green",
                                         "Elderly" = "purple"),
                          "monaco" = c("Naive CD8 T cells" = "deepskyblue",
                                       "Central memory CD8 T cells" = "seagreen",
                                       "Effector memory CD8 T cells" = "darkgoldenrod",
                                       "Terminal effector CD8 T cells" = "plum3"))

temp_df <- select(ordered_meta_data, c('AgeGroup2','monaco'))
rownames(temp_df) <- NULL

ha = HeatmapAnnotation(df = temp_df,
                       show_annotation_name = TRUE,
                       col = annotation_colors)



my_data <- expt.obj.downsample@assays$RNA@scale.data
my_data <- my_data[, rownames(ordered_meta_data)]

col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))

test <- Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  use_raster = TRUE,
  # raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha
)

generate_figs(as.grob(test), paste0('./plots/', experiment, '_cs_explore_heatmap_test'), c(15, 12))


match(common, rownames(as.matrix(GetAssayData(atlas))))
common <- intersect(rownames(cartex_200_weights), rownames(expt.obj))
common <- intersect(rownames(sig_weights), rownames(atlas))
match(common, rownames(expt.obj))

my_data[match(common, rownames(expt.obj)),]

HeatmapWithMultiGroups <- function(atlas, gene_list, anno_colors, sort_by_group){
  meta_data <- atlas@meta.data
  # sort_by_group <- names(anno_colors)[anno_i]
  ordered_meta_data <- meta_data[order(meta_data[[sort_by_group]]), ]
  
  common <- intersect(gene_list, rownames(atlas))
  
  scale_data <- atlas@assays$RNA@scale.data
  scale_data <- scale_data[match(common, rownames(atlas)), rownames(ordered_meta_data)]
  
  ordered_meta_data <- select(ordered_meta_data, names(anno_colors))
  rownames(ordered_meta_data) <- NULL
  
  heatmap_anno = HeatmapAnnotation(df = ordered_meta_data,
                         show_annotation_name = TRUE,
                         col = annotation_colors)
  
  col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  
  heatmap_plt <- Heatmap(scale_data, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                         column_order = NULL, show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE,
                         use_raster = TRUE, bottom_annotation = NULL, top_annotation = heatmap_anno)
  
  
}

test <- HeatmapWithMultiGroups(expt.obj.downsample, rownames(cartex_200_weights), annotation_colors, 'AgeGroup2')
generate_figs(as.grob(test), paste0('./plots/', experiment, '_cs_explore_heatmap_test'), c(15, 12))



HeatmapSubset <- function(atlas, gene_list, meta_group){
  common <- intersect(gene_list, rownames(atlas))
  expr <- t(as.matrix(GetAssayData(atlas))[match(common, rownames(as.matrix(GetAssayData(atlas)))),])
  expr <- scale(expr)
  # cluster_anno <- atlas@meta.data[[meta_group]]
  
  
  
  ## make the black color map to 0. the yellow map to highest and the purle map to the lowest
  col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  
  Heatmap(expr, name = "Expression",  
          # column_split = factor(cluster_anno),
          cluster_columns = TRUE,
          show_column_dend = FALSE,
          cluster_column_slices = TRUE,
          column_title_gp = gpar(fontsize = 8),
          column_gap = unit(0.5, "mm"),
          cluster_rows = TRUE,
          show_row_dend = FALSE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 4),
          column_title_rot = 90,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
          show_column_names = FALSE,
          use_raster = TRUE,
          raster_quality = 4)
}

test <- HeatmapSubset(expt.obj.downsample, rownames(cartex_630_weights), 'monaco')

factor(expt.obj@meta.data[["AgeGroup2"]])
factor(expt.obj@meta.data$AgeGroup2)

factor(expt.obj.downsample@meta.data$AgeGroup2)

common <- intersect(rownames(cartex_630_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
expr <- scale(expr)

length(factor(expt.obj@meta.data[["AgeGroup2"]]))
dim(expr)




