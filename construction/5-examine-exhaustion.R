setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")


#------------
# Data
#------------
expression_data="./data/2021-04-26-cartex-experiment-log2-TPM-data.csv"
differential_analysis_results="./data/2021-04-26-cartex-gene-signature-3317-genes.csv"
samplefile="./data/CARTEx-experiment-metadata.csv"
precomputed_cluster5_genes = "./data/2021-04-30-cartex-gene-clusters.csv"

#------------------------------------------------------------------------------------------------
# Selecting 3317 CARTEx gene-list (from supervised analysis) to perform cluster analysis
#-------------------------------------------------------------------------------------------------
data <- read.csv(expression_data, header = TRUE, row.names = 1)
cart_exhaustion_genelist <- read.csv(differential_analysis_results, header = TRUE, row.names = 1)
cart_exhaustion_genelist <- cart_exhaustion_genelist$x
cart_exhaustion_genelist <- gsub("-", ".", cart_exhaustion_genelist)
missing_changed <- c("X43162", "X43346", "X43167", "X43351", "X43168")
cart_exhaustion_genelist <- c(cart_exhaustion_genelist, missing_changed)
select = data[rownames(data) %in% cart_exhaustion_genelist, ]
sampleTable <- read.csv(samplefile,header=T, row.names=1)
sampleTable=sampleTable[order(rownames(sampleTable)),]
all(rownames(sampleTable)==colnames(select)) # check if the samplenames and order match

#------------
# Heatmap
#------------
mat <- t(scale(t(select))) # mean-center data
ann <- data.frame(Timepoint = as.factor(sampleTable$Timepoint), CAR = sampleTable$CAR)
colAnn <- HeatmapAnnotation(
  df = ann,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'), 
  annotation_legend_param = list(
    CAR = list(nrow = 2, title = 'CAR',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Timepoint = list(nrow = 3,title = 'Timepoint',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold'))
  ))
# Color heatmap
myCol <- colorRampPalette(c('blue', 'white', 'red'))(100)
myBreaks <- seq(-3, 3, length.out = 100)


#-------------------------------------
# CLUSTER ANALYSIS
#-------------------------------------
pamClusters <- cluster::pam(mat[,-c(10:12)], k = 5)
pamClusters$clustering <- paste0('C', pamClusters$clustering)


cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))

rownames(mat)
index_CARTEx_630 <- rownames(mat) %in% rownames(cartex_630_weights)
index_CARTEx_200 <- rownames(mat) %in% rownames(cartex_200_weights)
index_CARTEx_84 <- rownames(mat) %in% rownames(cartex_84_weights)
index_NKlike <- rownames(mat) %in% NK_like
index_Wherry <- rownames(mat) %in% Wherry_Tex
index_BBD <- rownames(mat) %in% BBD_Tex
index_PD1 <- rownames(mat) %in% PD1_Tex

intersect_CARTEx_630 <- rownames(mat)[index_CARTEx_630]
intersect_CARTEx_200 <- rownames(mat)[index_CARTEx_200]
intersect_CARTEx_84 <- rownames(mat)[index_CARTEx_84]
intersect_NKlike <- rownames(mat)[index_NKlike]
intersect_Wherry <- rownames(mat)[index_Wherry]
intersect_BBD <- rownames(mat)[index_BBD]
intersect_PD1 <- rownames(mat)[index_PD1]

intersect(intersect_CARTEx_630, intersect_NKlike)
intersect(intersect_Wherry, intersect_NKlike)

# as.data.frame(rownames(mat))

mat_CARTEx_630 <- mat[index_CARTEx_630,]
mat_CARTEx_200 <- mat[index_CARTEx_200,]
mat_CARTEx_84 <- mat[index_CARTEx_84,]
mat_NKlike <- mat[index_NKlike,]
mat_Wherry <- mat[intersect_Wherry,]
mat_BBD <- mat[intersect_BBD,]
mat_PD1 <- mat[intersect_PD1,]


# plot
set.seed(123)
split = c('CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','Controls','Controls','Controls','HA','HA','HA','HA','HA','HA','HA','HA','HA')
col_order = c('CD19_Donor76_Day11', 'CD19_Donor86_Day11', 'CD19_Donor90_Day11', 'CD19_Donor76_Day15', 'CD19_Donor86_Day15', 'CD19_Donor90_Day15', 'CD19_Donor76_Day21', 'CD19_Donor86_Day21', 'CD19_Donor90_Day21',
              'HA_Donor76_Day11', 'HA_Donor86_Day11', 'HA_Donor90_Day11', 'HA_Donor76_Day15', 'HA_Donor86_Day15', 'HA_Donor90_Day15', 'HA_Donor76_Day21', 'HA_Donor86_Day21', 'HA_Donor90_Day21',
              'Control_Donor76_Day0', 'Control_Donor86_Day0', 'Control_Donor90_Day0')

hmap_CARTEx_630 <- Heatmap(mat_CARTEx_630, name = 'Zscore', column_split = split, column_order = col_order,
                       row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_630), './plots/plot_heatmap_CARTEx_630', c(10,6))

hmap_CARTEx_200 <- Heatmap(mat_CARTEx_200, name = 'Zscore', column_split = split, column_order = col_order,
                           row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_200), './plots/plot_heatmap_CARTEx_200', c(10,6))

hmap_CARTEx_84 <- Heatmap(mat_CARTEx_84, name = 'Zscore', column_split = split, column_order = col_order,
                           row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_84), './plots/plot_heatmap_CARTEx_84', c(10,6))

hmap_NKlike <- Heatmap(mat_NKlike, name = 'Zscore', column_split = split, column_order = col_order,
                row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_NKlike), './plots/plot_heatmap_NKlike', c(10,6))

hmap_Wherry <- Heatmap(mat_Wherry, name = 'Zscore', column_split = split, column_order = col_order,
                       row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_Wherry), './plots/plot_heatmap_Wherry', c(10,6))

hmap_BBD <- Heatmap(mat_BBD, name = 'Zscore', column_split = split, column_order = col_order,
                       row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_BBD), './plots/plot_heatmap_BBD', c(10,6))

hmap_PD1 <- Heatmap(mat_PD1, name = 'Zscore', column_split = split, column_order = col_order,
                    row_gap = unit(1.5, "mm"), cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                    cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                    row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                    top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_PD1), './plots/plot_heatmap_PD1', c(10,6))













