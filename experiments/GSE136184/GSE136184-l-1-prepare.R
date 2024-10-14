#!/usr/bin/env Rscript

### Script name: GSE136184-l-1-prepare.R
### Description: prepare seurat object for GSE136184 (longitudinal study)
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
######################################## Load data and filter ######################################
####################################################################################################

expt.obj <- readRDS(paste('./data/GSE136184_seu_obj2.Rds', sep = ''))
expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

expt.obj[['identifier']] <- experiment

vlnplot_quality_control_standard_original <- ViolinPlotQC(expt.obj, c('nFeature_RNA','nCount_RNA', "percent.mt"), c(200, NA, NA), c(6000, NA, 10), 'identifier', 3)
generate_figs(vlnplot_quality_control_standard_original, paste('./plots/', experiment, '_l_prepare_vlnplot_quality_control_standard_original', sep = ''), c(8, 5))


# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_l_allgenes.csv', sep = ''))

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes
# also compute the percentage of CARTEx genes detected

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "/cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "/cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "/cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

expt.obj@meta.data$percent.CARTEx_630 <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA')[1:length(Cells(expt.obj))]
expt.obj@meta.data$PFSD.CARTEx_630 <- PercentageFeatureSetDetected(expt.obj, rownames(cartex_630_weights)) 

expt.obj@meta.data$percent.CARTEx_200 <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA')[1:length(Cells(expt.obj))]
expt.obj@meta.data$PFSD.CARTEx_200 <- PercentageFeatureSetDetected(expt.obj, rownames(cartex_200_weights)) 

expt.obj@meta.data$percent.CARTEx_84 <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA')[1:length(Cells(expt.obj))]
expt.obj@meta.data$PFSD.CARTEx_84 <- PercentageFeatureSetDetected(expt.obj, rownames(cartex_84_weights)) 

# capture counts before CD8+ T cell filter
qc_review <- dim(expt.obj) # [genes, cells]

# Filter for CD8A +  cells, removing CD4+
CD8A_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0))
expt.obj <- subset(expt.obj,cells=pos_ids)

#- Add meta.data
expt.obj@meta.data$Sex <- plyr::mapvalues(x = expt.obj@meta.data$Code,
                                                          from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                          to = c('M', 'F', 'M', 'F', 'F', 'M', 'F', 'M', 'F', 'F', 'M', 'M', 'F', 'M', 'M', 'F', 'F', 'M', 'M', 'F', 'F', 'M', 'F', 'M'))
expt.obj@meta.data$Sex <- factor(expt.obj@meta.data$Sex, levels = c('M', 'F'))

expt.obj@meta.data <- mutate(expt.obj@meta.data, Age = case_when(
  Code == 'sc_d1' ~ 0,
  Code == 'sc_d2' ~ 0,
  Code == 'sc_d3' ~ 26,
  Code == 'sc_d4' ~ 28,
  Code == 'sc_d5' ~ 37,
  Code == 'sc_d6' ~ 46,
  Code == 'sc_d7' ~ 53,
  Code == 'sc_d8' ~ 64,
  Code == 'sc_d9' ~ 68,
  Code == 'sc_d10' ~ 76,
  Code == 'sc_d11' ~ 77,
  Code == 'sc_d12' ~ 83,
  Code == 'sc_d13' ~ 83,
  Code == 'sc_d14' ~ 84,
  Code == 'sc_d15' ~ 88,
  Code == 'sc_d16' ~ 90,
  Code == 'L1' & visit == 1 ~ 30,
  Code == 'L1' & visit == 2 ~ 39,
  Code == 'L2' & visit == 1 ~ 32,
  Code == 'L2' & visit == 2 ~ 40,
  Code == 'L3' & visit == 1 ~ 48,
  Code == 'L3' & visit == 2 ~ 56,
  Code == 'L4' & visit == 1 ~ 54,
  Code == 'L4' & visit == 2 ~ 65,
  Code == 'L5' & visit == 1 ~ 60,
  Code == 'L5' & visit == 2 ~ 65,
  Code == 'L5' & visit == 3 ~ 70,
  Code == 'L6' & visit == 1 ~ 60,
  Code == 'L6' & visit == 2 ~ 71,
  Code == 'L7' & visit == 1 ~ 69,
  Code == 'L7' & visit == 2 ~ 78,
  Code == 'L8' & visit == 1 ~ 69,
  Code == 'L8' & visit == 2 ~ 77
  ))

expt.obj@meta.data <- mutate(expt.obj@meta.data, AgeGroup = case_when(
  Age < 30 ~ 'Young',
  Age < 70 ~ 'Middle',
  Age >= 69 ~ 'Old'
  ))
  
expt.obj@meta.data$AgeGroup <- factor(expt.obj@meta.data$AgeGroup, levels = c('Young', 'Middle', 'Old'))

expt.obj@meta.data <- mutate(expt.obj@meta.data, AgeGroup2 = case_when(
  Age == 0 ~ 'Newborn',
  Age < 30 ~ 'Under 30',
  Age < 50 ~ 'Under 50',
  Age < 70 ~ 'Under 70',
  Age >= 70 ~ 'Elderly'
  ))

expt.obj@meta.data$AgeGroup2 <- factor(expt.obj@meta.data$AgeGroup2, levels = c('Newborn', 'Under 30', 'Under 50', 'Under 70', 'Elderly'))

expt.obj@meta.data$Cohort <- plyr::mapvalues(x = expt.obj@meta.data$Code,
                                                      from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                      to = c('Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal'))
expt.obj@meta.data$Cohort <- factor(expt.obj@meta.data$Cohort, levels = c('Cross-sectional', 'Longitudinal'))

expt.obj@meta.data$Code <- factor(expt.obj@meta.data$Code, levels = c("L1", "L2", "L3", "L4", "L5", "L6", "L7"))

expt.obj@meta.data$visit <- factor(expt.obj@meta.data$visit, levels = c(1,2,3))


Idents(expt.obj) <- "Cohort"
expt.obj <- subset(expt.obj, idents = c("Longitudinal"))

# check levels() for metadata
# check anything NULL: unique(expt.obj@meta.data[['identity']])
check_levels(expt.obj)

# Examine features before quality control
vlnplot_quality_control_standard_pre <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'identifier', ncol=3)
generate_figs(vlnplot_quality_control_standard_pre, paste('./plots/', experiment, '_l_prepare_vlnplot_quality_control_standard_pre', sep = ''))

vlnplot_quality_control_CARTEx_pre <- VlnPlot(object = expt.obj, features = c('percent.CARTEx_630','percent.CARTEx_200', 'percent.CARTEx_84', 'PFSD.CARTEx_630', 'PFSD.CARTEx_200', 'PFSD.CARTEx_84'), group.by = 'identifier', ncol=3)
#### generate_figs(vlnplot_quality_control_CARTEx_pre, paste('./plots/', experiment, '_l_prepare_vlnplot_quality_control_standard_pre', sep = ''))

# capture counts before quality filter
qc_review <- rbind(qc_review, dim(expt.obj)) # [genes, cells]

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Examine features after quality control
vlnplot_quality_control_standard_post <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'identifier', ncol=3)
generate_figs(vlnplot_quality_control_standard_post, paste('./plots/', experiment, '_l_prepare_vlnplot_quality_control_standard_post', sep = ''))

vlnplot_quality_control_CARTEx_post <- VlnPlot(object = expt.obj, features = c('percent.CARTEx_630','percent.CARTEx_200', 'percent.CARTEx_84', 'PFSD.CARTEx_630', 'PFSD.CARTEx_200', 'PFSD.CARTEx_84'), group.by = 'identifier', ncol=3)

# capture counts after quality filter
qc_review <- rbind(qc_review, dim(expt.obj)) # [genes, cells]
rownames(qc_review) <- c("All", "preQC", "preQC")
colnames(qc_review) <- c("genes", "cells")
write.csv(qc_review, paste('./data/', experiment, '_l_prepare_qc_review.csv', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_varfeats <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_varfeats, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_l_prepare_varplt_labeled', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 30)

inspect_elbow <- ElbowPlot(expt.obj)
generate_figs(inspect_elbow, paste('./plots/', experiment, '_l_prepare_inspect_elbow', sep = ''), c(5, 4))

# select dimensionality based on elbowplot analysis
dim.max <- ChoosePC(expt.obj)
expt.obj <- FindNeighbors(expt.obj, dims = 1:dim.max)

# Examine a range of clustering resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.2)
expt.obj <- FindClusters(expt.obj, resolution = resolution.range)

inspect_clustering <- clustree(expt.obj, prefix = "RNA_snn_res.")
generate_figs(inspect_clustering, paste('./plots/', experiment, '_l_prepare_inspect_clustering', sep = ''), c(12, 8))

# Select clustering resolution based on clustree analysis
expt.obj@meta.data$seurat_clusters <- expt.obj@meta.data$RNA_snn_res.0.4
expt.obj@meta.data$seurat_clusters <- factor(expt.obj@meta.data$seurat_clusters, levels = SortNumStrList(unique(expt.obj@meta.data$seurat_clusters), shift = TRUE))

print("Calculating UMAPs...")
expt.obj <- RunUMAP(expt.obj, dims = 1:dim.max)

print("Calculating diffusion maps...")
diffusion_map <- RunDiffusion(expt.obj, k_int = 10)
expt.obj <- diffusion_map$atlas
saveRDS(diffusion_map$dmap, file = paste('./data/', experiment, '_l_dmap.rds', sep = ''))
rm(diffusion_map)

print("Saving Seurat object...")
saveRDS(expt.obj, file = paste('./data/', experiment, '_l.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_l.rds', sep = ''))

head(expt.obj)

# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_l_prepare_umap_seurat_clusters', sep = ''))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_l_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 10))

umap_age_group <- DimPlot(expt.obj, reduction = "umap", group.by = "AgeGroup", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group, paste('./plots/', experiment, '_l_prepare_umap_age_group', sep = ''))

umap_age_group_2_cols <- colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$AgeGroup2)))
umap_age_group_2 <- DimPlot(expt.obj, reduction = "umap", group.by = "AgeGroup2", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = umap_age_group_2_cols) + 
  theme(plot.title = element_blank()) + scale_color_manual(labels=c("U50", "U70", "E"), values = umap_age_group_2_cols)
generate_figs(umap_age_group_2, paste('./plots/', experiment, '_l_prepare_umap_age_group_2', sep = ''), c(2.9, 2))

umap_age_group_2_highlight <- DimPlotHighlightIdents(expt.obj, AgeGroup2, 'umap', 'blue', 0.1, 3)
generate_figs(umap_age_group_2_highlight, paste('./plots/', experiment, '_l_prepare_umap_age_group_2_highlight', sep = ''), c(12, 10))

umap_sex <- DimPlot(expt.obj, reduction = "umap", group.by = "Sex", shuffle = TRUE, seed = 123)
generate_figs(umap_sex, paste('./plots/', experiment, '_l_prepare_umap_sex', sep = ''))

umap_sex_highlight <- DimPlotHighlightIdents(expt.obj, Sex, 'umap', 'blue', 0.1, 2)
generate_figs(umap_sex_highlight, paste('./plots/', experiment, '_l_prepare_umap_sex_highlight', sep = ''), c(6, 5))

umap_visit <- DimPlot(expt.obj, reduction = "umap", group.by = "visit", shuffle = TRUE, seed = 123)
generate_figs(umap_visit, paste('./plots/', experiment, '_l_prepare_umap_visit', sep = ''))

umap_visit_highlight <- DimPlotHighlightIdents(expt.obj, visit, 'umap', 'blue', 0.1, 2)
generate_figs(umap_visit_highlight, paste('./plots/', experiment, '_l_prepare_umap_visit_highlight', sep = ''), c(6, 5))

umap_code <- DimPlot(expt.obj, reduction = "umap", group.by = "Code", shuffle = TRUE, seed = 123)
generate_figs(umap_code, paste('./plots/', experiment, '_l_prepare_umap_code', sep = ''))

umap_code_highlight <- DimPlotHighlightIdents(expt.obj, Code, 'umap', 'blue', 0.1, 2)
generate_figs(umap_code_highlight, paste('./plots/', experiment, '_l_prepare_umap_code_highlight', sep = ''), c(6, 5))

umap_cohort <- DimPlot(expt.obj, reduction = "umap", group.by = "Cohort", shuffle = TRUE, seed = 123)
generate_figs(umap_cohort, paste('./plots/', experiment, '_l_prepare_umap_cohort', sep = ''))

umap_cohort_highlight <- DimPlotHighlightIdents(expt.obj, Cohort, 'umap', 'blue', 0.1, 2)
generate_figs(umap_cohort_highlight, paste('./plots/', experiment, '_l_prepare_umap_cohort_highlight', sep = ''), c(6, 5))

# generate diffusion maps for metadata

dmap_seurat_clusters <- DimPlot(expt.obj, reduction = "dm", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(dmap_seurat_clusters, paste('./plots/', experiment, '_l_prepare_dmap_seurat_clusters', sep = ''), c(6, 5))

dmap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'dm', 'blue', 0.1, 4)
generate_figs(dmap_seurat_clusters_highlight, paste('./plots/', experiment, '_l_prepare_dmap_seurat_clusters_highlight', sep = ''), c(12, 10))

dmap_age_group <- DimPlot(expt.obj, reduction = "dm", group.by = "AgeGroup", shuffle = TRUE, seed = 123)
generate_figs(dmap_age_group, paste('./plots/', experiment, '_l_prepare_dmap_age_group', sep = ''), c(6, 5))

dmap_age_group_2 <- DimPlot(expt.obj, reduction = "dm", group.by = "AgeGroup2", shuffle = TRUE, seed = 123)
generate_figs(dmap_age_group_2, paste('./plots/', experiment, '_l_prepare_dmap_age_group_2', sep = ''), c(6, 5))

dmap_visit <- DimPlot(expt.obj, reduction = "dm", group.by = "visit", shuffle = TRUE, seed = 123)
generate_figs(dmap_visit, paste('./plots/', experiment, '_l_prepare_dmap_visit', sep = ''), c(6, 5))


####################################################################################################
###################################### Seurat cluster analysis #####################################
####################################################################################################

# examine metadata split by Seurat clusters

barplot_age_group_2_seurat_clusters <- BarPlotStackSplit(expt.obj, 'AgeGroup2', 'seurat_clusters')
generate_figs(barplot_age_group_2_seurat_clusters, paste('./plots/', experiment, '_l_prepare_barplot_age_group_2_seurat_clusters', sep = ''), c(8,4))

barplot_cohort_seurat_clusters <- BarPlotStackSplit(expt.obj, 'Cohort', 'seurat_clusters')
generate_figs(barplot_cohort_seurat_clusters, paste('./plots/', experiment, '_l_prepare_barplot_cohort_seurat_clusters', sep = ''), c(8,4))


# identify markers for each Seurat cluster
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers

expt.markers <- FindAllMarkers(expt.obj, only.pos = TRUE)
saveRDS(expt.markers, paste(file='./data/', experiment, '_l_seurat_markers.rds', sep = ''))
write.csv(expt.markers, paste(file='./data/', experiment, '_l_seurat_markers.csv', sep = ''))
# expt.markers <- readRDS(paste('./data/', experiment, '_seurat_markers.rds', sep = ''))

expt.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# visualize; downsampling is necessary for DoHeatmap()
# https://github.com/satijalab/seurat/issues/2724
expt.markers.top10 <- expt.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup()
cluster_markers_heatmap <- DoHeatmap(subset(expt.obj, downsample = 100), features = expt.markers.top10$gene) + NoLegend()
generate_figs(cluster_markers_heatmap, paste('./plots/', experiment, '_l_prepare_cluster_markers_heatmap', sep = ''), c(25, 20))

for (i in unique(expt.markers$cluster)){
  print(paste("Cluster:", i))
  print(expt.markers.top10[expt.markers.top10$cluster == i,])
}

# extract markers
cluster_markers <- data.frame()
for (i in unique(expt.markers$cluster)){
  print(paste("Cluster:", i))
  putative_markers = expt.markers.top10[expt.markers.top10$cluster == i,]$gene
  print(paste(unlist(putative_markers), collapse=', '))
  cluster_markers <- rbind(cluster_markers, cbind(i, paste(unlist(putative_markers), collapse=', ')))
  cat("\n")
}
colnames(cluster_markers) <- c('Cluster', 'Markers')
write.csv(cluster_markers, paste('./data/', experiment, '_l_prepare_cluster_markers.csv', sep = ''), row.names=FALSE)

# https://www.nature.com/articles/s12276-023-01105-x

# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()



