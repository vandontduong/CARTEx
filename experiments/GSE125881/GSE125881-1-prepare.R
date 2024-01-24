#!/usr/bin/env Rscript

### Script name: GSE125881-1-prepare.R
### Description: prepare seurat object for GSE125881
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

counts_data <- read.csv(file = "./data/GSE125881_raw.expMatrix.csv", header = TRUE, row.names = 1)
# expt.obj <- CreateSeuratObject(counts = counts_data, project = "Kinetics", min.cells = 3, min.features = 200)
# https://github.com/satijalab/seurat/issues/7545
expt.obj <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(counts_data), sparse = T), project = "Kinetics", min.cells = 3, min.features = 200)
class(expt.obj@assays$RNA$counts)

expt.obj <- NormalizeData(expt.obj)

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

expt.obj[['identifier']] <- experiment

vlnplot_quality_control_standard_original <- ViolinPlotQC(expt.obj, c('nFeature_RNA','nCount_RNA', "percent.mt"), c(200, NA, NA), c(6000, NA, 10), 'identifier', 3)
generate_figs(vlnplot_quality_control_standard_original, paste('./plots/', experiment, '_prepare_vlnplot_quality_control_standard_original', sep = ''), c(8, 5))


# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

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

# Filter for CD8A + cells, removing CD4+
CD4_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
expt.obj <- subset(expt.obj,cells=pos_ids)

#- Add meta.data
expt.obj@meta.data$SampleID <- sapply(rownames(expt.obj@meta.data), function(x) str_split(x, '[.]')[[1]][2])
expt.obj@meta.data$SampleID <- factor(expt.obj@meta.data$SampleID, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'))

expt.obj@meta.data$Patient <- plyr::mapvalues(x = expt.obj@meta.data$SampleID,
                                                          from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                          to = c('CLL-1', 'CLL-1', 'CLL-1', 'CLL-1', 'NHL-9', 'NHL-9', 'NHL-9', 'CLL-2', 'CLL-2', 'CLL-2', 'CLL-2', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-9'))
expt.obj@meta.data$Patient <- factor(expt.obj@meta.data$Patient, levels = c('CLL-1', 'CLL-2', 'NHL-9', 'NHL-10'))

expt.obj@meta.data$TimePoint <- plyr::mapvalues(x = expt.obj@meta.data$SampleID,
                                                            from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                            to = c('IP', 'd21', 'd38', 'd112', 'd12', 'd29', 'd102', 'IP', 'd12', 'd29', 'd83', 'IP', 'd12', 'd28', 'd89', 'IP'))
expt.obj@meta.data$TimePoint <- factor(expt.obj@meta.data$TimePoint, levels = c('IP', 'd12', 'd21', 'd28', 'd29', 'd38', 'd83', 'd89', 'd102', 'd112'))

expt.obj@meta.data$Group <- plyr::mapvalues(x = expt.obj@meta.data$SampleID,
                                                        from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                        to = c('IP', 'ExpansionPeak', 'Contraction', 'Late', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP'))
expt.obj@meta.data$Group <- factor(expt.obj@meta.data$Group, levels = c('IP', 'ExpansionPeak', 'Contraction', 'Late'))

expt.obj@meta.data$Disease <- plyr::mapvalues(x = expt.obj@meta.data$SampleID,
                                                          from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                          to = c('CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'NHL', 'NHL'))
expt.obj@meta.data$Disease <- factor(expt.obj@meta.data$Disease, levels = c('CLL', 'NHL'))

# check levels() for metadata
# check anything NULL: unique(expt.obj@meta.data[['identity']])
check_levels(expt.obj)

# Examine features before quality control
vlnplot_quality_control_standard_pre <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'identifier', ncol=3)
generate_figs(vlnplot_quality_control_standard_pre, paste('./plots/', experiment, '_prepare_vlnplot_quality_control_standard_pre', sep = ''))

vlnplot_quality_control_CARTEx_pre <- VlnPlot(object = expt.obj, features = c('percent.CARTEx_630','percent.CARTEx_200', 'percent.CARTEx_84', 'PFSD.CARTEx_630', 'PFSD.CARTEx_200', 'PFSD.CARTEx_84'), group.by = 'identifier', ncol=3)
#### generate_figs(vlnplot_quality_control_CARTEx_pre, paste('./plots/', experiment, '_prepare_vlnplot_quality_control_standard_pre', sep = ''))

# capture counts before quality filter
qc_review <- rbind(qc_review, dim(expt.obj)) # [genes, cells]

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Examine features after quality control
vlnplot_quality_control_standard_post <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'identifier', ncol=3)
generate_figs(vlnplot_quality_control_standard_post, paste('./plots/', experiment, '_prepare_vlnplot_quality_control_standard_post', sep = ''))

vlnplot_quality_control_CARTEx_post <- VlnPlot(object = expt.obj, features = c('percent.CARTEx_630','percent.CARTEx_200', 'percent.CARTEx_84', 'PFSD.CARTEx_630', 'PFSD.CARTEx_200', 'PFSD.CARTEx_84'), group.by = 'identifier', ncol=3)

# capture counts after quality filter
qc_review <- rbind(qc_review, dim(expt.obj)) # [genes, cells]
rownames(qc_review) <- c("All", "preQC", "preQC")
colnames(qc_review) <- c("genes", "cells")
write.csv(qc_review, paste('./data/', experiment, '_prepare_qc_review.csv', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_varfeats <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_varfeats, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 30)

inspect_elbow <- ElbowPlot(expt.obj)
generate_figs(inspect_elbow, paste('./plots/', experiment, '_prepare_inspect_elbow', sep = ''), c(5, 4))

# select dimensionality based on elbowplot analysis
dim.max <- ChoosePC(expt.obj)
expt.obj <- FindNeighbors(expt.obj, dims = 1:dim.max)

# Examine a range of clustering resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.2)
expt.obj <- FindClusters(expt.obj, resolution = resolution.range)

inspect_clustering <- clustree(expt.obj, prefix = "RNA_snn_res.")
generate_figs(inspect_clustering, paste('./plots/', experiment, '_prepare_inspect_clustering', sep = ''), c(12, 8))

# Select clustering resolution based on clustree analysis
expt.obj@meta.data$seurat_clusters <- expt.obj@meta.data$RNA_snn_res.0.6
expt.obj@meta.data$seurat_clusters <- factor(expt.obj@meta.data$seurat_clusters, levels = SortNumStrList(unique(expt.obj@meta.data$seurat_clusters), shift = TRUE))

print("Calculating UMAPs...")
expt.obj <- RunUMAP(expt.obj, dims = 1:dim.max)

print("Calculating diffusion maps...")
expt.obj <- RunDiffusion(expt.obj, k_int = 10)

print("Saving Seurat object...")
saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

head(expt.obj)


# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''), c(6, 5))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 10))

umap_patient <- DimPlot(expt.obj, reduction = "umap", group.by = "Patient", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_prepare_umap_patient', sep = ''), c(6.5, 5))

umap_patient_highlight <- DimPlotHighlightIdents(expt.obj, Patient, 'umap', 'blue', 0.1, 2)
generate_figs(umap_patient_highlight, paste('./plots/', experiment, '_prepare_umap_patient_highlight', sep = ''), c(14, 14))

# umap_timepoint_cols <- c('IP' = 'red', 'd12' = 'violetred', 'd21' = 'violetred', 'd28' = 'violet', 'd29' = 'violet', 'd38' = 'violet', 'd83' = 'purple', 'd89' = 'purple', 'd102' = 'purple', 'd112' = 'purple')
# https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
umap_timepoint_cols <- colorRampPalette(c("red","violetred","violet","purple"))(length(unique(expt.obj@meta.data$TimePoint)))
umap_timepoint <- DimPlot(expt.obj, reduction = "umap", group.by = "TimePoint", shuffle = TRUE, seed = 123, cols = umap_timepoint_cols)
generate_figs(umap_timepoint, paste('./plots/', experiment, '_prepare_umap_timepoint', sep = ''), c(6.5, 5))

umap_timepoint_highlight <- DimPlotHighlightIdents(expt.obj, TimePoint, 'umap', 'blue', 0.1, 5)
generate_figs(umap_timepoint_highlight, paste('./plots/', experiment, '_prepare_umap_timepoint_highlight', sep = ''), c(25, 14))

umap_group_cols <- c('IP' = 'red', 'ExpansionPeak' = 'violetred', 'Contraction' = 'violet', 'Late' = 'purple')
umap_group <- DimPlot(expt.obj, reduction = "umap", group.by = "Group", shuffle = TRUE, seed = 123, cols = umap_group_cols)
generate_figs(umap_group, paste('./plots/', experiment, '_prepare_umap_group', sep = ''), c(6.5, 5))

umap_group_highlight <- DimPlotHighlightIdents(expt.obj, Group, 'umap', 'blue', 0.1, 2)
generate_figs(umap_group_highlight, paste('./plots/', experiment, '_prepare_umap_group_highlight', sep = ''), c(14, 14))

umap_disease <- DimPlot(expt.obj, reduction = "umap", group.by = "Disease", shuffle = TRUE, seed = 123)
generate_figs(umap_disease, paste('./plots/', experiment, '_prepare_umap_disease', sep = ''))

umap_disease_highlight <- DimPlotHighlightIdents(expt.obj, Disease, 'umap', 'blue', 0.1, 2)
generate_figs(umap_disease_highlight, paste('./plots/', experiment, '_prepare_umap_disease_highlight', sep = ''), c(15, 8))

# generate diffusion maps for metadata

dmap_seurat_clusters <- DimPlot(expt.obj, reduction = "dm", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(dmap_seurat_clusters, paste('./plots/', experiment, '_prepare_dmap_seurat_clusters', sep = ''), c(6, 5))

dmap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'dm', 'blue', 0.1, 4)
generate_figs(dmap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_dmap_seurat_clusters_highlight', sep = ''), c(12, 10))

dmap_patient <- DimPlot(expt.obj, reduction = "dm", group.by = "Patient", shuffle = TRUE, seed = 123)
generate_figs(dmap_patient, paste('./plots/', experiment, '_prepare_dmap_patient', sep = ''), c(6.5, 5))

dmap_patient_highlight <- DimPlotHighlightIdents(expt.obj, Patient, 'dm', 'blue', 0.1, 2)
generate_figs(dmap_patient_highlight, paste('./plots/', experiment, '_prepare_dmap_patient_highlight', sep = ''), c(14, 14))

dmap_timepoint_cols <- colorRampPalette(c("red","violetred","violet","purple"))(length(unique(expt.obj@meta.data$TimePoint)))
dmap_timepoint <- DimPlot(expt.obj, reduction = "dm", group.by = "TimePoint", shuffle = TRUE, seed = 123, cols = dmap_timepoint_cols)
generate_figs(dmap_timepoint, paste('./plots/', experiment, '_prepare_dmap_timepoint', sep = ''), c(6.5, 5))

dmap_timepoint_highlight <- DimPlotHighlightIdents(expt.obj, TimePoint, 'dm', 'blue', 0.1, 5)
generate_figs(dmap_timepoint_highlight, paste('./plots/', experiment, '_prepare_dmap_timepoint_highlight', sep = ''), c(25, 14))

dmap_group_cols <- c('IP' = 'red', 'ExpansionPeak' = 'violetred', 'Contraction' = 'violet', 'Late' = 'purple')
dmap_group <- DimPlot(expt.obj, reduction = "dmap", group.by = "Group", shuffle = TRUE, seed = 123, cols = dmap_group_cols)
generate_figs(dmap_group, paste('./plots/', experiment, '_prepare_dmap_group', sep = ''), c(6.5, 5))

dmap_group_highlight <- DimPlotHighlightIdents(expt.obj, Group, 'dm', 'blue', 0.1, 2)
generate_figs(dmap_group_highlight, paste('./plots/', experiment, '_prepare_dmap_group_highlight', sep = ''), c(14, 14))

dmap_disease <- DimPlot(expt.obj, reduction = "dm", group.by = "Disease", shuffle = TRUE, seed = 123)
generate_figs(dmap_disease, paste('./plots/', experiment, '_prepare_dmap_disease', sep = ''))

dmap_disease_highlight <- DimPlotHighlightIdents(expt.obj, Disease, 'dm', 'blue', 0.1, 2)
generate_figs(dmap_disease_highlight, paste('./plots/', experiment, '_prepare_dmap_disease_highlight', sep = ''), c(15, 8))



####################################################################################################
###################################### Seurat cluster analysis #####################################
####################################################################################################

# examine metadata split by Seurat clusters

barplot_group_seurat_clusters <- BarPlotStackSplit(expt.obj, 'Group', 'seurat_clusters')
generate_figs(barplot_group_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_group_seurat_clusters', sep = ''), c(8,4))

barplot_timepoint_seurat_clusters <- BarPlotStackSplit(expt.obj, 'TimePoint', 'seurat_clusters')
generate_figs(barplot_timepoint_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_timepoint_seurat_clusters', sep = ''), c(8,4))

barplot_patient_seurat_clusters <- BarPlotStackSplit(expt.obj, 'Patient', 'seurat_clusters')
generate_figs(barplot_patient_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_patient_seurat_clusters', sep = ''), c(8,4))

barplot_disease_seurat_clusters <- BarPlotStackSplit(expt.obj, 'Disease', 'seurat_clusters')
generate_figs(barplot_disease_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_disease_seurat_clusters', sep = ''), c(8,4))


# identify markers for each Seurat cluster
# https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers

expt.markers <- FindAllMarkers(expt.obj, only.pos = TRUE)
saveRDS(expt.markers, paste(file='./data/', experiment, '_seurat_markers.rds', sep = ''))
write.csv(expt.markers, paste(file='./data/', experiment, '_seurat_markers.csv', sep = ''))
# expt.markers <- readRDS(paste('./data/', experiment, '_seurat_markers.rds', sep = ''))

expt.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# visualize; downsampling is necessary for DoHeatmap()
# https://github.com/satijalab/seurat/issues/2724
expt.markers.top10 <- expt.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup()
cluster_markers_heatmap <- DoHeatmap(subset(expt.obj, downsample = 100), features = expt.markers.top10$gene) + NoLegend()
generate_figs(cluster_markers_heatmap, paste('./plots/', experiment, '_prepare_cluster_markers_heatmap', sep = ''), c(25, 20))

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
write.csv(cluster_markers, paste('./data/', experiment, '_prepare_cluster_markers.csv', sep = ''), row.names=FALSE)

# https://www.nature.com/articles/s12276-023-01105-x

# cluster 0 - DUSP2 is a T cell checkpoint # https://www.nature.com/articles/s41590-019-0577-9
# cluster 4 - GZMB, IL17RB indicates TC1, TC17?
# cluster 5 - NKG7 indicates NK-like

# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()
