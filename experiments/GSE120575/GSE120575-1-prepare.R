#!/usr/bin/env Rscript

### Script name: GSE120575-1-prepare.R
### Description: prepare seurat object for GSE120575
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE120575'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

express <- read.csv("./data/2019-04-29-GSE120575-dat.csv", row.names=1, header=TRUE)
metadata <- read.csv("./data/GSE120575-patient-data-formatted.csv", row.names=1, header=TRUE)
texpress <- read.tcsv("./data/2019-04-29-GSE120575-dat.csv", row.names=1, header=TRUE)
modmetadata <- read.csv("./data/GSE120575-patient-data-formatted-modified.csv", row.names=1, header=TRUE)
# write.csv(texpress, "./data/2019-04-29-GSE120575-dat-transpose.csv")

# expt.obj <- CreateSeuratObject(counts = texpress, meta.data = modmetadata, min.cells = 3, min.features = 200)
# https://github.com/satijalab/seurat/issues/7545
expt.obj <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(texpress), sparse = T), meta.data = modmetadata, min.cells = 3, min.features = 200)
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

# assemble metadata
orig_id <- data.frame()
for (i in LETTERS[1:8]){for (j in 1:12){orig_id <- rbind(orig_id, paste0(i,j))}}
expt.obj@meta.data$orig.ident <- factor(expt.obj@meta.data$orig.ident, levels = unlist(unname(orig_id)))

expt.obj[["characteristics_response"]] <- modmetadata$characteristics_response
expt.obj$characteristics_response <- factor(expt.obj$characteristics_response, levels = c("Non-responder", "Responder"))
expt.obj[["characteristics_therapy"]] <- modmetadata$characteristics_therapy
expt.obj$characteristics_therapy <- factor(expt.obj$characteristics_therapy, levels = c("anti-CTLA4","anti-PD1","anti-CTLA4+PD1"))

sum(is.na(expt.obj[["characteristics_response"]]))

expt.obj@meta.data$Timepoint <- gsub("_.*$", "", modmetadata$characteristics_patinet_ID_Prebaseline_Post_ontreatment)
expt.obj@meta.data$Timepoint <- factor(expt.obj@meta.data$Timepoint, levels = c("Pre","Post"))
expt.obj@meta.data$Timepoint_Response <- paste(expt.obj@meta.data$Timepoint, expt.obj@meta.data$characteristics_response, sep = "_")
expt.obj@meta.data$Timepoint_Response <- factor(expt.obj@meta.data$Timepoint_Response, levels = c("Pre_Non-responder","Pre_Responder","Post_Non-responder","Post_Responder"))

# check levels() for metadata
# check anything NULL: unique(expt.obj@meta.data[['identity']])
check_levels(expt.obj)

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
expt.obj <- subset(expt.obj,cells=pos_ids)

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
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled', sep = ''), c(5,4))

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
expt.obj@meta.data$seurat_clusters <- expt.obj@meta.data$RNA_snn_res.0.4
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
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 8))

umap_orig_ident <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_ident, paste('./plots/', experiment, '_prepare_umap_orig_ident', sep = ''), c(10, 6))

umap_orig_ident_highlight <- DimPlotHighlightIdents(expt.obj, orig.ident, 'umap', 'blue', 0.1, 12)
generate_figs(umap_orig_ident_highlight, paste('./plots/', experiment, '_prepare_umap_orig_ident_highlight', sep = ''), c(30, 20))

umap_response <- DimPlot(expt.obj, reduction = "umap", group.by = "characteristics_response", shuffle = TRUE, seed = 123)
generate_figs(umap_response, paste('./plots/', experiment, '_prepare_umap_response', sep = ''))

umap_response_highlight <- DimPlotHighlightIdents(expt.obj, characteristics_response, 'umap', 'blue', 0.1, 2)
generate_figs(umap_response_highlight, paste('./plots/', experiment, '_prepare_umap_response_highlight', sep = ''), c(10, 6))

umap_therapy <- DimPlot(expt.obj, reduction = "umap", group.by = "characteristics_therapy", shuffle = TRUE, seed = 123)
generate_figs(umap_therapy, paste('./plots/', experiment, '_prepare_umap_therapy', sep = ''))

umap_therapy_highlight <- DimPlotHighlightIdents(expt.obj, characteristics_therapy, 'umap', 'blue', 0.1, 2)
generate_figs(umap_therapy_highlight, paste('./plots/', experiment, '_prepare_umap_therapy_highlight', sep = ''), c(8, 8))

umap_timepoint <- DimPlot(expt.obj, reduction = "umap", group.by = "Timepoint", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint, paste('./plots/', experiment, '_prepare_umap_timepoint', sep = ''))

umap_timepoint_highlight <- DimPlotHighlightIdents(expt.obj, Timepoint, 'umap', 'blue', 0.1, 2)
generate_figs(umap_timepoint_highlight, paste('./plots/', experiment, '_prepare_umap_timepoint_highlight', sep = ''), c(10, 6))

umap_timepoint_response <- DimPlot(expt.obj, reduction = "umap", group.by = "Timepoint_Response", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint_response, paste('./plots/', experiment, '_prepare_umap_timepoint_response', sep = ''))

umap_timepoint_response_highlight <- DimPlotHighlightIdents(expt.obj, Timepoint_Response, 'umap', 'blue', 0.1, 2)
generate_figs(umap_timepoint_response_highlight, paste('./plots/', experiment, '_prepare_umap_timepoint_response_highlight', sep = ''), c(8, 8))


####################################################################################################
###################################### Seurat cluster analysis #####################################
####################################################################################################

# examine metadata split by Seurat clusters

barplot_response_seurat_clusters <- BarPlotStackSplit(expt.obj, 'characteristics_response', 'seurat_clusters')
generate_figs(barplot_response_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_response_seurat_clusters', sep = ''), c(8,4))

barplot_therapy_seurat_clusters <- BarPlotStackSplit(expt.obj, 'characteristics_therapy', 'seurat_clusters')
generate_figs(barplot_therapy_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_therapy_seurat_clusters', sep = ''), c(8,4))

barplot_timepoint_seurat_clusters <- BarPlotStackSplit(expt.obj, 'Timepoint', 'seurat_clusters')
generate_figs(barplot_timepoint_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_timepoint_seurat_clusters', sep = ''), c(8,4))


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


# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()

