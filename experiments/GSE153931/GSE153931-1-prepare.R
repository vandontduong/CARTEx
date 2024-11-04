#!/usr/bin/env Rscript

### Script name: GSE153931-1-prepare.R
### Description: prepare seurat object for GSE153931
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE153931'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

countsData <- read.csv(file = "./data/GSE153931_cd8_t24_processed_data_umi_counts.txt", header = TRUE, row.names = 1, sep = ',')
annotatedData <- read.csv(file = "./data/GSE153931_cd8_t24_processed_data_annotations.txt", header = TRUE, row.names = 1, sep = '\t')

# assemble metadata

expt.obj <- CreateSeuratObject(counts = countsData, min.cells = 3, min.features = 200)
expt.obj@meta.data$orig.ident <- annotatedData$orig.ident
expt.obj@meta.data$origlib <- annotatedData$origlib
expt.obj@meta.data$orig.project_id <- annotatedData$orig.project_id
expt.obj@meta.data$orig.n_donors <- annotatedData$orig.n_donors
expt.obj@meta.data$percent.mt <- annotatedData$percent.mt
expt.obj@meta.data$orig.project_group <- annotatedData$orig.project_group
expt.obj@meta.data$orig.virus <- annotatedData$orig.virus
expt.obj@meta.data$orig.virus2 <- annotatedData$orig.virus2
expt.obj@meta.data$orig.asthma <- annotatedData$orig.asthma
expt.obj@meta.data$orig.peptide <- annotatedData$orig.peptide
expt.obj@meta.data$orig.donor <- annotatedData$orig.donor
expt.obj@meta.data$orig.donors <- annotatedData$orig.donors
expt.obj@meta.data$orig.severity <- annotatedData$orig.severity
expt.obj@meta.data$orig.severity_score <- annotatedData$orig.severity_score
expt.obj@meta.data$orig.severity_score <- factor(expt.obj@meta.data$orig.severity_score, levels = c('Mild', 'Moderate', 'Severe', 'Critical', 'void'))
expt.obj@meta.data$orig.severity_x <- annotatedData$orig.severity_x
expt.obj@meta.data$orig.severity_x <- factor(expt.obj@meta.data$orig.severity_x, levels = c('Mild', 'Moderate', 'Severe', 'void'))
expt.obj@meta.data$drug <- annotatedData$drug
expt.obj@meta.data$orig.diabetes <- annotatedData$orig.diabetes
expt.obj@meta.data$age <- annotatedData$age
expt.obj@meta.data$orig.sex <- annotatedData$orig.sex
expt.obj@meta.data$cardiovascular <- annotatedData$cardiovascular
expt.obj@meta.data$tumour <- annotatedData$tumour
expt.obj@meta.data$orig.hospital <- annotatedData$orig.hospital
expt.obj@meta.data$orig.unit <- annotatedData$orig.unit

expt.obj@meta.data$severity_mod <- plyr::mapvalues(x = expt.obj@meta.data$orig.severity_x,
                                              from = c('Mild', 'Moderate', 'Severe', 'void'),
                                              to = c('Mild', 'Severe', 'Severe', 'void'))
expt.obj@meta.data$severity_mod <- factor(expt.obj@meta.data$severity_mod, levels = c('Mild', 'Severe', 'void'))

# eliminate void donors
expt.obj <- SetIdent(expt.obj, value = "orig.donor")
expt.obj <- subset(expt.obj,cells=colnames(expt.obj)[Idents(expt.obj) !="void"])


# contingency tables
table(expt.obj@meta.data$orig.virus, expt.obj@meta.data$orig.virus2)



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
rownames(qc_review) <- c("All", "preQC", "postQC")
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
expt.obj@meta.data$seurat_clusters <- expt.obj@meta.data$RNA_snn_res.0.4
expt.obj@meta.data$seurat_clusters <- factor(expt.obj@meta.data$seurat_clusters, levels = SortNumStrList(unique(expt.obj@meta.data$seurat_clusters), shift = TRUE))

print("Calculating UMAPs...")
expt.obj <- RunUMAP(expt.obj, dims = 1:dim.max)

print("Calculating diffusion maps...")
diffusion_map <- RunDiffusion(expt.obj, k_int = 10)
expt.obj <- diffusion_map$atlas
saveRDS(diffusion_map$dmap, file = paste('./data/', experiment, '_dmap.rds', sep = ''))
rm(diffusion_map)

print("Saving Seurat object...")
saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

head(expt.obj)


# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''), c(6, 5))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 10))

umap_patient <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.donor", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_prepare_umap_patient', sep = ''), c(6.5, 5))

umap_orig_virus <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.virus", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c('orangered', 'seagreen', 'darkblue')) + theme(plot.title = element_blank())
generate_figs(umap_orig_virus, paste('./plots/', experiment, '_prepare_umap_orig_virus', sep = ''), c(2.8,2))

umap_orig_virus2 <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.virus2", shuffle = TRUE, seed = 123, pt.size = 0.1)
generate_figs(umap_orig_virus2, paste('./plots/', experiment, '_prepare_umap_orig_virus2', sep = ''), c(3,2))

umap_orig_severity <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.severity", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_severity, paste('./plots/', experiment, '_prepare_umap_orig_severity', sep = ''), c(6.5, 5))

umap_orig_severity_score <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.severity_score", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c('cadetblue', 'steelblue', 'indianred', 'firebrick', 'grey')) + theme(plot.title = element_blank())
generate_figs(umap_orig_severity_score, paste('./plots/', experiment, '_prepare_umap_orig_severity_score', sep = ''), c(6.5, 5))

umap_orig_severity_x <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.severity_x", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c('cadetblue', 'steelblue', 'indianred', 'grey')) + theme(plot.title = element_blank())
generate_figs(umap_orig_severity_x, paste('./plots/', experiment, '_prepare_umap_orig_severity_x', sep = ''), c(3.2,2))

umap_orig_peptide <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.peptide", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_peptide, paste('./plots/', experiment, '_prepare_umap_orig_peptide', sep = ''), c(6.5, 5))

umap_orig_hospital <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.hospital", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_hospital, paste('./plots/', experiment, '_prepare_umap_orig_hospital', sep = ''), c(6.5, 5))

umap_orig_unit <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.unit", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_unit, paste('./plots/', experiment, '_prepare_umap_orig_unit', sep = ''), c(6.5, 5))

umap_orig_diabetes <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.diabetes", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_diabetes, paste('./plots/', experiment, '_prepare_umap_orig_diabetes', sep = ''), c(6.5, 5))

umap_orig_asthma <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.asthma", shuffle = TRUE, seed = 123)
generate_figs(umap_orig_asthma, paste('./plots/', experiment, '_prepare_umap_orig_asthma', sep = ''), c(6.5, 5))

umap_cardiovascular <- DimPlot(expt.obj, reduction = "umap", group.by = "cardiovascular", shuffle = TRUE, seed = 123)
generate_figs(umap_cardiovascular, paste('./plots/', experiment, '_prepare_umap_cardiovascular', sep = ''), c(6.5, 5))

umap_tumour <- DimPlot(expt.obj, reduction = "umap", group.by = "tumour", shuffle = TRUE, seed = 123)
generate_figs(umap_tumour, paste('./plots/', experiment, '_prepare_umap_tumour', sep = ''), c(6.5, 5))


umap_severity_mod <- DimPlot(expt.obj, reduction = "umap", group.by = "severity_mod", shuffle = TRUE, seed = 123, pt.size = 0.1, cols = c('cadetblue', 'indianred', 'grey')) + theme(plot.title = element_blank())
generate_figs(umap_severity_mod, paste('./plots/', experiment, '_prepare_umap_severity_mod', sep = ''), c(2.8,2))


####################################################################################################
###################################### Seurat cluster analysis #####################################
####################################################################################################

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









