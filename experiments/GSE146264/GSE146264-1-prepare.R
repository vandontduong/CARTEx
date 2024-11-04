#!/usr/bin/env Rscript

### Script name: GSE146264-1-prepare.R
### Description: prepare seurat object for GSE146264
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE146264'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

library(org.Hs.eg.db) ## Annotation DB

counts_data <- data.frame(fread("./data/GSE146264_CountTable_SingleCellPsoCd8s_20200220.tsv", sep = "\t"), row.names = 1)
counts_data$ENSEMBL <- gsub('\\..+$', '', rownames(counts_data))

ensemble <- rownames(counts_data)
ensemble <- gsub('\\..+$', '', ensemble)
gene_symbol <- select(org.Hs.eg.db,keys = ensemble,columns = "SYMBOL",keytype = "ENSEMBL")

mapping <- match(counts_data$ENSEMBL,gene_symbol$ENSEMBL)
counts_data$Gene <- gene_symbol$SYMBOL[mapping]
counts_data <- counts_data[!duplicated(counts_data$Gene),]
counts_data <- counts_data[!is.na(counts_data$Gene),]
rownames(counts_data) <- counts_data$Gene
counts_data <- counts_data[1:(ncol(counts_data)-2)]

expt.obj <- CreateSeuratObject(counts = counts_data, project = "Psoraisis", min.cells = 3, min.features = 200)

expt.obj <- NormalizeData(expt.obj)

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT.")
grep("^MT.", rownames(expt.obj[["RNA"]]),value = T) # normally MT- or mt-
# https://github.com/satijalab/seurat/issues/1582

expt.obj[['identifier']] <- experiment

vlnplot_quality_control_standard_original <- ViolinPlotQC(expt.obj, c('nFeature_RNA','nCount_RNA', "percent.mt"), c(200, NA, NA), c(6000, NA, 10), 'identifier', 3)
generate_figs(vlnplot_quality_control_standard_original, paste('./plots/', experiment, '_prepare_vlnplot_quality_control_standard_original', sep = ''), c(8, 5))


#- Add meta.data
expt.obj@meta.data$orig.ident <- gsub('(.*?)\\.', '', colnames(expt.obj))
expt.obj@meta.data$subject <- gsub('(?:[^.]+\\.){2}([^.]+).*', '\\1', colnames(expt.obj))


expt.obj@meta.data$PASI <- as.numeric(plyr::mapvalues(x = expt.obj@meta.data$subject,
                                                                  from = c('P5', 'P9', 'P10', 'P11', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'C4', 'C5', 'C6', 'C7', 'C13', 'C14'),
                                                                  to = c(4.8, 14.9, 16.1, 8.2, 4.8, 15.2, 17.5, 16.2, 16, 14.6, 9.1, 0, 0, 0, 0, 0, 0)))
expt.obj@meta.data$PASIi <- integerize(expt.obj@meta.data$PASI)

# PASI score > 10 is considered extensive psoriasis
expt.obj@meta.data$Severity <- ifelse(expt.obj@meta.data$PASI > 10, 'Extensive', ifelse(expt.obj@meta.data$PASI != 0, 'Moderate', 'Healthy'))
expt.obj@meta.data$Severity <- factor(expt.obj@meta.data$Severity, levels = c('Healthy', 'Moderate', 'Extensive'))

expt.obj@meta.data$PGA <- plyr::mapvalues(x = expt.obj@meta.data$subject,
                                                      from = c('P5', 'P9', 'P10', 'P11', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'C4', 'C5', 'C6', 'C7', 'C13', 'C14'),
                                                      to = c(2, 3, 4, 3, 2, 4, 4, 3, 4, 4, 4, 0, 0, 0, 0, 0, 0))
# expt.obj@meta.data$PGA <- factor(expt.obj@meta.data$PGA, levels = c(0, 2, 3, 4))

expt.obj@meta.data$Onset <- as.numeric(plyr::mapvalues(x = expt.obj@meta.data$subject,
                                                                   from = c('P5', 'P9', 'P10', 'P11', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'C4', 'C5', 'C6', 'C7', 'C13', 'C14'),
                                                                   to = c(32, 5, 16, 30, 14, 36, 26, 14, 56, 25, 13, 0, 0, 0, 0, 0, 0)))

expt.obj@meta.data$Duration <- as.numeric(plyr::mapvalues(x = expt.obj@meta.data$subject,
                                                                      from = c('P5', 'P9', 'P10', 'P11', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'C4', 'C5', 'C6', 'C7', 'C13', 'C14'),
                                                                      to = c(8, 30, 12, 35, 9, 15, 12, 24, 7, 30, 18, 0, 0, 0, 0, 0, 0)))


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

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''), c(3, 2))

umap_orig_ident <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_orig_ident, paste('./plots/', experiment, '_prepare_umap_orig_ident', sep = ''), c(3, 2))

umap_subject <- DimPlot(expt.obj, reduction = "umap", group.by = "subject", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_subject, paste('./plots/', experiment, '_prepare_umap_subject', sep = ''), c(3, 2))

umap_PASI <- DimPlot(expt.obj, reduction = "umap", group.by = "PASI", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_PASI, paste('./plots/', experiment, '_prepare_umap_PASI', sep = ''), c(3, 2))

umap_PGA <- DimPlot(expt.obj, reduction = "umap", group.by = "PGA", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_PGA, paste('./plots/', experiment, '_prepare_umap_PGA', sep = ''), c(3, 2))

umap_Onset <- DimPlot(expt.obj, reduction = "umap", group.by = "Onset", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_Onset, paste('./plots/', experiment, '_prepare_umap_Onset', sep = ''), c(3, 2))

umap_Duration <- DimPlot(expt.obj, reduction = "umap", group.by = "Duration", shuffle = TRUE, seed = 123) + theme(plot.title = element_blank())
generate_figs(umap_Duration, paste('./plots/', experiment, '_prepare_umap_Duration', sep = ''), c(3, 2))

umap_Severity <- DimPlot(expt.obj, reduction = "umap", group.by = "Severity", shuffle = TRUE, seed = 123, cols = c('seagreen', 'steelblue', 'firebrick'), pt.size = 0.1) + 
  theme(plot.title = element_blank()) + scale_color_manual(labels=c("C", "M", "S"), values = c('seagreen', 'steelblue', 'firebrick'))
generate_figs(umap_Severity, paste('./plots/', experiment, '_prepare_umap_Severity', sep = ''), c(3, 2))


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




