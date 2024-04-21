#!/usr/bin/env Rscript

### Script name: GSE151511-1-prepare.R
### Description: prepare seurat object for GSE151511
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE151511'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# for 10X outputs, need to MANUALLY remove the prefix in the folder such that it has "barcodes.tsv" instead of "ac01_barcodes.tsv"
# uploaders changed the file names, so needs to be converted back

# read raw data, after manually renaming the files
pt01 <- Read10X("./data/filtered_gene_bc_matrices/ac01")
pt02 <- Read10X("./data/filtered_gene_bc_matrices/ac02")
pt03 <- Read10X("./data/filtered_gene_bc_matrices/ac03")
pt04 <- Read10X("./data/filtered_gene_bc_matrices/ac04")
pt05 <- Read10X("./data/filtered_gene_bc_matrices/ac05")
pt06 <- Read10X("./data/filtered_gene_bc_matrices/ac06")
pt07 <- Read10X("./data/filtered_gene_bc_matrices/ac07")
pt08 <- Read10X("./data/filtered_gene_bc_matrices/ac08")
pt09 <- Read10X("./data/filtered_gene_bc_matrices/ac09")
pt10 <- Read10X("./data/filtered_gene_bc_matrices/ac10")
pt11 <- Read10X("./data/filtered_gene_bc_matrices/ac11")
pt12 <- Read10X("./data/filtered_gene_bc_matrices/ac12")
pt13 <- Read10X("./data/filtered_gene_bc_matrices/ac13")
pt14 <- Read10X("./data/filtered_gene_bc_matrices/ac14")
pt15 <- Read10X("./data/filtered_gene_bc_matrices/ac15")
pt16 <- Read10X("./data/filtered_gene_bc_matrices/ac16")
pt17 <- Read10X("./data/filtered_gene_bc_matrices/ac17")
pt18 <- Read10X("./data/filtered_gene_bc_matrices/ac18")
pt19 <- Read10X("./data/filtered_gene_bc_matrices/ac19")
pt20 <- Read10X("./data/filtered_gene_bc_matrices/ac20")
pt21 <- Read10X("./data/filtered_gene_bc_matrices/ac21")
pt22 <- Read10X("./data/filtered_gene_bc_matrices/ac22")
pt23 <- Read10X("./data/filtered_gene_bc_matrices/ac23")
pt24 <- Read10X("./data/filtered_gene_bc_matrices/ac24")

genes_lists <- list(pt01 = unlist(pt01@Dimnames[1]), pt02 = unlist(pt02@Dimnames[1]), pt03 = unlist(pt03@Dimnames[1]), pt04 = unlist(pt04@Dimnames[1]),
                    pt05 = unlist(pt05@Dimnames[1]), pt06 = unlist(pt06@Dimnames[1]), pt07 = unlist(pt07@Dimnames[1]), pt08 = unlist(pt08@Dimnames[1]),
                    pt09 = unlist(pt09@Dimnames[1]), pt10 = unlist(pt10@Dimnames[1]), pt11 = unlist(pt11@Dimnames[1]), pt12 = unlist(pt12@Dimnames[1]),
                    pt13 = unlist(pt13@Dimnames[1]), pt14 = unlist(pt14@Dimnames[1]), pt15 = unlist(pt15@Dimnames[1]), pt16 = unlist(pt16@Dimnames[1]),
                    pt17 = unlist(pt17@Dimnames[1]), pt18 = unlist(pt18@Dimnames[1]), pt19 = unlist(pt19@Dimnames[1]), pt20 = unlist(pt20@Dimnames[1]),
                    pt21 = unlist(pt21@Dimnames[1]), pt22 = unlist(pt22@Dimnames[1]), pt23 = unlist(pt23@Dimnames[1]), pt24 = unlist(pt24@Dimnames[1]))


pt01_obj <- CreateSeuratObject(counts = pt01, project = "pt01")
pt02_obj <- CreateSeuratObject(counts = pt02, project = "pt02")
pt03_obj <- CreateSeuratObject(counts = pt03, project = "pt03")
pt04_obj <- CreateSeuratObject(counts = pt04, project = "pt04")
pt05_obj <- CreateSeuratObject(counts = pt05, project = "pt05")
pt06_obj <- CreateSeuratObject(counts = pt06, project = "pt06")
pt07_obj <- CreateSeuratObject(counts = pt07, project = "pt07")
pt08_obj <- CreateSeuratObject(counts = pt08, project = "pt08")
pt09_obj <- CreateSeuratObject(counts = pt09, project = "pt09")
pt10_obj <- CreateSeuratObject(counts = pt10, project = "pt10")
pt11_obj <- CreateSeuratObject(counts = pt11, project = "pt11")
pt12_obj <- CreateSeuratObject(counts = pt12, project = "pt12")
pt13_obj <- CreateSeuratObject(counts = pt13, project = "pt13")
pt14_obj <- CreateSeuratObject(counts = pt14, project = "pt14")
pt15_obj <- CreateSeuratObject(counts = pt15, project = "pt15")
pt16_obj <- CreateSeuratObject(counts = pt16, project = "pt16")
pt17_obj <- CreateSeuratObject(counts = pt17, project = "pt17")
pt18_obj <- CreateSeuratObject(counts = pt18, project = "pt18")
pt19_obj <- CreateSeuratObject(counts = pt19, project = "pt19")
pt20_obj <- CreateSeuratObject(counts = pt20, project = "pt20")
pt21_obj <- CreateSeuratObject(counts = pt21, project = "pt21")
pt22_obj <- CreateSeuratObject(counts = pt22, project = "pt22")
pt23_obj <- CreateSeuratObject(counts = pt23, project = "pt23")
pt24_obj <- CreateSeuratObject(counts = pt24, project = "pt24")

expt.obj <- merge(pt01_obj, y = c(pt02_obj, pt03_obj, pt04_obj, pt05_obj, pt06_obj, pt07_obj, pt08_obj, pt09_obj,
                                       pt10_obj, pt11_obj, pt12_obj, pt13_obj, pt14_obj, pt15_obj, pt16_obj, pt17_obj,
                                       pt18_obj, pt19_obj, pt20_obj, pt21_obj, pt22_obj, pt23_obj, pt24_obj), 
                       add.cell.ids = c("pt01", "pt02", "pt03", "pt04", "pt05", "pt06", "pt07", "pt08", "pt09", "pt10", "pt11", "pt12", "pt13",
                                        "pt14", "pt15", "pt16", "pt17", "pt18", "pt19", "pt20", "pt21", "pt22", "pt23", "pt24"), project = "CART")

rm(list = c('pt01_obj', 'pt02_obj', 'pt03_obj', 'pt04_obj', 'pt05_obj', 'pt06_obj', 'pt07_obj', 'pt08_obj', 'pt09_obj',
            'pt10_obj', 'pt11_obj', 'pt12_obj', 'pt13_obj', 'pt14_obj', 'pt15_obj', 'pt16_obj', 'pt17_obj',
            'pt18_obj', 'pt19_obj', 'pt20_obj', 'pt21_obj', 'pt22_obj', 'pt23_obj', 'pt24_obj'))

rm(list = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09',
            'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17',
            'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'))

expt.obj <- JoinLayers(expt.obj)
expt.obj <- NormalizeData(expt.obj)
expt.obj <- UpdateSeuratObject(expt.obj)

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

#- Add meta.data (should have added earlier)
expt.obj@meta.data$ClinicalResponse <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                                       from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                       to = c('CR', 'PD', 'PD', 'PD', 'CR', 'NE', 'CR', 'CR', 'CR', 'CR', 'PD', 'CR', 'PD', 'CR', 'PD', 'CR', 'PD', 'PD', 'PD', 'PR', 'PD', 'PD', 'PD', 'PD'))
expt.obj@meta.data$ClinicalResponse <- factor(expt.obj@meta.data$ClinicalResponse, levels = c('CR', 'PR', 'PD', 'NE'))
expt.obj@meta.data$EMR <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                          from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                          to = c('GOOD', 'NE', 'BAD', 'NE', 'GOOD', 'BAD', 'GOOD', 'GOOD', 'NE', 'GOOD', 'BAD', 'NE', 'BAD', 'NE', 'NE', 'GOOD', 'GOOD', 'BAD', 'BAD', 'BAD', 'GOOD', 'BAD', 'BAD', 'NE'))
expt.obj@meta.data$EMR <- factor(expt.obj@meta.data$EMR, levels = c('GOOD', 'BAD', 'NE'))

expt.obj@meta.data$Responder <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                                from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                to = c('R', 'NR', 'NR', 'NR', 'R', 'NE', 'R', 'R', 'R', 'R', 'NR', 'R', 'NR', 'R', 'NR', 'R', 'NR', 'NR', 'NR', 'R', 'NR', 'NR', 'NR', 'NR'))
expt.obj@meta.data$Responder <- factor(expt.obj@meta.data$Responder, levels = c('R', 'NR', 'NE'))

expt.obj@meta.data$CRS <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                          from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                          to = c(2, 2, 4, 3, 1, 2, 1, 2, 2, 2, 2, 3, 1, 2, 2, 2, 1, 1, 1, 1, 3, 1, 2, 2))
expt.obj@meta.data$CRS <- factor(expt.obj@meta.data$CRS, levels = c(1, 2, 3, 4))

expt.obj@meta.data$ICANS <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                            from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                            to = c(3, 0, 4, 3, 0, 4, 3, 0, 3, 0, 3, 3, 0, 2, 3, 1, 0, 2, 2, 0, 3, 3, 3, 0))
expt.obj@meta.data$ICANS <- factor(expt.obj@meta.data$ICANS, levels = c(0, 1, 2, 3, 4))

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

# print("Calculating diffusion maps...")
# diffusion_map <- RunDiffusion(expt.obj, k_int = 10)
# expt.obj <- diffusion_map$atlas
# saveRDS(diffusion_map$dmap, file = paste('./data/', experiment, '_dmap.rds', sep = ''))
# rm(diffusion_map)

print("Saving Seurat object...")
saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

head(expt.obj)


# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''), c(6, 5))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 10))

umap_patient <- DimPlot(expt.obj, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_prepare_umap_patient', sep = ''), c(6.5,5))

umap_clinicalresponse <- DimPlot(expt.obj, reduction = "umap", group.by = "ClinicalResponse", shuffle = TRUE, seed = 123, cols = c("seagreen", "slateblue", "firebrick", "lightgrey"))
generate_figs(umap_clinicalresponse, paste('./plots/', experiment, '_prepare_umap_clinicalresponse', sep = ''), c(6,5))

umap_responder <- DimPlot(expt.obj, reduction = "umap", group.by = "Responder", shuffle = TRUE, seed = 123, cols = c("seagreen", "firebrick", "lightgrey"))
generate_figs(umap_responder, paste('./plots/', experiment, '_prepare_umap_responder', sep = ''), c(6.5,5))

umap_EMR <- DimPlot(expt.obj, reduction = "umap", group.by = "EMR", shuffle = TRUE, seed = 123, cols = c("seagreen", "firebrick", "lightgrey"))
generate_figs(umap_EMR, paste('./plots/', experiment, '_prepare_umap_EMR', sep = ''), c(6.5,5))



####################################################################################################
###################################### Seurat cluster analysis #####################################
####################################################################################################

# examine metadata split by Seurat clusters

barplot_patient_seurat_clusters <- BarPlotStackSplit(expt.obj, 'orig.ident', 'seurat_clusters')
generate_figs(barplot_patient_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_patient_seurat_clusters', sep = ''), c(8,4))

barplot_clinical_response_seurat_clusters <- BarPlotStackSplit(expt.obj, 'ClinicalResponse', 'seurat_clusters')
generate_figs(barplot_clinical_response_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_clinical_response_seurat_clusters', sep = ''), c(8,4))

barplot_EMR_seurat_clusters <- BarPlotStackSplit(expt.obj, 'EMR', 'seurat_clusters')
generate_figs(barplot_EMR_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_EMR_seurat_clusters', sep = ''), c(8,4))

barplot_CRS_seurat_clusters <- BarPlotStackSplit(expt.obj, 'CRS', 'seurat_clusters')
generate_figs(barplot_CRS_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_CRS_seurat_clusters', sep = ''), c(8,4))

barplot_ICANS_seurat_clusters <- BarPlotStackSplit(expt.obj, 'ICANS', 'seurat_clusters')
generate_figs(barplot_ICANS_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_ICANS_seurat_clusters', sep = ''), c(8,4))


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



















