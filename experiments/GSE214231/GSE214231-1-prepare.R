#!/usr/bin/env Rscript

### Script name: GSE214231-1-prepare.R
### Description: prepare seurat object for GSE214231
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE214231'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################


### Make function to delete genes to avoid clustering biases
geneTCR <- function(Data = seurat){
  TCR.features <- c(grep("^TRAV", rownames(Data)),
                    grep("^TRAJ", rownames(Data)),
                    grep("^TRBV", rownames(Data)),
                    grep("^TRBD", rownames(Data)),
                    grep("^TRBJ", rownames(Data)),
                    grep("^XIST", rownames(Data))) ## also added this gene cause it was differentiation sex of donors
  return(TCR.features)
}


# 1.1: Loading of seurat objects and addition of metadata annotations 

### input: raw file data directory, name of sample, CAR annotation .csv
### output: list with seurat objects
Samples <- c("D1_2DomLib",
             "D1_2DomLib_WT",
             "D2_2DomLib_WT",
             "D3_2DomLib_1",
             "D3_2DomLib_2",
             "D3_2DomLib_WT")
Annotation <- c("GSM6601767_D1_2DomLib_barcodes_assigned.csv",
                "GSM6601768_D1_2DomLib_WT_barcodes_assigned.csv",
                "GSM6601769_D2_2DomLib_WT_barcodes_assigned.csv",
                "GSM6601770_D3_2DomLib_1_barcodes_assigned.csv",
                "GSM6601771_D3_2DomLib_2_barcodes_assigned.csv",
                "GSM6601772_D3_2DomLib_WT_barcodes_assigned.csv")
raw <- list()
for (j in 1:length(Samples)) {
  
  ## loading Count matrix
  # object <- Read10X(data.dir = Samples[j])
  object <- Read10X(data.dir = paste0("./data/", Samples[j]))
  
  ### Delete TCR genes to avoid clustering by clones
  features_to_delete <- geneTCR(Data = object)
  object <- object[-features_to_delete,]
  
  ## make seurat object
  object  <- CreateSeuratObject(object, project= Samples[j])
  ### assign CAR variants
  CAR_annotation <- read.csv(paste0("./data/", Annotation[j])) %>% select( -X)
  object  <- AddMetaData(
    object = object ,
    metadata = CAR_annotation$CAR,
    col.name = "CAR_Variant"
  )
  
  raw[[j]] <- object
  
}
names(raw) <- Samples 

# 1.2: Loading of control samples (incorporation into list of seurat objects in 1.1)
control_Samples <- c("D2_Unstim_28z",
                     "D3_hT")
Annotation <- c("unstim-CD28-CD3Z",
                "hT")
for (j in 1:length(control_Samples)) {
  
  ## loading Count matrix
  object <- Read10X(data.dir = paste0("./data/", control_Samples[j]))
  
  ### Delete TCR genes to avoid clustering by clones (function described above)
  features_to_delete <- geneTCR(Data = object)
  object <- object[-features_to_delete,]
  
  ## make seurat object
  object  <- CreateSeuratObject(object, project= control_Samples[j])
  object<- AddMetaData(
    object = object,
    metadata = Annotation[j],
    col.name = "CAR_Variant")
  
  raw[[j+ length(Samples)]] <- object
  
}
names(raw) <- c(Samples, control_Samples)

# 2.1: scRNAseq QC and data filtering 

### Merge seurat objects
expt.obj <- merge(raw[[1]], y = c(raw[c(2:length(raw))]), add.cell.ids = c(names(raw)))
#Plot number of cells per sample before QC
p1 <-  ggplot(expt.obj@meta.data, aes(orig.ident))+geom_bar(stat="count")+ theme(text = element_text(),axis.text.x = element_text(angle=45, hjust=1)) + geom_text(stat='count', aes(label=..count..), vjust=-1)

# join layers
expt.obj <- JoinLayers(expt.obj)

expt.obj <- NormalizeData(expt.obj)

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

expt.obj[['identifier']] <- experiment

# Name edits
expt.obj$CAR_Variant <-  gsub(' ','', expt.obj$CAR_Variant)
expt.obj$CAR_Variant <-  gsub('HVEM','HVEML', expt.obj$CAR_Variant)
expt.obj$CAR_Variant <-  gsub('4-1BB','41BB', expt.obj$CAR_Variant)
### Filter out cells that dont have a CAR assigned
Domain_A <- c("CD28", "41BB", "CD4", "CD150", "HVEML", "CD30", "CD84", "CD357", "CD244", "CD226", "FCRL6","FCRL1",
              "CD223", "DR3","TIM1")
Domain_B <- c("CD3Z",	"LMP2",	"K1",	"FCGR2A",	"FCER1G",	"DAP12",	"CD79B",	"CD79A",	"CD3G",	"CD3E",	"CD3D",	"GP")
variants <- paste(rep(Domain_A, each = length(Domain_B)), Domain_B, sep = "-")

## To include control groups
variants <- append(variants, "hT") 
variants <- append(variants, "unstim-CD28-CD3Z") 
expt.obj <- subset(expt.obj, cells = (rownames(expt.obj@meta.data[expt.obj@meta.data$CAR_Variant %in% variants,])))


### Filter out cells that dont express CD45 and Express HER2 (To get rid of SKBR3 cells)
# expt.obj <- subset(expt.obj, cells = WhichCells(expt.obj, expression = PTPRC > 0))
# expt.obj <- subset(expt.obj, cells = WhichCells(expt.obj, expression = ERBB2 < 2.5))
PTPRC_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["PTPRC",]
ERBB2_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["ERBB2",]
expt.obj <- subset(expt.obj, cells = names(which(PTPRC_expression > 0 & ERBB2_expression < 2.5)))


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

# meta data

table(expt.obj@meta.data$CAR_Variant)
# negative controls include:
# ht: co-cultured TCR-negative T cells (generated by Cas9-induced NHEJ in TRAC locus)
# unstim-CD28-CD3Z: no tumor co-culture of 28z CAR T cells

# Identify positions of unique cases
# unique_cases <- expt.obj@meta.data$CAR_Variant %in% c("hT", "unstim-CD28-CD3Z")
list_CAR_Variant <- plyr::mapvalues(x = expt.obj@meta.data$CAR_Variant,
                                            from = c("hT", "unstim-CD28-CD3Z"),
                                            to = c("hT-hT", "unstim-unstim"))

split_values <- strsplit(list_CAR_Variant, "-")

expt.obj@meta.data$CAR_Domain_A <- sapply(split_values, `[`, 1)
table(expt.obj@meta.data$CAR_Domain_A)
expt.obj@meta.data$CAR_Domain_A <- factor(expt.obj@meta.data$CAR_Domain_A, levels = c('41BB', 'CD28', 'CD4', 'CD150', 'CD226', 'HVEML', 'CD30', 'CD84', 'CD244', 'CD357', 'DR3', 'TIM1', 'FCRL6', 'FCRL1', 'CD223', 'hT', 'unstim'))

expt.obj@meta.data$CAR_Domain_B <- sapply(split_values, `[`, 2)
table(expt.obj@meta.data$CAR_Domain_B)
expt.obj@meta.data$CAR_Domain_B <- factor(expt.obj@meta.data$CAR_Domain_B, levels = c('CD3Z', 'LMP2', 'K1', 'FCGR2A', 'FCER1G', 'DAP12', 'CD79B', 'CD79A', 'CD3G', 'CD3E', 'CD3D', 'GP', 'hT', 'unstim'))

expt.obj@meta.data$CAR_Control <- rep("CAR", length(expt.obj@meta.data$CAR_Variant))
expt.obj@meta.data$CAR_Control[expt.obj@meta.data$CAR_Variant %in% c("hT")] <- "hT"
expt.obj@meta.data$CAR_Control[expt.obj@meta.data$CAR_Variant %in% c("unstim-CD28-CD3Z")] <- "unstim"
table(expt.obj@meta.data$CAR_Control)
expt.obj@meta.data$CAR_Control <- factor(expt.obj@meta.data$CAR_Control, levels = c('CAR', 'hT', 'unstim'))


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

umap_CAR_Variant <- DimPlot(expt.obj, reduction = "umap", group.by = "CAR_Variant", shuffle = TRUE, seed = 123)
generate_figs(umap_CAR_Variant, paste('./plots/', experiment, '_prepare_umap_CAR_Variant', sep = ''), c(6.5, 5))

umap_CAR_Domain_A <- DimPlot(expt.obj, reduction = "umap", group.by = "CAR_Domain_A", shuffle = TRUE, seed = 123)
generate_figs(umap_CAR_Domain_A, paste('./plots/', experiment, '_prepare_umap_CAR_Domain_A', sep = ''), c(6.5, 5))

umap_CAR_Domain_B <- DimPlot(expt.obj, reduction = "umap", group.by = "CAR_Domain_B", shuffle = TRUE, seed = 123)
generate_figs(umap_CAR_Domain_B, paste('./plots/', experiment, '_prepare_umap_CAR_Domain_B', sep = ''), c(6.5, 5))

umap_CAR_Control <- DimPlot(expt.obj, reduction = "umap", group.by = "CAR_Control", shuffle = TRUE, seed = 123)
generate_figs(umap_CAR_Control, paste('./plots/', experiment, '_prepare_umap_CAR_Control', sep = ''), c(6.5, 5))

umap_CAR_Control_highlight <- DimPlotHighlightIdents(expt.obj, CAR_Control, 'umap', 'blue', 0.1, 3)
generate_figs(umap_CAR_Control_highlight, paste('./plots/', experiment, '_prepare_umap_CAR_Control_highlight', sep = ''), c(12, 5))






