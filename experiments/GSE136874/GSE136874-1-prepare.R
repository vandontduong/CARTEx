### Script name: GSE136874-1-prepare.R
### Description: prepare seurat object for GSE136874
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

CART_CD19.data <- Read10X_h5("./data/CD19-28z-filtered_gene_bc_matrices_h5.h5")
CART_GD2.data <- Read10X_h5("./data/GD2-28z-filtered_gene_bc_matrices_h5.h5")
CART_CD19 <- CreateSeuratObject(counts = CART_CD19.data, project = "CD19", min.cells = 3, min.features = 200)
CART_GD2 <- CreateSeuratObject(counts = CART_GD2.data, project = "GD2", min.cells = 3, min.features = 200)
expt.obj <- merge(CART_CD19, y = CART_GD2, add.cell.ids = c("CD19", "GD2"), project = "CART")
rm(list = c('CART_CD19.data', 'CART_GD2.data', 'CART_CD19', 'CART_GD2'))

expt.obj <- JoinLayers(expt.obj)
expt.obj <- NormalizeData(expt.obj)

# expt.obj <- SCTransform(expt.obj)
# expt.obj <- RunPCA(expt.obj, npcs = 30, verbose = F)
# expt.obj <- IntegrateLayers(expt.obj, method = RPCAIntegration, normalization.method = "SCT")

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

#- Add meta.data
expt.obj@meta.data$CAR <- expt.obj@meta.data$orig.ident
expt.obj@meta.data$CAR <- factor(expt.obj@meta.data$CAR, levels = c('CD19', 'GD2'))

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
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

head(expt.obj)


# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 3)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(10, 8))

umap_CAR <- DimPlot(expt.obj, reduction = "umap", group.by = "CAR", shuffle = TRUE, seed = 123)
generate_figs(umap_CAR, paste('./plots/', experiment, '_prepare_umap_CAR', sep = ''))

umap_CAR_highlight <- DimPlotHighlightIdents(expt.obj, CAR, 'umap', 'blue', 0.1, 2)
generate_figs(umap_CAR_highlight, paste('./plots/', experiment, '_prepare_umap_CAR_highlight', sep = ''), c(6,4))


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


####################################################################################################
######################################## Cell cycle analysis #######################################
####################################################################################################
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = paste(PATH_CELLANNOTATE, "nestorawa_forcellcycle_expressionMatrix.txt", sep = ''), header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

expt.obj <- CellCycleScoring(expt.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(expt.obj[[]])

# Visualize the distribution of cell cycle markers across
ridgeplt <- RidgePlot(expt.obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
generate_figs(ridgeplt, paste('./plots/', experiment, '_prepare_ridgeplt', sep = ''))

expt.obj <- RunPCA(expt.obj, features = c(s.genes, g2m.genes))
# DimPlot(expt.obj)

# regress cell cycle
expt.obj <- ScaleData(expt.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(expt.obj))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(expt.obj), nfeatures.print = 10)

expt.obj@meta.data$Phase <- factor(expt.obj@meta.data$Phase, levels = c('G1', 'S', 'G2M'))

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "Phase")]
phase_data <- md[, .N, by = c("Phase")]
setorder(phase_data, cols = "Phase")
phase_data$percent <- round(100*phase_data$N / sum(phase_data$N), digits = 1)
phase_data$label <- paste(phase_data$Phase, " (", phase_data$percent,"%)", sep = "")
phase_data$cols <- hcl.colors(n = nrow(phase_data), palette = "Temps")
barplot_phase <- ggplot(data= phase_data, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = phase_data$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_prepare_barplot_phase', sep = ''))

umap_phase <- DimPlot(expt.obj, group.by = "Phase", cols = phase_data$cols, shuffle = TRUE, seed = 123)
generate_figs(umap_phase, paste('./plots/', experiment, '_prepare_umap_phase', sep = ''))

umap_phase_highlight <- DimPlotHighlightIdents(expt.obj, Phase, 'umap', 'blue', 0.1, 3)
generate_figs(umap_phase_highlight, paste('./plots/', experiment, '_prepare_umap_phase_highlight', sep = ''), c(15, 6))


saveRDS(expt.obj, file = paste('./data/', experiment, '_cellcycle.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))

####################################################################################################
######################################## Cell type annotation ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/introduction.html
# https://github.com/dviraran/SingleR/issues/150
# https://support.bioconductor.org/p/9136971/
# https://rdrr.io/github/LTLA/celldex/man/MonacoImmuneData.html
# https://rdrr.io/github/LTLA/celldex/man/DatabaseImmuneCellExpressionData.html

test_assay <- LayerData(expt.obj) # LayerData() is the updated function; GetAssayData() was depreciated

# azimuth.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_cellcycle.rds")
azimuth.obj <- readRDS(paste(PATH_CELLANNOTATE, "AzimuthPBMC.rds", sep = ''))
azimuth.obj.md <- azimuth.obj@meta.data %>% as.data.table
azimuth.obj.md[, .N, by = c("celltype.l1", "celltype.l2")]
azimuth.obj <- SetIdent(azimuth.obj, value = "celltype.l1")
azimuth.obj <- subset(azimuth.obj, idents = c("CD8 T"))
# downsample
azimuth.obj <- SetIdent(azimuth.obj, value = "celltype.l2")
azimuth.obj <- subset(azimuth.obj, downsample = 100)
azimuth_assay <- LayerData(azimuth.obj)

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html#subsetting-summarizedexperiment
monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells')
monaco.obj <- monaco.obj[, monaco.index]
unique(monaco.obj$label.main)

dice.obj <- DatabaseImmuneCellExpressionData(ensembl=F)
dice.obj.md <- dice.obj@colData %>% as.data.table
dice.obj.md[, .N, by = c("label.main", "label.fine")]
dice.index <- dice.obj$label.main %in% c('T cells, CD8+')
dice.obj <- dice.obj[, dice.index]
# downsample
dice.obj@assays@data@listData$logcounts <- downsampleMatrix(dice.obj@assays@data@listData$logcounts, prop = 0.05)
unique(dice.obj$label.main)

azimuth.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = azimuth_assay, labels = azimuth.obj$celltype.l2)
table(azimuth.predictions$labels)
saveRDS(azimuth.predictions, paste(file='./data/', experiment, '_azimuth_predictions.rds', sep = ''))
write.csv(azimuth.predictions, paste(file='./data/', experiment, '_azimuth_predictions.csv', sep = ''))
# azimuth.predictions <- readRDS(paste('./data/', experiment, '_azimuth_predictions.rds', sep = ''))

expt.obj[["azimuth"]] <- azimuth.predictions$labels
unique(expt.obj[["azimuth"]])
expt.obj@meta.data$azimuth <- factor(expt.obj@meta.data$azimuth, levels = c('CD8 Naive', 'CD8 Proliferating'))

# umap_predicted_azimuth <- DimPlot(expt.obj, reduction = "umap", group.by = "azimuth", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_azimuth <- DimPlot(expt.obj, reduction = "umap", group.by = "azimuth")
generate_figs(umap_predicted_azimuth, paste('./plots/', experiment, '_umap_predicted_azimuth', sep = ''))

umap_predicted_azimuth_highlight <- DimPlotHighlightIdents(expt.obj, azimuth, 'umap', 'blue', 0.1, 2)
generate_figs(umap_predicted_azimuth_highlight, paste('./plots/', experiment, '_prepare_umap_predicted_azimuth_highlight', sep = ''), c(12, 8))


monaco.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = monaco.obj, labels = monaco.obj$label.fine)
table(monaco.predictions$labels)
saveRDS(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))
write.csv(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.csv', sep = ''))
# monaco.predictions <- readRDS(paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))

expt.obj[["monaco"]] <- monaco.predictions$labels
unique(expt.obj[["monaco"]])
expt.obj@meta.data$monaco <- factor(expt.obj@meta.data$monaco, levels = c('Naive CD8 T cells', 'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Terminal effector CD8 T cells'))

# umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco")
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_umap_predicted_monaco', sep = ''))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.csv', sep = ''))
# dice.predictions <- readRDS(paste('./data/', experiment, '_dice_predictions.rds', sep = ''))

expt.obj[["dice"]] <- dice.predictions$labels
unique(expt.obj[["dice"]])
expt.obj@meta.data$dice <- factor(expt.obj@meta.data$dice, levels = c('T cells, CD8+, naive', 'T cells, CD8+, naive, stimulated'))

# umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice")
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_umap_predicted_dice', sep = ''))

umap_predicted_dice_highlight <- DimPlotHighlightIdents(expt.obj, dice, 'umap', 'blue', 0.1, 2)
generate_figs(umap_predicted_dice_highlight, paste('./plots/', experiment, '_prepare_umap_predicted_dice_highlight', sep = ''), c(12, 8))


featureplot_Tcell_markers <- FeaturePlot(expt.obj, features = c("CD4", "CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_featureplot_Tcell_markers', sep = ''))


# BAR CHARTS of CELL TYPES
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

barplot_azimuth_seurat_clusters <- BarPlotStackSplit(expt.obj, 'azimuth', 'seurat_clusters')
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_azimuth_seurat_clusters', sep = ''), c(8,4))

barplot_monaco_seurat_clusters <- BarPlotStackSplit(expt.obj, 'monaco', 'seurat_clusters')
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_monaco_seurat_clusters', sep = ''), c(8,4))

barplot_dice_seurat_clusters <- BarPlotStackSplit(expt.obj, 'dice', 'seurat_clusters')
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_prepare_barplot_dice_seurat_clusters', sep = ''), c(8,4))

barplot_azimuth_CAR <- BarPlotStackSplit(expt.obj, 'azimuth', 'CAR')
generate_figs(barplot_azimuth_CAR, paste('./plots/', experiment, '_prepare_barplot_azimuth_CAR', sep = ''), c(8,4))

barplot_monaco_CAR <- BarPlotStackSplit(expt.obj, 'monaco', 'CAR')
generate_figs(barplot_monaco_CAR, paste('./plots/', experiment, '_prepare_barplot_monaco_CAR', sep = ''), c(8,4))

barplot_dice_CAR <- BarPlotStackSplit(expt.obj, 'dice', 'CAR')
generate_figs(barplot_dice_CAR, paste('./plots/', experiment, '_prepare_barplot_dice_CAR', sep = ''), c(8,4))

barplot_phase_CAR <- BarPlotStackSplit(expt.obj, 'Phase', 'CAR')
generate_figs(barplot_phase_CAR, paste('./plots/', experiment, '_prepare_barplot_phase_CAR', sep = ''), c(8,4))




saveRDS(expt.obj, file = paste('./data/', experiment, '_annotated.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))


####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

expt.obj <- AddModuleScore(expt.obj, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
expt.obj@meta.data$Activation <- scale(expt.obj@meta.data$State1)
expt.obj@meta.data$Anergy <- scale(expt.obj@meta.data$State2)
expt.obj@meta.data$Stemness <- scale(expt.obj@meta.data$State3)
expt.obj@meta.data$Senescence <- scale(expt.obj@meta.data$State4)

expt.obj@meta.data$Activationi <- integerize(expt.obj@meta.data$Activation)
expt.obj@meta.data$Anergyi <- integerize(expt.obj@meta.data$Anergy)
expt.obj@meta.data$Stemnessi <- integerize(expt.obj@meta.data$Stemness)
expt.obj@meta.data$Senescencei <- integerize(expt.obj@meta.data$Senescence)

expt.obj@meta.data$State1 <- NULL
expt.obj@meta.data$State2 <- NULL
expt.obj@meta.data$State3 <- NULL
expt.obj@meta.data$State4 <- NULL

# examine other signatures
NK_like = rownames(read.csv("../../signatures/NK-like-dysfunction.csv", row.names=1, header = TRUE))
Wherry_Tex = rownames(read.csv("../../signatures/Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", row.names = 1, header = TRUE))
BBD_Tex = rownames(read.csv("../../signatures/Selli_2023_Blood_TBBDex.csv", row.names = 1, header = TRUE))
PD1_Tex = rownames(read.csv("../../signatures/Cai_2020_Pathology_PD1_Tex.csv", row.names = 1, header = TRUE))

expt.obj <- AddModuleScore(expt.obj, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

expt.obj@meta.data$NKlike_Tex <- scale(expt.obj@meta.data$Signature1)
expt.obj@meta.data$LCMV_Tex <- scale(expt.obj@meta.data$Signature2)
expt.obj@meta.data$BBD_Tex <- scale(expt.obj@meta.data$Signature3)
expt.obj@meta.data$PD1_Tex <- scale(expt.obj@meta.data$Signature4)

expt.obj@meta.data$NKlike_Texi <- integerize(expt.obj@meta.data$NKlike_Tex)
expt.obj@meta.data$LCMV_Texi <- integerize(expt.obj@meta.data$LCMV_Tex)
expt.obj@meta.data$BBD_Texi <- integerize(expt.obj@meta.data$BBD_Tex)
expt.obj@meta.data$PD1_Texi <- integerize(expt.obj@meta.data$PD1_Tex)

expt.obj@meta.data$Signature1 <- NULL
expt.obj@meta.data$Signature2 <- NULL
expt.obj@meta.data$Signature3 <- NULL
expt.obj@meta.data$Signature4 <- NULL



saveRDS(expt.obj, file = paste('./data/', experiment, '_scored.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

head(expt.obj)


# UMAP of scores
umap_CARTEx_630i <- DimPlot(expt.obj, group.by = "CARTEx_630i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_630i, paste('./plots/', experiment, '_umap_CARTEx_630i', sep = ''))

umap_CARTEx_200i <- DimPlot(expt.obj, group.by = "CARTEx_200i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_200i, paste('./plots/', experiment, '_umap_CARTEx_200i', sep = ''))

umap_CARTEx_84i <- DimPlot(expt.obj, group.by = "CARTEx_84i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_84i, paste('./plots/', experiment, '_umap_CARTEx_84i', sep = ''))

umap_sig_activationi <- DimPlot(expt.obj, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_umap_sig_activationi', sep = ''))

umap_sig_anergyi <- DimPlot(expt.obj, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_umap_sig_anergyi', sep = ''))

umap_sig_stemnessi <- DimPlot(expt.obj, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_umap_sig_stemnessi', sep = ''))

umap_sig_senescencei <- DimPlot(expt.obj, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_umap_sig_senescencei', sep = ''))

# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(expt.obj, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(expt.obj, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(expt.obj, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(expt.obj, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(expt.obj, features = c("Senescence"), order = TRUE) + fix.sc

generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_umap_sig_activation', sep = ''))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_umap_sig_anergy', sep = ''))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_umap_sig_stemness', sep = ''))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_umap_sig_senescence', sep = ''))



####################################################################################################
####################################### CARTEx representation ######################################
####################################################################################################

cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
all.genes <- rownames(expt.obj)

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes

expt.obj@meta.data$CARTEx_630_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_630_weights))) / length(rownames(cartex_630_weights))
expt.obj@meta.data$CARTEx_200_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_200_weights))) / length(rownames(cartex_200_weights))
expt.obj@meta.data$CARTEx_84_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_84_weights))) / length(rownames(cartex_84_weights))

QCscatter_CARTEx_630_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_630", feature2 = "CARTEx_630_countsproportion")
QCscatter_CARTEx_200_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_200", feature2 = "CARTEx_200_countsproportion")
QCscatter_CARTEx_84_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_84", feature2 = "CARTEx_84_countsproportion")

generate_figs(QCscatter_CARTEx_630_countsproportion, paste('./plots/', experiment, '_QCscatter_CARTEx_630_countsproportion', sep = ''))
generate_figs(QCscatter_CARTEx_200_countsproportion, paste('./plots/', experiment, '_QCscatter_CARTEx_200_countsproportion', sep = ''))
generate_figs(QCscatter_CARTEx_84_countsproportion, paste('./plots/', experiment, '_QCscatter_CARTEx_84_countsproportion', sep = ''))

# do it for signatures as well

expt.obj@meta.data$Activation_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, activation_sig), assay = 'RNA') * length(intersect(all.genes, activation_sig)) / length(activation_sig)
expt.obj@meta.data$Anergy_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, anergy_sig), assay = 'RNA') * length(intersect(all.genes, anergy_sig)) / length(anergy_sig)
expt.obj@meta.data$Senescence_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, senescence_sig), assay = 'RNA') * length(intersect(all.genes, senescence_sig)) / length(senescence_sig)
expt.obj@meta.data$Stemness_countsproportion <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, stemness_sig), assay = 'RNA') * length(intersect(all.genes, stemness_sig)) / length(stemness_sig)


QCscatter_Activation_countsproportion <- FeatureScatter(expt.obj, feature1 = "Activation", feature2 = "Activation_countsproportion")
QCscatter_Anergy_countsproportion <- FeatureScatter(expt.obj, feature1 = "Anergy", feature2 = "Anergy_countsproportion")
QCscatter_Senescence_countsproportion <- FeatureScatter(expt.obj, feature1 = "Senescence", feature2 = "Senescence_countsproportion")
QCscatter_Stemness_countsproportion <- FeatureScatter(expt.obj, feature1 = "Stemness", feature2 = "Stemness_countsproportion")


generate_figs(QCscatter_Activation_countsproportion, paste('./plots/', experiment, '_QCscatter_Activation_countsproportion', sep = ''))
generate_figs(QCscatter_Anergy_countsproportion, paste('./plots/', experiment, '_QCscatter_Anergy_countsproportion', sep = ''))
generate_figs(QCscatter_Senescence_countsproportion, paste('./plots/', experiment, '_QCscatter_Senescence_countsproportion', sep = ''))
generate_figs(QCscatter_Stemness_countsproportion, paste('./plots/', experiment, '_QCscatter_Stemness_countsproportion', sep = ''))



####################################################################################################
############################################ CD8 subtypes ##########################################
####################################################################################################

# Paul, Michael St, and Pamela S. Ohashi. "The roles of CD8+ T cell subsets in antitumor immunity." Trends in cell biology 30.9 (2020): 695-704.
library(tensorflow)
library(cellassign)
# https://bioinformatics.stackexchange.com/questions/8830/import-gene-list-to-seurat-to-define-cell-types

# https://compbionju.github.io/scPlant/reference/AutoAnnotate_CellAssign.html

# https://www.rdocumentation.org/packages/cellassign/versions/0.99.2








