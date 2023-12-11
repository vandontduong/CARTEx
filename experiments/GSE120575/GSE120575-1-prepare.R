### Script name: GSE120575-1-prepare.R
### Description: prepare seurat object for GSE120575
### Author: Vandon Duong

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

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

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

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'characteristics_response', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_prepare_vlnplot_quality', sep = ''))

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

saveRDS(expt.obj, paste('./data/', experiment, '.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

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




####################################################################################################
######################################## Cell cycle analysis #######################################
####################################################################################################
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "../cellannotate/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)

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

expt.obj <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))


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
azimuth.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/cellannotate/AzimuthPBMC.rds")
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
umap_predicted_azimuth <- DimPlot(expt.obj, reduction = "umap", group.by = "azimuth", shuffle = TRUE, seed = 123)
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
umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123)
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_prepare_umap_predicted_monaco', sep = ''))

umap_predicted_monaco_highlight <- DimPlotHighlightIdents(expt.obj, monaco, 'umap', 'blue', 0.1, 2)
generate_figs(umap_predicted_monaco_highlight, paste('./plots/', experiment, '_prepare_umap_predicted_monaco_highlight', sep = ''), c(22, 20))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.csv', sep = ''))
# dice.predictions <- readRDS(paste('./data/', experiment, '_dice_predictions.rds', sep = ''))

expt.obj[["dice"]] <- dice.predictions$labels
unique(expt.obj[["dice"]])
expt.obj@meta.data$dice <- factor(expt.obj@meta.data$dice, levels = c('T cells, CD8+, naive', 'T cells, CD8+, naive, stimulated'))

# umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice", shuffle = TRUE, seed = 123)
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_umap_predicted_dice', sep = ''))

umap_predicted_dice_highlight <- DimPlotHighlightIdents(expt.obj, dice, 'umap', 'blue', 0.1, 2)
generate_figs(umap_predicted_dice_highlight, paste('./plots/', experiment, '_prepare_umap_predicted_dice_highlight', sep = ''), c(12, 8))


featureplot_Tcell_markers <- FeaturePlot(expt.obj, features = c("CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_prepare_featureplot_Tcell_markers', sep = ''))


# BAR CHARTS of CELL TYPES
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'characteristics_response')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_response <- ggplot(md_temp, aes(x = characteristics_response, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_response, paste('./plots/', experiment, '_prepare_barplot_azimuth_response', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'characteristics_response')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_response <- ggplot(md_temp, aes(x = characteristics_response, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_response, paste('./plots/', experiment, '_prepare_barplot_monaco_response', sep = ''))

md_temp <- md[, .N, by = c('dice', 'characteristics_response')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_response <- ggplot(md_temp, aes(x = characteristics_response, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_response, paste('./plots/', experiment, '_prepare_barplot_dice_response', sep = ''))

md_temp <- md[, .N, by = c('Phase', 'characteristics_response')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_phase_response <- ggplot(md_temp, aes(x = characteristics_response, y = N, fill = Phase)) + geom_col(position = "fill")
generate_figs(barplot_phase_response, paste('./plots/', experiment, '_prepare_barplot_phase_response', sep = ''))



# label CD8s based on having 
# azimuth: CD8 Naive; CD8 TCM; CD8 TEM; CD8 Proliferating
# monaco: Naive CD8 T cells; Effector memory CD8 T cells; Central memory CD8 T cells; Terminal effector CD8 T cells
# dice: T cells, CD8+, naive; T cells, CD8+, naive, stimulated"

# CD8 T cell detected in any reference
expt.obj$CD8Tref_1 <- with(expt.obj, ifelse((expt.obj@meta.data$azimuth == "CD8 Naive") | (expt.obj@meta.data$azimuth == "CD8 TEM") | (expt.obj@meta.data$azimuth == "CD8 Proliferating") | 
                                              (expt.obj@meta.data$monaco == "Naive CD8 T cells") | (expt.obj@meta.data$monaco == "Effector CD8 T cells") | (expt.obj@meta.data$monaco == "Central Memory CD8 T cells") | (expt.obj@meta.data$monaco == "Terminal effector CD8 T cells") |
                                              (expt.obj@meta.data$dice == "T cells, CD8+, naive") | (expt.obj@meta.data$dice == "T cells, CD8+, naive, stimulated") , 
                                            'CD8 T cell', 'Other'))

dplyr::count(expt.obj@meta.data, CD8Tref_1, sort = TRUE)

umap_predicted_CD8Tref_1 <- DimPlot(expt.obj, reduction = "umap", group.by = "CD8Tref_1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_CD8Tref_1, paste('./plots/', experiment, '_umap_predicted_CD8Tref_1', sep = ''))


# CD8 T cell detected in at least 2 references
expt.obj$CD8Tref_2 <- with(expt.obj, ifelse(expt.obj@meta.data$azimuth %in% c("CD8 Naive", "CD8 TEM", "CD8 Proliferating") +
                                              expt.obj@meta.data$monaco %in% c("Naive CD8 T cells", "Effector CD8 T cells", "Central Memory CD8 T cells", "Terminal effector CD8 T cells") +
                                              expt.obj@meta.data$dice %in% c("T cells, CD8+, naive", "T cells, CD8+, naive, stimulated") > 1,
                                            'CD8 T cell', 'Other'))

dplyr::count(expt.obj@meta.data, CD8Tref_2, sort = TRUE)

umap_predicted_CD8Tref_2 <- DimPlot(expt.obj, reduction = "umap", group.by = "CD8Tref_2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_CD8Tref_2, paste('./plots/', experiment, '_umap_predicted_CD8Tref_2', sep = ''))


# QC: Examine CARTEx representation at single-cell resolution

cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
all.genes <- rownames(expt.obj)

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes

CARTEx_630_cp <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_630_weights))) / length(rownames(cartex_630_weights))
expt.obj@meta.data$CARTEx_630_countsproportion <- CARTEx_630_cp[1:length(Cells(expt.obj))]
CARTEx_200_cp <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_200_weights))) / length(rownames(cartex_200_weights))
expt.obj@meta.data$CARTEx_200_countsproportion <- CARTEx_200_cp[1:length(Cells(expt.obj))]
CARTEx_84_cp <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_84_weights))) / length(rownames(cartex_84_weights))
expt.obj@meta.data$CARTEx_84_countsproportion <- CARTEx_84_cp[1:length(Cells(expt.obj))]




saveRDS(expt.obj, file = paste('./data/', experiment, '_annotated.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))





####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_630 <- Z(scores)
expt.obj@meta.data$CARTEx_630i <- integerize(expt.obj@meta.data$CARTEx_630)

expt.obj@meta.data$CARTEx_630_repcount <- rowSums(expr > 0)
expt.obj@meta.data$CARTEx_630_reppercent <- expt.obj@meta.data$CARTEx_630_repcount / 630

scatterplot_CARTEx_630_representation <- FeatureScatter(expt.obj, feature1 = 'CARTEx_630', feature2 = 'CARTEx_630_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_630_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_630_representation', sep = ''))

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_200 <- Z(scores)
expt.obj@meta.data$CARTEx_200i <- integerize(expt.obj@meta.data$CARTEx_200)

expt.obj@meta.data$CARTEx_200_repcount <- rowSums(expr > 0)
expt.obj@meta.data$CARTEx_200_reppercent <- expt.obj@meta.data$CARTEx_200_repcount / 200

scatterplot_CARTEx_200_representation <- FeatureScatter(expt.obj, feature1 = 'CARTEx_200', feature2 = 'CARTEx_200_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_200_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_200_representation', sep = ''))

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_84 <- Z(scores)
expt.obj@meta.data$CARTEx_84i <- integerize(expt.obj@meta.data$CARTEx_84)

expt.obj@meta.data$CARTEx_84_repcount <- rowSums(expr > 0)
expt.obj@meta.data$CARTEx_84_reppercent <- expt.obj@meta.data$CARTEx_84_repcount / 84

scatterplot_CARTEx_84_representation <- FeatureScatter(expt.obj, feature1 = 'CARTEx_84', feature2 = 'CARTEx_84_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_84_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_84_representation', sep = ''))

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


saveRDS(expt.obj, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(expt.obj)



















