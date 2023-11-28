# prepare seurat object
# Vandon Duong

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR)
library(scuttle)
library(data.table)
library(patchwork)

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

counts_data <- read.csv(file = "./data/GSE125881_raw.expMatrix.csv", header = TRUE, row.names = 1)
# expt.obj <- CreateSeuratObject(counts = counts_data, project = "Kinetics", min.cells = 3, min.features = 200)
# https://github.com/satijalab/seurat/issues/7545
expt.obj <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(counts_data), sparse = T), project = "Kinetics", min.cells = 3, min.features = 200)
class(expt.obj@assays$RNA$counts)

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

# Filter for CD8A +  cells, removing CD4+
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
for (i in colnames(expt.obj@meta.data)){
  j = levels(expt.obj@meta.data[[i]])
  if (is.null(j) == TRUE){
    print(paste(i, 'levels() is NULL'))
  } else {
    print(paste(i, 'levels() is not NULL'))
    print(j)
  }
  cat("\n")}

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'orig.ident', ncol=3)
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

saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

head(expt.obj)

# Generate UMAPs for metadata

umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_prepare_umap_seurat_clusters', sep = ''))

umap_seurat_clusters_highlight <- DimPlotHighlightIdents(expt.obj, seurat_clusters, 'umap', 'blue', 0.1, 4)
generate_figs(umap_seurat_clusters_highlight, paste('./plots/', experiment, '_prepare_umap_seurat_clusters_highlight', sep = ''), c(12, 10))

umap_patient <- DimPlot(expt.obj, reduction = "umap", group.by = "Patient", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_prepare_umap_patient', sep = ''))

umap_patient_highlight <- DimPlotHighlightIdents(expt.obj, Patient, 'umap', 'blue', 0.1, 2)
generate_figs(umap_patient_highlight, paste('./plots/', experiment, '_prepare_umap_patient_highlight', sep = ''), c(22, 20))

umap_timepoint <- DimPlot(expt.obj, reduction = "umap", group.by = "TimePoint", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint, paste('./plots/', experiment, '_prepare_umap_timepoint', sep = ''))

umap_timepoint_highlight <- DimPlotHighlightIdents(expt.obj, TimePoint, 'umap', 'blue', 0.1, 5)
generate_figs(umap_timepoint_highlight, paste('./plots/', experiment, '_prepare_umap_timepoint_highlight', sep = ''), c(20, 8))

umap_group <- DimPlot(expt.obj, reduction = "umap", group.by = "Group", shuffle = TRUE, seed = 123)
generate_figs(umap_group, paste('./plots/', experiment, '_prepare_umap_group', sep = ''))

umap_group_highlight <- DimPlotHighlightIdents(expt.obj, Group, 'umap', 'blue', 0.1, 2)
generate_figs(umap_group_highlight, paste('./plots/', experiment, '_prepare_umap_group_highlight', sep = ''), c(22, 20))

umap_disease <- DimPlot(expt.obj, reduction = "umap", group.by = "Disease", shuffle = TRUE, seed = 123)
generate_figs(umap_disease, paste('./plots/', experiment, '_prepare_umap_disease', sep = ''))

umap_disease_highlight <- DimPlotHighlightIdents(expt.obj, Disease, 'umap', 'blue', 0.1, 2)
generate_figs(umap_disease_highlight, paste('./plots/', experiment, '_prepare_umap_disease_highlight', sep = ''), c(15, 8))



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

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_barplot_dice_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('azimuth', 'Group')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_group <- ggplot(md_temp, aes(x = Group, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_group, paste('./plots/', experiment, '_barplot_azimuth_group', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'Group')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_group <- ggplot(md_temp, aes(x = Group, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_group, paste('./plots/', experiment, '_barplot_monaco_group', sep = ''))

md_temp <- md[, .N, by = c('dice', 'Group')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_group <- ggplot(md_temp, aes(x = Group, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_group, paste('./plots/', experiment, '_barplot_dice_group', sep = ''))

md_temp <- md[, .N, by = c('Phase', 'Group')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_phase_group <- ggplot(md_temp, aes(x = Group, y = N, fill = Phase)) + geom_col(position = "fill")
generate_figs(barplot_phase_group, paste('./plots/', experiment, '_barplot_phase_group', sep = ''))




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

expt.obj@meta.data$Activationi <- factor(expt.obj@meta.data$Activationi, levels = c(5,4,3,2,1,0,-1,-2,-3,-4,-5))
expt.obj@meta.data$Anergyi <- factor(expt.obj@meta.data$Anergyi, levels = c(5,4,3,2,1,0,-1,-2,-3,-4,-5))
expt.obj@meta.data$Stemnessi <- factor(expt.obj@meta.data$Stemnessi, levels = c(5,4,3,2,1,0,-1,-2,-3,-4,-5))
expt.obj@meta.data$Senescencei <- factor(expt.obj@meta.data$Senescencei, levels = c(5,4,3,2,1,0,-1,-2,-3,-4,-5))

expt.obj@meta.data$State1 <- NULL
expt.obj@meta.data$State2 <- NULL
expt.obj@meta.data$State3 <- NULL
expt.obj@meta.data$State4 <- NULL

saveRDS(expt.obj, file = paste('./data/', experiment, '_scored.rds', sep = ''))

head(expt.obj)














