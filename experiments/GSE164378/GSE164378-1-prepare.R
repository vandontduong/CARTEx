# prepare seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE164378'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR)
library(scuttle)
library(data.table)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

generate_figs = function(figure_object, file_name){
  ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
  ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  return (paste("generating figure for ", file_name))
}

integerize = function(score){
  score_mod = round(score)
  score_mod[score_mod < -4] <- -5
  score_mod[score_mod > 4] <- 5
  return (score_mod)
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# this is the downsampled version provided by SeuratData
# pbmcref <- readRDS(paste('./data/pbmcref.rds', sep = '')) 

pbmcref.matrix <- ReadMtx(mtx = './data/GSM5008737_RNA_3P-matrix.mtx', features = './data/GSM5008737_RNA_3P-features.tsv', cells = './data/GSM5008737_RNA_3P-barcodes.tsv')
pbmcref <- CreateSeuratObject(counts = pbmcref.matrix)
pbmcref.metadata <- read.table(file="./data/GSE164378_sc.meta.data_3P.csv", sep = ",", row.names = 1, header=T)

# https://github.com/emmanuelaaaaa/CyTOF_scRNA_integration/blob/58f6c371bdfaa2049f1abfe818c70fa53d4a79f2/Rscripts/data_preprocessing/PBMC_REF_Unvaccinated_RNA_batch_correction.R#L20
pbmcref <- AddMetaData(pbmcref, metadata = list(pbmcref.metadata$nCount_RNA, pbmcref.metadata$nFeature_RNA, pbmcref.metadata$orig.ident, pbmcref.metadata$lane, pbmcref.metadata$donor, pbmcref.metadata$time, pbmcref.metadata$celltype.l1, pbmcref.metadata$celltype.l2, pbmcref.metadata$celltype.l3, pbmcref.metadata$Phase, pbmcref.metadata$Batch), 
                       col.name = c("nCount_RNA", "nFeature_RNA", "orig.ident", "lane", "donor", "time", "celltype.l1", "celltype.l2", "celltype.l3", "Phase_ref", "Batch"))

pbmcref[["percent.mt"]] <- PercentageFeatureSet(pbmcref, pattern = "^MT-")

# Extract only the un-vaccinated (time 0) cells and doublets
pbmcref$donor_time <- paste(pbmcref$donor,"_",pbmcref$time,sep = "")
Idents(pbmcref) <- "donor_time"
pbmcref <- subset(x = pbmcref, idents = c("P1_0", "P2_0","P3_0","P4_0","P5_0","P6_0","P7_0","P8_0"))

Idents(pbmcref) <- "celltype.l2"
pbmcref <- subset(pbmcref, idents = "Doublet", invert=T)

# Quality filter
pbmcref <- subset(pbmcref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
pbmcref <- NormalizeData(pbmcref)

# extract genes
all.genes <- rownames(pbmcref)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = pbmcref, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'orig.ident', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality', sep = ''))

# Variable features and initial UMAP analysis
pbmcref <- FindVariableFeatures(pbmcref, selection.method = "vst", nfeatures = 2000)
top10_pbmcref <- head(VariableFeatures(pbmcref), 10)
varplt <- VariableFeaturePlot(pbmcref)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_pbmcref, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled', sep = ''))

pbmcref <- ScaleData(pbmcref, features = all.genes)
pbmcref <- RunPCA(pbmcref, features = VariableFeatures(object = pbmcref), npcs = 40)

pbmcref <- FindNeighbors(pbmcref, dims = 1:10)
pbmcref <- FindClusters(pbmcref, resolution = 0.5)
pbmcref <- RunUMAP(pbmcref, dims = 1:10)

saveRDS(pbmcref, file = paste('./data/', experiment, '.rds', sep = ''))

pbmcref <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

umap_seurat_clusters <- DimPlot(pbmcref, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))


####################################################################################################
######################################## Cell cycle analysis #######################################
####################################################################################################

pbmcref@meta.data$Phase_ref <- pbmcref@meta.data$Phase

# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "../cellannotate/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmcref <- CellCycleScoring(pbmcref, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(pbmcref[[]])

# Visualize the distribution of cell cycle markers across
ridgeplt <- RidgePlot(pbmcref, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
generate_figs(ridgeplt, paste('./plots/', experiment, '_ridgeplt', sep = ''))

pbmcref <- RunPCA(pbmcref, features = c(s.genes, g2m.genes))
# DimPlot(pbmcref)

# regress cell cycle
pbmcref <- ScaleData(pbmcref, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmcref))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
pbmcref <- RunPCA(pbmcref, features = VariableFeatures(pbmcref), nfeatures.print = 10)

saveRDS(pbmcref, file = paste('./data/', experiment, '_cellcycle.rds', sep = ''))

pbmcref <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))

# how many cells match (TRUE) / mismatch (FALSE)?
summary(pbmcref@meta.data$Phase == pbmcref@meta.data$Phase_ref)


####################################################################################################
######################################## Cell type annotation ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/introduction.html
# https://github.com/dviraran/SingleR/issues/150
# https://support.bioconductor.org/p/9136971/
# https://rdrr.io/github/LTLA/celldex/man/MonacoImmuneData.html
# https://rdrr.io/github/LTLA/celldex/man/DatabaseImmuneCellExpressionData.html

test_assay <- LayerData(pbmcref) # LayerData() is the updated function; GetAssayData() was depreciated

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html#subsetting-summarizedexperiment
monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
unique(monaco.obj$label.main)

dice.obj <- DatabaseImmuneCellExpressionData(ensembl=F)
dice.obj.md <- dice.obj@colData %>% as.data.table
dice.obj.md[, .N, by = c("label.main", "label.fine")]
# downsample
dice.obj@assays@data@listData$logcounts <- downsampleMatrix(dice.obj@assays@data@listData$logcounts, prop = 0.05)
unique(dice.obj$label.main)


monaco.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = monaco.obj, labels = monaco.obj$label.fine)
table(monaco.predictions$labels)
saveRDS(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))
write.csv(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.csv', sep = ''))
monaco.predictions <- readRDS(paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))

pbmcref[["monaco"]] <- monaco.predictions$labels

umap_predicted_monaco <- DimPlot(pbmcref, reduction = "umap", group.by = "monaco", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_monaco <- DimPlot(pbmcref, reduction = "umap", group.by = "monaco")
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_umap_predicted_monaco', sep = ''))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.csv', sep = ''))
dice.predictions <- readRDS(paste('./data/', experiment, '_dice_predictions.rds', sep = ''))

pbmcref[["dice"]] <- dice.predictions$labels

umap_predicted_dice <- DimPlot(pbmcref, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_dice <- DimPlot(pbmcref, reduction = "umap", group.by = "dice")
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_umap_predicted_dice', sep = ''))


featureplot_Tcell_markers <- FeaturePlot(pbmcref, features = c("CD4", "CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_featureplot_Tcell_markers', sep = ''))


pbmcref[["azimuth"]] <- pbmcref@meta.data$celltype.l2

# BAR CHARTS of CELL TYPES
md <- pbmcref@meta.data %>% as.data.table
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

md_temp <- md[, .N, by = c('Phase', 'azimuth')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_phase_azimuth <- ggplot(md_temp, aes(x = azimuth, y = N, fill = Phase)) + geom_col(position = "fill")
generate_figs(barplot_phase_azimuth, paste('./plots/', experiment, '_barplot_phase_azimuth', sep = ''))




# label CD8s based on having 
# azimuth: CD8 Naive; CD8 TCM; CD8 TEM; CD8 Proliferating
# monaco: Naive CD8 T cells; Effector memory CD8 T cells; Central memory CD8 T cells; Terminal effector CD8 T cells
# dice: T cells, CD8+, naive; T cells, CD8+, naive, stimulated"

# CD8 T cell detected in any reference
pbmcref$CD8Tref_1 <- with(pbmcref, ifelse((pbmcref@meta.data$azimuth == "CD8 Naive") | (pbmcref@meta.data$azimuth == "CD8 TEM") | (pbmcref@meta.data$azimuth == "CD8 Proliferating") | 
                                              (pbmcref@meta.data$monaco == "Naive CD8 T cells") | (pbmcref@meta.data$monaco == "Effector CD8 T cells") | (pbmcref@meta.data$monaco == "Central Memory CD8 T cells") | (pbmcref@meta.data$monaco == "Terminal effector CD8 T cells") |
                                              (pbmcref@meta.data$dice == "T cells, CD8+, naive") | (pbmcref@meta.data$dice == "T cells, CD8+, naive, stimulated") , 
                                            'CD8 T cell', 'Other'))

dplyr::count(pbmcref@meta.data, CD8Tref_1, sort = TRUE)

umap_predicted_CD8Tref_1 <- DimPlot(pbmcref, reduction = "umap", group.by = "CD8Tref_1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_CD8Tref_1, paste('./plots/', experiment, '_umap_predicted_CD8Tref_1', sep = ''))


# CD8 T cell detected in at least 2 references
pbmcref$CD8Tref_2 <- with(pbmcref, ifelse(pbmcref@meta.data$azimuth %in% c("CD8 Naive", "CD8 TEM", "CD8 Proliferating") +
                                              pbmcref@meta.data$monaco %in% c("Naive CD8 T cells", "Effector CD8 T cells", "Central Memory CD8 T cells", "Terminal effector CD8 T cells") +
                                              pbmcref@meta.data$dice %in% c("T cells, CD8+, naive", "T cells, CD8+, naive, stimulated") > 1,
                                            'CD8 T cell', 'Other'))

dplyr::count(pbmcref@meta.data, CD8Tref_2, sort = TRUE)

umap_predicted_CD8Tref_2 <- DimPlot(pbmcref, reduction = "umap", group.by = "CD8Tref_2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_CD8Tref_2, paste('./plots/', experiment, '_umap_predicted_CD8Tref_2', sep = ''))


# QC: Examine CARTEx representation at single-cell resolution

cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
all.genes <- rownames(pbmcref)

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes

CARTEx_630_cp <- PercentageFeatureSet(pbmcref, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_630_weights))) / length(rownames(cartex_630_weights))
pbmcref@meta.data$CARTEx_630_countsproportion <- CARTEx_630_cp[1:length(Cells(pbmcref))]
CARTEx_200_cp <- PercentageFeatureSet(pbmcref, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_200_weights))) / length(rownames(cartex_200_weights))
pbmcref@meta.data$CARTEx_200_countsproportion <- CARTEx_200_cp[1:length(Cells(pbmcref))]
CARTEx_84_cp <- PercentageFeatureSet(pbmcref, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_84_weights))) / length(rownames(cartex_84_weights))
pbmcref@meta.data$CARTEx_84_countsproportion <- CARTEx_84_cp[1:length(Cells(pbmcref))]



# pseudo-bulk reference
# https://bioconductor.org/books/release/SingleRBook/sc-mode.html#pseudo-bulk-aggregation



saveRDS(pbmcref, file = paste('./data/', experiment, '_annotated.rds', sep = ''))

pbmcref <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))



####################################################################################################
######################################## Filter on CD8 T cells ######################################
####################################################################################################

pbmcref <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))
Idents(pbmcref) <- "celltype.l2"
pbmcref.CD8pos <- subset(x = pbmcref, idents = c("CD8 TEM", "CD8 TCM", "CD8 Naive", "CD8 Proliferating"))


test_assay <- LayerData(pbmcref.CD8pos) # LayerData() is the updated function; GetAssayData() was depreciated

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html#subsetting-summarizedexperiment
monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
# monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells', 'CD4+ T cells', 'T cells')
monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells')
monaco.obj <- monaco.obj[, monaco.index]
unique(monaco.obj$label.main)

dice.obj <- DatabaseImmuneCellExpressionData(ensembl=F)
dice.obj.md <- dice.obj@colData %>% as.data.table
dice.obj.md[, .N, by = c("label.main", "label.fine")]
# dice.index <- dice.obj$label.main %in% c('T cells, CD4+', 'T cells, CD8+')
dice.index <- dice.obj$label.main %in% c('T cells, CD8+')
dice.obj <- dice.obj[, dice.index]
# downsample
dice.obj@assays@data@listData$logcounts <- downsampleMatrix(dice.obj@assays@data@listData$logcounts, prop = 0.05)
unique(dice.obj$label.main)

monaco.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = monaco.obj, labels = monaco.obj$label.fine)
table(monaco.predictions$labels)
saveRDS(monaco.predictions, file=paste('./data/', experiment, '_CD8pos_monaco_predictions.rds', sep = ''))
write.csv(monaco.predictions, file=paste('./data/', experiment, '_CD8pos_monaco_predictions.csv', sep = ''))
monaco.predictions <- readRDS(paste('./data/', experiment, '_CD8pos_monaco_predictions.rds', sep = ''))

pbmcref.CD8pos[["monaco"]] <- monaco.predictions$labels

umap_predicted_monaco <- DimPlot(pbmcref.CD8pos, reduction = "umap", group.by = "monaco", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_monaco <- DimPlot(pbmcref.CD8pos, reduction = "umap", group.by = "monaco")
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_CD8pos_umap_predicted_monaco', sep = ''))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_CD8pos_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_CD8pos_dice_predictions.csv', sep = ''))
dice.predictions <- readRDS(paste('./data/', experiment, '_CD8pos_dice_predictions.rds', sep = ''))

pbmcref.CD8pos[["dice"]] <- dice.predictions$labels

umap_predicted_dice <- DimPlot(pbmcref.CD8pos, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_dice <- DimPlot(pbmcref.CD8pos, reduction = "umap", group.by = "dice")
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_CD8pos_umap_predicted_dice', sep = ''))


featureplot_Tcell_markers <- FeaturePlot(pbmcref.CD8pos, features = c("CD4", "CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_CD8pos_featureplot_Tcell_markers', sep = ''))


pbmcref.CD8pos[["azimuth"]] <- pbmcref.CD8pos@meta.data$celltype.l2



# QC: Examine CARTEx representation at single-cell resolution

cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
all.genes <- rownames(pbmcref.CD8pos)

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes

CARTEx_630_cp <- PercentageFeatureSet(pbmcref.CD8pos, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_630_weights))) / length(rownames(cartex_630_weights))
pbmcref.CD8pos@meta.data$CARTEx_630_countsproportion <- CARTEx_630_cp[1:length(Cells(pbmcref.CD8pos))]
CARTEx_200_cp <- PercentageFeatureSet(pbmcref.CD8pos, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_200_weights))) / length(rownames(cartex_200_weights))
pbmcref.CD8pos@meta.data$CARTEx_200_countsproportion <- CARTEx_200_cp[1:length(Cells(pbmcref.CD8pos))]
CARTEx_84_cp <- PercentageFeatureSet(pbmcref.CD8pos, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_84_weights))) / length(rownames(cartex_84_weights))
pbmcref.CD8pos@meta.data$CARTEx_84_countsproportion <- CARTEx_84_cp[1:length(Cells(pbmcref.CD8pos))]



saveRDS(pbmcref.CD8pos, file = paste('./data/', experiment, '_CD8pos_annotated.rds', sep = ''))

pbmcref.CD8pos <- readRDS(paste('./data/', experiment, '_CD8pos_annotated.rds', sep = ''))



# BAR CHARTS of CELL TYPES
md <- pbmcref.CD8pos@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_CD8pos_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_CD8pos_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_CD8pos_barplot_dice_seurat_clusters', sep = ''))



####################################################################################################
######################################## Filter on all T cells ######################################
####################################################################################################

pbmcref <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))
Idents(pbmcref) <- "celltype.l1"
pbmcref.allT <- subset(x = pbmcref, idents = c("CD4 T", "CD8 T", "other T"))


test_assay <- LayerData(pbmcref.allT) # LayerData() is the updated function; GetAssayData() was depreciated

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html#subsetting-summarizedexperiment
monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells', 'CD4+ T cells', 'T cells')
monaco.obj <- monaco.obj[, monaco.index]
unique(monaco.obj$label.main)

dice.obj <- DatabaseImmuneCellExpressionData(ensembl=F)
dice.obj.md <- dice.obj@colData %>% as.data.table
dice.obj.md[, .N, by = c("label.main", "label.fine")]
dice.index <- dice.obj$label.main %in% c('T cells, CD4+', 'T cells, CD8+')
dice.obj <- dice.obj[, dice.index]
# downsample
dice.obj@assays@data@listData$logcounts <- downsampleMatrix(dice.obj@assays@data@listData$logcounts, prop = 0.05)
unique(dice.obj$label.main)

monaco.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = monaco.obj, labels = monaco.obj$label.fine)
table(monaco.predictions$labels)
saveRDS(monaco.predictions, file=paste('./data/', experiment, '_allT_monaco_predictions.rds', sep = ''))
write.csv(monaco.predictions, file=paste('./data/', experiment, '_allT_monaco_predictions.csv', sep = ''))
monaco.predictions <- readRDS(paste('./data/', experiment, '_allT_monaco_predictions.rds', sep = ''))

pbmcref.allT[["monaco"]] <- monaco.predictions$labels

umap_predicted_monaco <- DimPlot(pbmcref.allT, reduction = "umap", group.by = "monaco", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_monaco <- DimPlot(pbmcref.allT, reduction = "umap", group.by = "monaco")
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_allT_umap_predicted_monaco', sep = ''))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_allT_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_allT_dice_predictions.csv', sep = ''))
dice.predictions <- readRDS(paste('./data/', experiment, '_allT_dice_predictions.rds', sep = ''))

pbmcref.allT[["dice"]] <- dice.predictions$labels

umap_predicted_dice <- DimPlot(pbmcref.allT, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# umap_predicted_dice <- DimPlot(pbmcref.allT, reduction = "umap", group.by = "dice")
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_allT_umap_predicted_dice', sep = ''))


featureplot_Tcell_markers <- FeaturePlot(pbmcref.allT, features = c("CD4", "CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_allT_featureplot_Tcell_markers', sep = ''))


pbmcref.allT[["azimuth"]] <- pbmcref.allT@meta.data$celltype.l2



# QC: Examine CARTEx representation at single-cell resolution

cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
all.genes <- rownames(pbmcref.allT)

# Calculate the percentage of all counts that belong to a given set of features
# i.e. compute the percentage of transcripts that map to CARTEx genes

CARTEx_630_cp <- PercentageFeatureSet(pbmcref.allT, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_630_weights))) / length(rownames(cartex_630_weights))
pbmcref.allT@meta.data$CARTEx_630_countsproportion <- CARTEx_630_cp[1:length(Cells(pbmcref.allT))]
CARTEx_200_cp <- PercentageFeatureSet(pbmcref.allT, features = intersect(all.genes, rownames(cartex_200_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_200_weights))) / length(rownames(cartex_200_weights))
pbmcref.allT@meta.data$CARTEx_200_countsproportion <- CARTEx_200_cp[1:length(Cells(pbmcref.allT))]
CARTEx_84_cp <- PercentageFeatureSet(pbmcref.allT, features = intersect(all.genes, rownames(cartex_84_weights)), assay = 'RNA') * length(intersect(all.genes, rownames(cartex_84_weights))) / length(rownames(cartex_84_weights))
pbmcref.allT@meta.data$CARTEx_84_countsproportion <- CARTEx_84_cp[1:length(Cells(pbmcref.allT))]



saveRDS(pbmcref.allT, file = paste('./data/', experiment, '_allT_annotated.rds', sep = ''))

pbmcref.allT <- readRDS(paste('./data/', experiment, '_allT_annotated.rds', sep = ''))



# BAR CHARTS of CELL TYPES
md <- pbmcref.allT@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_allT_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_allT_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_allT_barplot_dice_seurat_clusters', sep = ''))
















####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

pbmcref.CD8pos <- readRDS(paste('./data/', experiment, '_CD8pos_annotated.rds', sep = ''))

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(pbmcref.CD8pos))
expr <- t(as.matrix(GetAssayData(pbmcref.CD8pos))[match(common, rownames(as.matrix(GetAssayData(pbmcref.CD8pos)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref.CD8pos@meta.data$CARTEx_630 <- Z(scores)
pbmcref.CD8pos@meta.data$CARTEx_630i <- integerize(pbmcref.CD8pos@meta.data$CARTEx_630)

pbmcref.CD8pos@meta.data$CARTEx_630_repcount <- rowSums(expr > 0)
pbmcref.CD8pos@meta.data$CARTEx_630_reppercent <- pbmcref.CD8pos@meta.data$CARTEx_630_repcount / 630

scatterplot_CARTEx_630_representation <- FeatureScatter(pbmcref.CD8pos, feature1 = 'CARTEx_630', feature2 = 'CARTEx_630_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_630_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_630_representation', sep = ''))

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(pbmcref.CD8pos))
expr <- t(as.matrix(GetAssayData(pbmcref.CD8pos))[match(common, rownames(as.matrix(GetAssayData(pbmcref.CD8pos)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref.CD8pos@meta.data$CARTEx_200 <- Z(scores)
pbmcref.CD8pos@meta.data$CARTEx_200i <- integerize(pbmcref.CD8pos@meta.data$CARTEx_200)

pbmcref.CD8pos@meta.data$CARTEx_200_repcount <- rowSums(expr > 0)
pbmcref.CD8pos@meta.data$CARTEx_200_reppercent <- pbmcref.CD8pos@meta.data$CARTEx_200_repcount / 200

scatterplot_CARTEx_200_representation <- FeatureScatter(pbmcref.CD8pos, feature1 = 'CARTEx_200', feature2 = 'CARTEx_200_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_200_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_200_representation', sep = ''))

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(pbmcref.CD8pos))
expr <- t(as.matrix(GetAssayData(pbmcref.CD8pos))[match(common, rownames(as.matrix(GetAssayData(pbmcref.CD8pos)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref.CD8pos@meta.data$CARTEx_84 <- Z(scores)
pbmcref.CD8pos@meta.data$CARTEx_84i <- integerize(pbmcref.CD8pos@meta.data$CARTEx_84)

pbmcref.CD8pos@meta.data$CARTEx_84_repcount <- rowSums(expr > 0)
pbmcref.CD8pos@meta.data$CARTEx_84_reppercent <- pbmcref.CD8pos@meta.data$CARTEx_84_repcount / 84

scatterplot_CARTEx_84_representation <- FeatureScatter(pbmcref.CD8pos, feature1 = 'CARTEx_84', feature2 = 'CARTEx_84_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_84_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_84_representation', sep = ''))

####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

pbmcref.CD8pos <- AddModuleScore(pbmcref.CD8pos, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
pbmcref.CD8pos@meta.data$Activation <- scale(pbmcref.CD8pos@meta.data$State1)
pbmcref.CD8pos@meta.data$Anergy <- scale(pbmcref.CD8pos@meta.data$State2)
pbmcref.CD8pos@meta.data$Stemness <- scale(pbmcref.CD8pos@meta.data$State3)
pbmcref.CD8pos@meta.data$Senescence <- scale(pbmcref.CD8pos@meta.data$State4)

pbmcref.CD8pos@meta.data$Activationi <- integerize(pbmcref.CD8pos@meta.data$Activation)
pbmcref.CD8pos@meta.data$Anergyi <- integerize(pbmcref.CD8pos@meta.data$Anergy)
pbmcref.CD8pos@meta.data$Stemnessi <- integerize(pbmcref.CD8pos@meta.data$Stemness)
pbmcref.CD8pos@meta.data$Senescencei <- integerize(pbmcref.CD8pos@meta.data$Senescence)

pbmcref.CD8pos@meta.data$State1 <- NULL
pbmcref.CD8pos@meta.data$State2 <- NULL
pbmcref.CD8pos@meta.data$State3 <- NULL
pbmcref.CD8pos@meta.data$State4 <- NULL

# examine other signatures
NK_like = rownames(read.csv("../../signatures/NK-like-dysfunction.csv", row.names=1, header = TRUE))
Wherry_Tex = rownames(read.csv("../../signatures/Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", row.names = 1, header = TRUE))
BBD_Tex = rownames(read.csv("../../signatures/Selli_2023_Blood_TBBDex.csv", row.names = 1, header = TRUE))
PD1_Tex = rownames(read.csv("../../signatures/Cai_2020_Pathology_PD1_Tex.csv", row.names = 1, header = TRUE))

pbmcref.CD8pos <- AddModuleScore(pbmcref.CD8pos, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

pbmcref.CD8pos@meta.data$NKlike_Tex <- scale(pbmcref.CD8pos@meta.data$Signature1)
pbmcref.CD8pos@meta.data$LCMV_Tex <- scale(pbmcref.CD8pos@meta.data$Signature2)
pbmcref.CD8pos@meta.data$BBD_Tex <- scale(pbmcref.CD8pos@meta.data$Signature3)
pbmcref.CD8pos@meta.data$PD1_Tex <- scale(pbmcref.CD8pos@meta.data$Signature4)

pbmcref.CD8pos@meta.data$NKlike_Texi <- integerize(pbmcref.CD8pos@meta.data$NKlike_Tex)
pbmcref.CD8pos@meta.data$LCMV_Texi <- integerize(pbmcref.CD8pos@meta.data$LCMV_Tex)
pbmcref.CD8pos@meta.data$BBD_Texi <- integerize(pbmcref.CD8pos@meta.data$BBD_Tex)
pbmcref.CD8pos@meta.data$PD1_Texi <- integerize(pbmcref.CD8pos@meta.data$PD1_Tex)

pbmcref.CD8pos@meta.data$Signature1 <- NULL
pbmcref.CD8pos@meta.data$Signature2 <- NULL
pbmcref.CD8pos@meta.data$Signature3 <- NULL
pbmcref.CD8pos@meta.data$Signature4 <- NULL

saveRDS(pbmcref.CD8pos, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(pbmcref.CD8pos)












# TransferIdent() is depreciated

data.to.transfer <- c("CARTEx_630", "CARTEx_630i", "CARTEx_630_repcount", "CARTEx_630_reppercent",
                      "CARTEx_200", "CARTEx_200i", "CARTEx_200_repcount", "CARTEx_200_reppercent",
                      "CARTEx_84", "CARTEx_84i", "CARTEx_84_repcount", "CARTEx_84_reppercent",
                      "Activation", "Anergy", "Stemness", "Senescence",
                      "Activationi", "Anergyi", "Stemnessi", "Senescencei",
                      "NKlike_Tex", "LCMV_Tex", "BBD_Tex", "PD1_Tex",
                      "NKlike_Texi", "LCMV_Texi", "BBD_Texi", "PD1_Texi")

for (metadatum in data.to.transfer) {
  pbmcref <- AddMetaData(object = pbmcref, metadata = pbmcref.CD8pos[[metadatum]], col.name = metadatum)
}

rm(pbmcref.CD8pos)

Idents(pbmcref) <- "celltype.l2"
saveRDS(pbmcref, file = paste('./data/', experiment, '_scored.rds', sep = ''))

head(pbmcref)





