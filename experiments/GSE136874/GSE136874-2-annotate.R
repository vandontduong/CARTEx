#!/usr/bin/env Rscript

### Script name: GSE136874-2-annotate.R
### Description: annotate seurat object for GSE136874
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Cell cycle analysis #######################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))

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
PCAplot <- DimPlot(expt.obj, reduction = 'pca')
generate_figs(PCAplot, paste('./plots/', experiment, '_prepare_PCAplot', sep = ''), c(6, 5))

# regress cell cycle
expt.obj <- ScaleData(expt.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(expt.obj))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(expt.obj), nfeatures.print = 10)
PCAplot_regressed <- DimPlot(expt.obj, reduction = 'pca')
generate_figs(PCAplot_regressed, paste('./plots/', experiment, '_prepare_PCAplot_regressed', sep = ''), c(6, 5))

expt.obj@meta.data$Phase <- factor(expt.obj@meta.data$Phase, levels = c('G1', 'S', 'G2M'))

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
# md <- expt.obj@meta.data %>% as.data.table
md <- setDT(expt.obj@meta.data, keep.rownames=TRUE)
md[, .N, by = c("orig.ident", "Phase")]
phase_data <- md[, .N, by = c("Phase")]
setorder(phase_data, cols = "Phase")
phase_data$percent <- round(100*phase_data$N / sum(phase_data$N), digits = 1)
phase_data$label <- paste(phase_data$Phase, " (", phase_data$percent,"%)", sep = "")
phase_data$cols <- hcl.colors(n = nrow(phase_data), palette = "Temps")
barplot_phase <- ggplot(data= phase_data, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = phase_data$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_prepare_barplot_phase', sep = ''))

umap_phase <- DimPlot(expt.obj, group.by = "Phase", cols = phase_data$cols, shuffle = TRUE, seed = 123)
generate_figs(umap_phase, paste('./plots/', experiment, '_prepare_umap_phase', sep = ''), c(6, 5))

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
# monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj <- readRDS(paste(PATH_CELLANNOTATE, "MonacoImmuneData.rds", sep = ''))
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells')
monaco.obj <- monaco.obj[, monaco.index]
unique(monaco.obj$label.main)

# dice.obj <- DatabaseImmuneCellExpressionData(ensembl=F)
dice.obj <- readRDS(paste(PATH_CELLANNOTATE, "DatabaseImmuneCellExpressionData.rds", sep = ''))
dice.obj.md <- dice.obj@colData %>% as.data.table
dice.obj.md[, .N, by = c("label.main", "label.fine")]
dice.index <- dice.obj$label.main %in% c('T cells, CD8+')
dice.obj <- dice.obj[, dice.index]
# downsample
dice.obj@assays@data@listData$logcounts <- downsampleMatrix(dice.obj@assays@data@listData$logcounts, prop = 0.05)
unique(dice.obj$label.main)

# SingleR() has this error
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
# solution is to downgrade matrixStats, according to https://github.com/satijalab/seurat/issues/7501#issuecomment-1854571904
# remotes::install_version("matrixStats", version="1.1.0")

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
generate_figs(umap_predicted_azimuth, paste('./plots/', experiment, '_umap_predicted_azimuth', sep = ''), c(6.5, 5))

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
umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, cols = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_prepare_umap_predicted_monaco', sep = ''), c(7.5, 5))

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
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_umap_predicted_dice', sep = ''), c(6.5, 5))

umap_predicted_dice_highlight <- DimPlotHighlightIdents(expt.obj, dice, 'umap', 'blue', 0.1, 2)
generate_figs(umap_predicted_dice_highlight, paste('./plots/', experiment, '_prepare_umap_predicted_dice_highlight', sep = ''), c(12, 8))


featureplot_Tcell_markers <- FeaturePlot(expt.obj, features = c("CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_prepare_featureplot_Tcell_markers', sep = ''))


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

barplot_monaco_CAR <- BarPlotStackSplit(expt.obj, 'monaco', 'CAR', color_set = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(barplot_monaco_CAR, paste('./plots/', experiment, '_prepare_barplot_monaco_CAR', sep = ''), c(6,4))

barplot_dice_CAR <- BarPlotStackSplit(expt.obj, 'dice', 'CAR')
generate_figs(barplot_dice_CAR, paste('./plots/', experiment, '_prepare_barplot_dice_CAR', sep = ''), c(8,4))

barplot_phase_CAR <- BarPlotStackSplit(expt.obj, 'Phase', 'CAR', color_set = hcl.colors(3, palette = "Temps"))
generate_figs(barplot_phase_CAR, paste('./plots/', experiment, '_prepare_barplot_phase_CAR', sep = ''), c(5,4))




saveRDS(expt.obj, file = paste('./data/', experiment, '_annotated.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))



####################################################################################################
########################################### Entropy scoring ########################################
####################################################################################################

# The ROGUE method was developed to robustly measure the purity of cell populations

rogue.res <- EntropyScore(expt.obj, 'seurat_clusters', 'CAR')
rogue_boxplot_seurat_clusters_CAR <- rogue.boxplot(rogue.res)
generate_figs(rogue_boxplot_seurat_clusters_CAR, paste('./plots/', experiment, '_prepare_rogue_boxplot_seurat_clusters_CAR', sep = ''), c(8,4))

rogue.res <- EntropyScore(expt.obj, 'monaco', 'CAR')
rogue_boxplot_monaco_CAR <- rogue.boxplot(rogue.res)
generate_figs(rogue_boxplot_monaco_CAR, paste('./plots/', experiment, '_prepare_rogue_boxplot_monaco_CAR', sep = ''), c(8,4))



####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj@meta.data$CARTEx_630 <- Z(SignatureScore(expt.obj, cartex_630_weights))
expt.obj@meta.data$CARTEx_630i <- integerize(expt.obj@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj@meta.data$CARTEx_200 <- Z(SignatureScore(expt.obj, cartex_200_weights))
expt.obj@meta.data$CARTEx_200i <- integerize(expt.obj@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj@meta.data$CARTEx_84 <- Z(SignatureScore(expt.obj, cartex_84_weights))
expt.obj@meta.data$CARTEx_84i <- integerize(expt.obj@meta.data$CARTEx_84)


####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = TRUE, row.names = 1))
senescence.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = TRUE, row.names = 1))

expt.obj <- AddModuleScore(expt.obj, features = list(activation.sig, anergy.sig, stemness.sig, senescence.sig), name="State", search = TRUE)

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
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))

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
umap_CARTEx_200 <- FeaturePlot(expt.obj, features = c("CARTEx_200"), order = TRUE) + fix.sc
umap_CARTEx_630 <- FeaturePlot(expt.obj, features = c("CARTEx_630"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(expt.obj, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(expt.obj, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(expt.obj, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(expt.obj, features = c("Senescence"), order = TRUE) + fix.sc
umap_NKlike_Tex <- FeaturePlot(expt.obj, features = c("NKlike_Tex"), order = TRUE) + fix.sc
umap_LCMV_Tex <- FeaturePlot(expt.obj, features = c("LCMV_Tex"), order = TRUE) + fix.sc
umap_BBD_Tex <- FeaturePlot(expt.obj, features = c("BBD_Tex"), order = TRUE) + fix.sc
umap_PD1_Tex <- FeaturePlot(expt.obj, features = c("PD1_Tex"), order = TRUE) + fix.sc

generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_prepare_umap_CARTEx_84', sep = ''), c(5.5, 5))
generate_figs(umap_CARTEx_200, paste('./plots/', experiment, '_prepare_umap_CARTEx_200', sep = ''), c(5.5, 5))
generate_figs(umap_CARTEx_630, paste('./plots/', experiment, '_prepare_umap_CARTEx_630', sep = ''), c(5.5, 5))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_prepare_umap_sig_activation', sep = ''), c(5.5, 5))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_prepare_umap_sig_anergy', sep = ''), c(5.5, 5))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_prepare_umap_sig_stemness', sep = ''), c(5.5, 5))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_prepare_umap_sig_senescence', sep = ''), c(5.5, 5))
generate_figs(umap_NKlike_Tex, paste('./plots/', experiment, '_prepare_umap_NKlike_Tex', sep = ''), c(5.5,5))
generate_figs(umap_LCMV_Tex, paste('./plots/', experiment, '_prepare_umap_LCMV_Tex', sep = ''), c(5.5,5))
generate_figs(umap_BBD_Tex, paste('./plots/', experiment, '_prepare_umap_BBD_Tex', sep = ''), c(5.5,5))
generate_figs(umap_PD1_Tex, paste('./plots/', experiment, '_prepare_umap_PD1_Tex', sep = ''), c(5.5,5))




dmap_CARTEx_84 <- FeaturePlot(expt.obj, reduction = 'dm', features = c("CARTEx_84"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_CARTEx_200 <- FeaturePlot(expt.obj, reduction = 'dm', features = c("CARTEx_200"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_CARTEx_630 <- FeaturePlot(expt.obj, reduction = 'dm', features = c("CARTEx_630"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_sig_activation <- FeaturePlot(expt.obj, reduction = 'dm', features = c("Activation"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_sig_anergy <- FeaturePlot(expt.obj, reduction = 'dm', features = c("Anergy"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_sig_stemness <- FeaturePlot(expt.obj, reduction = 'dm', features = c("Stemness"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))
dmap_sig_senescence <- FeaturePlot(expt.obj, reduction = 'dm', features = c("Senescence"), order = TRUE) + fix.sc + xlim(c(-0.2, 0.2)) + ylim(c(-0.15, 0.2))

generate_figs(dmap_CARTEx_84, paste('./plots/', experiment, '_prepare_dmap_CARTEx_84', sep = ''), c(5.5, 5))
generate_figs(dmap_CARTEx_200, paste('./plots/', experiment, '_prepare_dmap_CARTEx_200', sep = ''), c(5.5, 5))
generate_figs(dmap_CARTEx_630, paste('./plots/', experiment, '_prepare_dmap_CARTEx_630', sep = ''), c(5.5, 5))
generate_figs(dmap_sig_activation, paste('./plots/', experiment, '_prepare_dmap_sig_activation', sep = ''), c(5.5, 5))
generate_figs(dmap_sig_anergy, paste('./plots/', experiment, '_prepare_dmap_sig_anergy', sep = ''), c(5.5, 5))
generate_figs(dmap_sig_stemness, paste('./plots/', experiment, '_prepare_dmap_sig_stemness', sep = ''), c(5.5, 5))
generate_figs(dmap_sig_senescence, paste('./plots/', experiment, '_prepare_dmap_sig_senescence', sep = ''), c(5.5, 5))


saveRDS(expt.obj, file = paste('./data/', experiment, '_scored.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

head(expt.obj)

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# NT5E corresponds to CD73

vlnplot_CAR_exhaustion_markers <- VlnPlot(expt.obj, features = c('PDCD1', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'ENTPD1'), group.by = 'CAR', ncol = 3, cols = c('dodgerblue', 'indianred'), y.max = 3)

vlnplot_CAR_exhaustion_markers <- plot_grid(VlnPlot(expt.obj, features=c('PDCD1'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('HAVCR2'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('LAG3'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('CTLA4'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('TIGIT'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE),
                                            VlnPlot(expt.obj, features=c('ENTPD1'), group.by = 'CAR', cols = c('dodgerblue', 'indianred'), y.max = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + guides(fill=FALSE))


generate_figs(vlnplot_CAR_exhaustion_markers, paste('./plots/', experiment, '_prepare_vlnplot_CAR_exhaustion_markers', sep = ''), c(6,5))



# percentage of CARTEx detected

featplot_CARTEx_630_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_630', feature2 = 'CARTEx_630', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 630') + xlab('% detected of CARTEx 630') + xlim(c(0, 35)) + ylim(c(-3, 5))
featplot_CARTEx_200_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_200', feature2 = 'CARTEx_200', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 200') + xlab('% detected of CARTEx 200') + xlim(c(0, 35)) + ylim(c(-3, 5))
featplot_CARTEx_84_CAR <- FeatureScatter(expt.obj, feature1 = 'PFSD.CARTEx_84', feature2 = 'CARTEx_84', group.by = 'CAR', cols=c('dodgerblue', 'indianred'), shuffle = TRUE, seed = 123) + theme(legend.position = 'none') + ylab('CARTEx 84') + xlab('% detected of CARTEx 84') + xlim(c(0, 35)) + ylim(c(-3, 5))

featplot_CARTEx_combined_CAR <- (featplot_CARTEx_630_CAR | featplot_CARTEx_200_CAR | featplot_CARTEx_84_CAR)
generate_figs(featplot_CARTEx_combined_CAR, paste('./plots/', experiment, '_featplot_CARTEx_combined_CAR', sep = ''), c(10,5))






# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()



