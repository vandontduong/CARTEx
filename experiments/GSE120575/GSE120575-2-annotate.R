#!/usr/bin/env Rscript

### Script name: GSE120575-1-prepare.R
### Description: annotate seurat object for GSE120575
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE120575'
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
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "Phase")]
phase_data <- md[, .N, by = c("Phase")]
setorder(phase_data, cols = "Phase")
phase_data$percent <- round(100*phase_data$N / sum(phase_data$N), digits = 1)
phase_data$label <- paste(phase_data$Phase, " (", phase_data$percent,"%)", sep = "")
phase_data$cols <- hcl.colors(n = nrow(phase_data), palette = "Temps")
barplot_phase <- ggplot(data= phase_data, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = phase_data$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_prepare_barplot_phase', sep = ''))

umap_phase <- DimPlot(expt.obj, reduction = 'umap', group.by = "Phase", cols = phase_data$cols, shuffle = TRUE, seed = 123, pt.size = 0.1) + theme(plot.title = element_blank())
generate_figs(umap_phase, paste('./plots/', experiment, '_prepare_umap_phase', sep = ''), c(3, 2))

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
umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, cols = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'), pt.size = 0.1) + 
  theme(plot.title = element_blank()) + scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3'))
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_prepare_umap_predicted_monaco', sep = ''), c(2.9, 2))

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

barplot_azimuth_response <- BarPlotStackSplit(expt.obj, 'azimuth', 'characteristics_response')
generate_figs(barplot_azimuth_response, paste('./plots/', experiment, '_prepare_barplot_azimuth_response', sep = ''), c(8,4))

barplot_monaco_response <- BarPlotStackSplit(expt.obj, 'monaco', 'characteristics_response', color_set = c("deepskyblue", "seagreen", "darkgoldenrod", "plum3"))
generate_figs(barplot_monaco_response, paste('./plots/', experiment, '_prepare_barplot_monaco_response', sep = ''), c(6,4))

barplot_dice_response <- BarPlotStackSplit(expt.obj, 'dice', 'characteristics_response')
generate_figs(barplot_dice_response, paste('./plots/', experiment, '_prepare_barplot_dice_response', sep = ''), c(8,4))

barplot_phase_response <- BarPlotStackSplit(expt.obj, 'Phase', 'characteristics_response', color_set = hcl.colors(3, palette = "Temps"))
generate_figs(barplot_phase_response, paste('./plots/', experiment, '_prepare_barplot_phase_response', sep = ''), c(5,4))

barplot_monaco_response_slim <- barplot_monaco_response + theme(legend.position = "none")
generate_figs(barplot_monaco_response_slim, paste('./plots/', experiment, '_prepare_barplot_monaco_response_slim', sep = ''), c(1,1.75))

barplot_phase_response_slim <- barplot_phase_response + theme(legend.position = "none")
generate_figs(barplot_phase_response_slim, paste('./plots/', experiment, '_prepare_barplot_phase_response_slim', sep = ''), c(1,1.75))





saveRDS(expt.obj, file = paste('./data/', experiment, '_annotated.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))


####################################################################################################
########################################### Entropy scoring ########################################
####################################################################################################

# The ROGUE method was developed to robustly measure the purity of cell populations

# expt.obj <- UpdateSeuratObject(expt.obj)

# rogue.res <- EntropyScore(expt.obj, 'seurat_clusters', 'characteristics_response')
# rogue.res <- rogue.res[ , SortNumStrList(colnames(rogue.res), shift = FALSE)]
# rogue_boxplot_seurat_clusters_response <- rogue.boxplot(rogue.res)
# generate_figs(rogue_boxplot_seurat_clusters_response, paste('./plots/', experiment, '_prepare_rogue_boxplot_seurat_clusters_response', sep = ''), c(8,4))

# rogue.res <- EntropyScore(expt.obj, 'monaco', 'characteristics_response')
# rogue_boxplot_monaco_response <- rogue.boxplot(rogue.res)
# generate_figs(rogue_boxplot_monaco_response, paste('./plots/', experiment, '_prepare_rogue_boxplot_monaco_response', sep = ''), c(8,4))

####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

expt.obj <- ScoreSubroutine(expt.obj)



# UMAP of cell state scores
umap_sig_activationi <- DimPlot(expt.obj, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_prepare_umap_sig_activationi', sep = ''))

umap_sig_activationi_highlight <- DimPlotHighlightIdents(expt.obj, Activationi, 'umap', 'blue', 0.1, 4)
generate_figs(umap_sig_activationi_highlight, paste('./plots/', experiment, '_prepare_umap_sig_activationi_highlight', sep = ''), c(10, 8))

umap_sig_anergyi <- DimPlot(expt.obj, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_prepare_umap_sig_anergyi', sep = ''))

umap_sig_anergyi_highlight <- DimPlotHighlightIdents(expt.obj, Anergyi, 'umap', 'blue', 0.1, 4)
generate_figs(umap_sig_anergyi_highlight, paste('./plots/', experiment, '_prepare_umap_sig_anergyi_highlight', sep = ''), c(10, 8))

umap_sig_stemnessi <- DimPlot(expt.obj, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_prepare_umap_sig_stemnessi', sep = ''))

umap_sig_stemnessi_highlight <- DimPlotHighlightIdents(expt.obj, Stemnessi, 'umap', 'blue', 0.1, 4)
generate_figs(umap_sig_stemnessi_highlight, paste('./plots/', experiment, '_prepare_umap_sig_stemnessi_highlight', sep = ''), c(10, 8))

umap_sig_senescencei <- DimPlot(expt.obj, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_prepare_umap_sig_senescencei', sep = ''))

umap_sig_senescencei_highlight <- DimPlotHighlightIdents(expt.obj, Stemnessi, 'umap', 'blue', 0.1, 4)
generate_figs(umap_sig_senescencei_highlight, paste('./plots/', experiment, '_prepare_umap_sig_senescencei_highlight', sep = ''), c(10, 8))


# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(expt.obj, features = c("CARTEx_84"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_CARTEx_200 <- FeaturePlot(expt.obj, features = c("CARTEx_200"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_CARTEx_630 <- FeaturePlot(expt.obj, features = c("CARTEx_630"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_sig_activation <- FeaturePlot(expt.obj, features = c("Activation"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_sig_anergy <- FeaturePlot(expt.obj, features = c("Anergy"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_sig_stemness <- FeaturePlot(expt.obj, features = c("Stemness"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_sig_senescence <- FeaturePlot(expt.obj, features = c("Senescence"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_NKlike_Tex <- FeaturePlot(expt.obj, features = c("NKlike_Tex"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_LCMV_Tex <- FeaturePlot(expt.obj, features = c("LCMV_Tex"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_BBD_Tex <- FeaturePlot(expt.obj, features = c("BBD_Tex"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_PD1_Tex <- FeaturePlot(expt.obj, features = c("PD1_Tex"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())
umap_TSR <- FeaturePlot(expt.obj, features = c("TSR"), order = FALSE, pt.size = 0.1) + fix.sc + theme(plot.title = element_blank())

generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_prepare_umap_CARTEx_84', sep = ''), c(2.7, 2))
generate_figs(umap_CARTEx_200, paste('./plots/', experiment, '_prepare_umap_CARTEx_200', sep = ''), c(2.7, 2))
generate_figs(umap_CARTEx_630, paste('./plots/', experiment, '_prepare_umap_CARTEx_630', sep = ''), c(2.7, 2))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_prepare_umap_sig_activation', sep = ''), c(2.7, 2))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_prepare_umap_sig_anergy', sep = ''), c(2.7, 2))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_prepare_umap_sig_stemness', sep = ''), c(2.7, 2))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_prepare_umap_sig_senescence', sep = ''), c(2.7, 2))
generate_figs(umap_NKlike_Tex, paste('./plots/', experiment, '_prepare_umap_NKlike_Tex', sep = ''), c(2.7, 2))
generate_figs(umap_LCMV_Tex, paste('./plots/', experiment, '_prepare_umap_LCMV_Tex', sep = ''), c(2.7, 2))
generate_figs(umap_BBD_Tex, paste('./plots/', experiment, '_prepare_umap_BBD_Tex', sep = ''), c(2.7, 2))
generate_figs(umap_PD1_Tex, paste('./plots/', experiment, '_prepare_umap_PD1_Tex', sep = ''), c(2.7, 2))
generate_figs(umap_TSR, paste('./plots/', experiment, '_prepare_umap_TSR', sep = ''), c(2.7, 2))

saveRDS(expt.obj, file = paste('./data/', experiment, '_scored.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

head(expt.obj)




# report time
print("The script has completed...")
proc.time() - ptm

sessionInfo()




