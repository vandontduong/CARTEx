# prepare seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE151511'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR) # BiocManager::install("beachmat", force = TRUE)
library(beachmat)
# library(celldex)
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


pt01_obj <- CreateSeuratObject(counts = pt01, min.cells = 3, min.genes = 200, project = "pt01")
pt02_obj <- CreateSeuratObject(counts = pt02, min.cells = 3, min.genes = 200, project = "pt02")
pt03_obj <- CreateSeuratObject(counts = pt03, min.cells = 3, min.genes = 200, project = "pt03")
pt04_obj <- CreateSeuratObject(counts = pt04, min.cells = 3, min.genes = 200, project = "pt04")
pt05_obj <- CreateSeuratObject(counts = pt05, min.cells = 3, min.genes = 200, project = "pt05")
pt06_obj <- CreateSeuratObject(counts = pt06, min.cells = 3, min.genes = 200, project = "pt06")
pt07_obj <- CreateSeuratObject(counts = pt07, min.cells = 3, min.genes = 200, project = "pt07")
pt08_obj <- CreateSeuratObject(counts = pt08, min.cells = 3, min.genes = 200, project = "pt08")
pt09_obj <- CreateSeuratObject(counts = pt09, min.cells = 3, min.genes = 200, project = "pt09")
pt10_obj <- CreateSeuratObject(counts = pt10, min.cells = 3, min.genes = 200, project = "pt10")
pt11_obj <- CreateSeuratObject(counts = pt11, min.cells = 3, min.genes = 200, project = "pt11")
pt12_obj <- CreateSeuratObject(counts = pt12, min.cells = 3, min.genes = 200, project = "pt12")
pt13_obj <- CreateSeuratObject(counts = pt13, min.cells = 3, min.genes = 200, project = "pt13")
pt14_obj <- CreateSeuratObject(counts = pt14, min.cells = 3, min.genes = 200, project = "pt14")
pt15_obj <- CreateSeuratObject(counts = pt15, min.cells = 3, min.genes = 200, project = "pt15")
pt16_obj <- CreateSeuratObject(counts = pt16, min.cells = 3, min.genes = 200, project = "pt16")
pt17_obj <- CreateSeuratObject(counts = pt17, min.cells = 3, min.genes = 200, project = "pt17")
pt18_obj <- CreateSeuratObject(counts = pt18, min.cells = 3, min.genes = 200, project = "pt18")
pt19_obj <- CreateSeuratObject(counts = pt19, min.cells = 3, min.genes = 200, project = "pt19")
pt20_obj <- CreateSeuratObject(counts = pt20, min.cells = 3, min.genes = 200, project = "pt20")
pt21_obj <- CreateSeuratObject(counts = pt21, min.cells = 3, min.genes = 200, project = "pt21")
pt22_obj <- CreateSeuratObject(counts = pt22, min.cells = 3, min.genes = 200, project = "pt22")
pt23_obj <- CreateSeuratObject(counts = pt23, min.cells = 3, min.genes = 200, project = "pt23")
pt24_obj <- CreateSeuratObject(counts = pt24, min.cells = 3, min.genes = 200, project = "pt24")

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

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")
head(expt.obj)

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
                                                to = c('Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Responder', 'NE', 'Responder', 'Responder', 'Responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Non-responder'))
expt.obj@meta.data$Responder <- factor(expt.obj@meta.data$Responder, levels = c('Responder', 'Non-responder', 'NE'))

expt.obj@meta.data$CRS <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                          from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                          to = c(2, 2, 4, 3, 1, 2, 1, 2, 2, 2, 2, 3, 1, 2, 2, 2, 1, 1, 1, 1, 3, 1, 2, 2))
expt.obj@meta.data$CRS <- factor(expt.obj@meta.data$CRS, levels = c(1, 2, 3, 4))

expt.obj@meta.data$ICANS <- plyr::mapvalues(x = expt.obj@meta.data$orig.ident,
                                            from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                            to = c(3, 0, 4, 3, 0, 4, 3, 0, 3, 0, 3, 3, 0, 2, 3, 1, 0, 2, 2, 0, 3, 3, 3, 0))
expt.obj@meta.data$ICANS <- factor(expt.obj@meta.data$ICANS, levels = c(0, 1, 2, 3, 4))



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
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'orig.ident', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

saveRDS(expt.obj, file = paste('./data/', experiment, '.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))


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
generate_figs(ridgeplt, paste('./plots/', experiment, '_ridgeplt', sep = ''))

expt.obj <- RunPCA(expt.obj, features = c(s.genes, g2m.genes))
# DimPlot(expt.obj)

# regress cell cycle
expt.obj <- ScaleData(expt.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(expt.obj))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(expt.obj), nfeatures.print = 10)

# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "Phase")]
phase_data <- md[, .N, by = c("Phase")]
setorder(phase_data, cols = "Phase")
phase_data$percent <- round(100*phase_data$N / sum(phase_data$N), digits = 1)
phase_data$label <- paste(phase_data$Phase, " (", phase_data$percent,"%)", sep = "")
phase_data$cols <- hcl.colors(n = nrow(phase_data), palette = "Temps")
barplot_phase <- ggplot(data= phase_data, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = phase_data$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_barplot_phase', sep = ''))
umap_phase <- DimPlot(expt.obj, group.by = "Phase", cols = phase_data$cols)
generate_figs(umap_phase, paste('./plots/', experiment, '_umap_phase', sep = ''))


saveRDS(expt.obj, file = paste('./data/', experiment, '_cellcycle.rds', sep = ''))

rm(CART.combined)
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
azimuth.predictions <- readRDS(paste('./data/', experiment, '_azimuth_predictions.rds', sep = ''))

expt.obj[["azimuth"]] <- azimuth.predictions$labels

# umap_predicted_azimuth <- DimPlot(expt.obj, reduction = "umap", group.by = "azimuth", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_azimuth <- DimPlot(expt.obj, reduction = "umap", group.by = "azimuth")
generate_figs(umap_predicted_azimuth, paste('./plots/', experiment, '_umap_predicted_azimuth', sep = ''))


monaco.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = monaco.obj, labels = monaco.obj$label.fine)
table(monaco.predictions$labels)
saveRDS(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))
write.csv(monaco.predictions, file=paste('./data/', experiment, '_monaco_predictions.csv', sep = ''))
monaco.predictions <- readRDS(paste('./data/', experiment, '_monaco_predictions.rds', sep = ''))

expt.obj[["monaco"]] <- monaco.predictions$labels

# umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_monaco <- DimPlot(expt.obj, reduction = "umap", group.by = "monaco")
generate_figs(umap_predicted_monaco, paste('./plots/', experiment, '_umap_predicted_monaco', sep = ''))


dice.predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = dice.obj, labels = dice.obj$label.fine)
table(dice.predictions$labels)
saveRDS(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.rds', sep = ''))
write.csv(dice.predictions, file=paste('./data/', experiment, '_dice_predictions.csv', sep = ''))
dice.predictions <- readRDS(paste('./data/', experiment, '_dice_predictions.rds', sep = ''))

expt.obj[["dice"]] <- dice.predictions$labels

# umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_dice <- DimPlot(expt.obj, reduction = "umap", group.by = "dice")
generate_figs(umap_predicted_dice, paste('./plots/', experiment, '_umap_predicted_dice', sep = ''))


featureplot_Tcell_markers <- FeaturePlot(expt.obj, features = c("CD4", "CD8A", "CD8B", "PDCD1"))
generate_figs(featureplot_Tcell_markers, paste('./plots/', experiment, '_featureplot_Tcell_markers', sep = ''))


# BAR CHARTS of CELL TYPES
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'orig.ident')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_CAR <- ggplot(md_temp, aes(x = orig.ident, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_CAR, paste('./plots/', experiment, '_barplot_azimuth_CAR', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'orig.ident')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_CAR <- ggplot(md_temp, aes(x = orig.ident, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_CAR, paste('./plots/', experiment, '_barplot_monaco_CAR', sep = ''))

md_temp <- md[, .N, by = c('dice', 'orig.ident')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_CAR <- ggplot(md_temp, aes(x = orig.ident, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_CAR, paste('./plots/', experiment, '_barplot_dice_CAR', sep = ''))

md_temp <- md[, .N, by = c('Phase', 'orig.ident')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_phase_CAR <- ggplot(md_temp, aes(x = orig.ident, y = N, fill = Phase)) + geom_col(position = "fill")
generate_figs(barplot_phase_CAR, paste('./plots/', experiment, '_barplot_phase_CAR', sep = ''))




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







# Novelty score
expt.obj$log10GenesPerUMI <- log10(expt.obj$nFeature_RNA) / log10(expt.obj$nCount_RNA)



saveRDS(expt.obj, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(expt.obj)
























