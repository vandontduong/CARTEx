# prepare seurat object
# Vandon Duong

experiment = 'GSE136874'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR) # BiocManager::install("beachmat", force = TRUE)
library(beachmat)
# library(celldex)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
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

CART_CD19.data <- Read10X_h5("./data/CD19-28z-filtered_gene_bc_matrices_h5.h5")
CART_GD2.data <- Read10X_h5("./data/GD2-28z-filtered_gene_bc_matrices_h5.h5")
CART_CD19 <- CreateSeuratObject(counts = CART_CD19.data, project = "CD19", min.cells = 3, min.features = 200)
CART_GD2 <- CreateSeuratObject(counts = CART_GD2.data, project = "GD2", min.cells = 3, min.features = 200)
expt.obj <- merge(CART_CD19, y = CART_GD2, add.cell.ids = c("CD19", "GD2"), project = "CART")
rm(list = c('CART_CD19.data', 'CART_GD2.data', 'CART_CD19', 'CART_GD2'))

expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
expt.obj <- NormalizeData(expt.obj)

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

saveRDS(expt.obj, paste('./data/', experiment, '.rds', sep = ''))

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

# regress cell cycle
expt.obj <- ScaleData(expt.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(expt.obj))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(expt.obj), nfeatures.print = 10)

saveRDS(expt.obj, file = paste('./data/', experiment, '_cellcycle.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))

####################################################################################################
######################################## Cell type annotation ######################################
####################################################################################################

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_06_celltype.html
# https://satijalab.org/seurat/articles/covid_sctmapping

ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_cellcycle.rds")
ref.obj <- SetIdent(ref.obj, value = "celltype.l1")
ref.obj <- subset(ref.obj, idents = c("CD4 T", "CD8 T", "other T"))

ref.obj <- UpdateSeuratObject(ref.obj)
ref.obj <- UpdateSCTAssays(ref.obj)

expt.obj <- UpdateSeuratObject(expt.obj)
expt.obj <- UpdateSCTAssays(expt.obj)



anchor <- FindTransferAnchors(reference = ref.obj, query = expt.obj, reference.reduction = "pca", normalization.method = "LogNormalize", dims = 1:50, recompute.residuals = FALSE)

# https://github.com/satijalab/seurat/issues/4117
# https://github.com/satijalab/seurat/issues/5422
expt.obj <- MapQuery(anchorset = anchor, query = expt.obj, reference = ref.obj,  transferdata.args = list(k.weight = FALSE),
                     refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", celltype.l3 = "celltype.l3"),
                     reference.reduction = "pca")

umap_predicted_celltype_l1 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_celltype_l2 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
umap_predicted_celltype_l3 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l3", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

generate_figs(umap_predicted_celltype_l1, paste('./plots/', experiment, '_umap_predicted_celltype_l1', sep = ''))
generate_figs(umap_predicted_celltype_l2, paste('./plots/', experiment, '_umap_predicted_celltype_l2', sep = ''))
generate_figs(umap_predicted_celltype_l3, paste('./plots/', experiment, '_umap_predicted_celltype_l3', sep = ''))

umap_markers <- FeaturePlot(expt.obj, reduction = "umap", features = c("CD4", "CD8A", "CD8B", "PDCD1"), ncol = 2) + NoLegend()
generate_figs(umap_markers, paste('./plots/', experiment, '_umap_markers', sep = ''))


#### SingleR
# https://bioconductor.org/books/release/SingleRBook/introduction.html
# https://github.com/dviraran/SingleR/issues/150
# https://support.bioconductor.org/p/9136971/

test_assay <- LayerData(expt.obj) # LayerData() is the updated function; GetAssayData() was depreciated
ref_assay <- LayerData(ref.obj)

predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = ref_assay, labels = ref.obj$celltype.l2)
table(predictions$labels)

expt.obj[["celltype"]] <- predictions$labels

umap_predicted_celltype_SR <- DimPlot(expt.obj, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_celltype_SR, paste('./plots/', experiment, '_umap_predicted_celltype_SR', sep = ''))




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

expt.obj <- readRDS(paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(expt.obj)

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








