# prepare seurat object
# Vandon Duong

experiment = 'GSE164378'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)

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

pbmcref[["percent.MT"]] <- PercentageFeatureSet(pbmcref, pattern = "^MT-")

# Extract only the un-vaccinated (time 0) cells and doublets
pbmcref$donor_time <- paste(pbmcref$donor,"_",pbmcref$time,sep = "")
Idents(pbmcref) <- "donor_time"
pbmcref <- subset(x = pbmcref, idents = c("P1_0", "P2_0","P3_0","P4_0","P5_0","P6_0","P7_0","P8_0"))

Idents(pbmcref) <- "celltype.l2"
pbmcref <- subset(pbmcref, idents = "Doublet", invert=T)

# Quality filter
pbmcref <- subset(pbmcref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.MT < 10)

# extract genes
all.genes <- rownames(pbmcref)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = pbmcref, features = c('nFeature_RNA','nCount_RNA', "percent.MT"), group.by = 'orig.ident', ncol=3)
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
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(pbmcref))
expr <- t(as.matrix(GetAssayData(pbmcref))[match(common, rownames(as.matrix(GetAssayData(pbmcref)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref@meta.data$CARTEx_630 <- Z(scores)
pbmcref@meta.data$CARTEx_630i <- integerize(pbmcref@meta.data$CARTEx_630)

pbmcref@meta.data$CARTEx_630_repcount <- rowSums(expr > 0)
pbmcref@meta.data$CARTEx_630_reppercent <- pbmcref@meta.data$CARTEx_630_repcount / 630

scatterplot_CARTEx_630_representation <- FeatureScatter(pbmcref, feature1 = 'CARTEx_630', feature2 = 'CARTEx_630_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_630_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_630_representation', sep = ''))

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(pbmcref))
expr <- t(as.matrix(GetAssayData(pbmcref))[match(common, rownames(as.matrix(GetAssayData(pbmcref)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref@meta.data$CARTEx_200 <- Z(scores)
pbmcref@meta.data$CARTEx_200i <- integerize(pbmcref@meta.data$CARTEx_200)

pbmcref@meta.data$CARTEx_200_repcount <- rowSums(expr > 0)
pbmcref@meta.data$CARTEx_200_reppercent <- pbmcref@meta.data$CARTEx_200_repcount / 200

scatterplot_CARTEx_200_representation <- FeatureScatter(pbmcref, feature1 = 'CARTEx_200', feature2 = 'CARTEx_200_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_200_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_200_representation', sep = ''))

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(pbmcref))
expr <- t(as.matrix(GetAssayData(pbmcref))[match(common, rownames(as.matrix(GetAssayData(pbmcref)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
pbmcref@meta.data$CARTEx_84 <- Z(scores)
pbmcref@meta.data$CARTEx_84i <- integerize(pbmcref@meta.data$CARTEx_84)

pbmcref@meta.data$CARTEx_84_repcount <- rowSums(expr > 0)
pbmcref@meta.data$CARTEx_84_reppercent <- pbmcref@meta.data$CARTEx_84_repcount / 84

scatterplot_CARTEx_84_representation <- FeatureScatter(pbmcref, feature1 = 'CARTEx_84', feature2 = 'CARTEx_84_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_84_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_84_representation', sep = ''))

####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

pbmcref <- AddModuleScore(pbmcref, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
pbmcref@meta.data$Activation <- scale(pbmcref@meta.data$State1)
pbmcref@meta.data$Anergy <- scale(pbmcref@meta.data$State2)
pbmcref@meta.data$Stemness <- scale(pbmcref@meta.data$State3)
pbmcref@meta.data$Senescence <- scale(pbmcref@meta.data$State4)

pbmcref@meta.data$Activationi <- integerize(pbmcref@meta.data$Activation)
pbmcref@meta.data$Anergyi <- integerize(pbmcref@meta.data$Anergy)
pbmcref@meta.data$Stemnessi <- integerize(pbmcref@meta.data$Stemness)
pbmcref@meta.data$Senescencei <- integerize(pbmcref@meta.data$Senescence)

pbmcref@meta.data$State1 <- NULL
pbmcref@meta.data$State2 <- NULL
pbmcref@meta.data$State3 <- NULL
pbmcref@meta.data$State4 <- NULL

# examine other signatures
NK_like = rownames(read.csv("../../signatures/NK-like-dysfunction.csv", row.names=1, header = TRUE))
Wherry_Tex = rownames(read.csv("../../signatures/Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", row.names = 1, header = TRUE))
BBD_Tex = rownames(read.csv("../../signatures/Selli_2023_Blood_TBBDex.csv", row.names = 1, header = TRUE))
PD1_Tex = rownames(read.csv("../../signatures/Cai_2020_Pathology_PD1_Tex.csv", row.names = 1, header = TRUE))

pbmcref <- AddModuleScore(pbmcref, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

pbmcref@meta.data$NKlike_Tex <- scale(pbmcref@meta.data$Signature1)
pbmcref@meta.data$LCMV_Tex <- scale(pbmcref@meta.data$Signature2)
pbmcref@meta.data$BBD_Tex <- scale(pbmcref@meta.data$Signature3)
pbmcref@meta.data$PD1_Tex <- scale(pbmcref@meta.data$Signature4)

pbmcref@meta.data$NKlike_Texi <- integerize(pbmcref@meta.data$NKlike_Tex)
pbmcref@meta.data$LCMV_Texi <- integerize(pbmcref@meta.data$LCMV_Tex)
pbmcref@meta.data$BBD_Texi <- integerize(pbmcref@meta.data$BBD_Tex)
pbmcref@meta.data$PD1_Texi <- integerize(pbmcref@meta.data$PD1_Tex)

pbmcref@meta.data$Signature1 <- NULL
pbmcref@meta.data$Signature2 <- NULL
pbmcref@meta.data$Signature3 <- NULL
pbmcref@meta.data$Signature4 <- NULL

Idents(pbmcref) <- "celltype.l2"
saveRDS(pbmcref, file = paste('./data/', experiment, '_scored.rds', sep = ''))

head(pbmcref)





