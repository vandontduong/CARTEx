# prepare seurat object
# Vandon Duong

experiment = 'GSE136184'

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

CART.combined.CD8pos <- readRDS(paste('./data/GSE136184_seu_obj2.Rds', sep = ''))
CART.combined.CD8pos[["percent.MT"]] <- PercentageFeatureSet(CART.combined.CD8pos, pattern = "^MT-")

#- Add meta.data
CART.combined.CD8pos@meta.data$Sex <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$Code,
                                                          from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                          to = c('M', 'F', 'M', 'F', 'F', 'M', 'F', 'M', 'F', 'F', 'M', 'M', 'F', 'M', 'M', 'F', 'F', 'M', 'M', 'F', 'F', 'M', 'F', 'M'))
CART.combined.CD8pos@meta.data$Sex <- factor(CART.combined.CD8pos@meta.data$Sex, levels = c('M', 'F'))

CART.combined.CD8pos@meta.data <- mutate(CART.combined.CD8pos@meta.data, Age = case_when(
  Code == 'sc_d1' ~ 0,
  Code == 'sc_d2' ~ 0,
  Code == 'sc_d3' ~ 26,
  Code == 'sc_d4' ~ 28,
  Code == 'sc_d5' ~ 37,
  Code == 'sc_d6' ~ 46,
  Code == 'sc_d7' ~ 53,
  Code == 'sc_d8' ~ 64,
  Code == 'sc_d9' ~ 68,
  Code == 'sc_d10' ~ 76,
  Code == 'sc_d11' ~ 77,
  Code == 'sc_d12' ~ 83,
  Code == 'sc_d13' ~ 83,
  Code == 'sc_d14' ~ 84,
  Code == 'sc_d15' ~ 88,
  Code == 'sc_d16' ~ 90,
  Code == 'L1' & visit == 1 ~ 30,
  Code == 'L1' & visit == 2 ~ 39,
  Code == 'L2' & visit == 1 ~ 32,
  Code == 'L2' & visit == 2 ~ 40,
  Code == 'L3' & visit == 1 ~ 48,
  Code == 'L3' & visit == 2 ~ 56,
  Code == 'L4' & visit == 1 ~ 54,
  Code == 'L4' & visit == 2 ~ 65,
  Code == 'L5' & visit == 1 ~ 60,
  Code == 'L5' & visit == 2 ~ 65,
  Code == 'L5' & visit == 3 ~ 70,
  Code == 'L6' & visit == 1 ~ 60,
  Code == 'L6' & visit == 2 ~ 71,
  Code == 'L7' & visit == 1 ~ 69,
  Code == 'L7' & visit == 2 ~ 78,
  Code == 'L8' & visit == 1 ~ 69,
  Code == 'L8' & visit == 2 ~ 77
  ))

CART.combined.CD8pos@meta.data <- mutate(CART.combined.CD8pos@meta.data, AgeGroup = case_when(
  Age <= 30 ~ 'Young',
  Age <= 69 ~ 'Middle',
  Age > 69 ~ 'Old'
  ))
  
CART.combined.CD8pos@meta.data$AgeGroup <- factor(CART.combined.CD8pos@meta.data$AgeGroup, levels = c('Young', 'Middle', 'Old'))

CART.combined.CD8pos@meta.data <- mutate(CART.combined.CD8pos@meta.data, AgeGroup2 = case_when(
  Age == 0 ~ 'Newborn',
  Age < 30 ~ 'Under 30',
  Age < 50 ~ 'Under 50',
  Age < 70 ~ 'Under 70',
  Age >= 70 ~ 'Elderly'
  ))

CART.combined.CD8pos@meta.data$AgeGroup2 <- factor(CART.combined.CD8pos@meta.data$AgeGroup2, levels = c('Newborn', 'Under 30', 'Under 50', 'Under 70', 'Elderly'))

CART.combined.CD8pos@meta.data$Cohort <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$Code,
                                                      from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                      to = c('Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal'))
CART.combined.CD8pos@meta.data$Cohort <- factor(CART.combined.CD8pos@meta.data$Cohort, levels = c('Cross-sectional', 'Longitudinal'))

# Quality filter
CART.combined.CD8pos <- subset(CART.combined.CD8pos, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.MT < 10)

# extract genes
all.genes <- rownames(CART.combined.CD8pos)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = CART.combined.CD8pos, features = c('nFeature_RNA','nCount_RNA', "percent.MT"), group.by = 'orig.ident', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality', sep = ''))

# Variable features and initial UMAP analysis
CART.combined.CD8pos <- FindVariableFeatures(CART.combined.CD8pos, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(CART.combined.CD8pos), 10)
varplt <- VariableFeaturePlot(CART.combined.CD8pos)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled', sep = ''))

CART.combined.CD8pos <- ScaleData(CART.combined.CD8pos, features = all.genes)
CART.combined.CD8pos <- RunPCA(CART.combined.CD8pos, features = VariableFeatures(object = CART.combined.CD8pos), npcs = 40)

CART.combined.CD8pos <- FindNeighbors(CART.combined.CD8pos, dims = 1:10)
CART.combined.CD8pos <- FindClusters(CART.combined.CD8pos, resolution = 0.5)
CART.combined.CD8pos <- RunUMAP(CART.combined.CD8pos, dims = 1:10)

saveRDS(CART.combined.CD8pos, file = paste('./data/', experiment, '.rds', sep = ''))

CART.combined.CD8pos <- readRDS(paste('./data/', experiment, '.rds', sep = ''))


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

CART.combined.CD8pos <- CellCycleScoring(CART.combined.CD8pos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(CART.combined.CD8pos[[]])

# Visualize the distribution of cell cycle markers across
ridgeplt <- RidgePlot(CART.combined.CD8pos, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
generate_figs(ridgeplt, paste('./plots/', experiment, '_ridgeplt', sep = ''))

CART.combined.CD8pos <- RunPCA(CART.combined.CD8pos, features = c(s.genes, g2m.genes))
# DimPlot(CART.combined.CD8pos)

# regress cell cycle
CART.combined.CD8pos <- ScaleData(CART.combined.CD8pos, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CART.combined.CD8pos))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
CART.combined.CD8pos <- RunPCA(CART.combined.CD8pos, features = VariableFeatures(CART.combined.CD8pos), nfeatures.print = 10)

saveRDS(CART.combined.CD8pos, file = paste('./data/', experiment, '_CD8pos_cellcycle.rds', sep = ''))

CART.combined.CD8pos <- readRDS(paste('./data/', experiment, '_CD8pos_cellcycle.rds', sep = ''))


####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(CART.combined.CD8pos))
expr <- t(as.matrix(GetAssayData(CART.combined.CD8pos))[match(common, rownames(as.matrix(GetAssayData(CART.combined.CD8pos)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
CART.combined.CD8pos@meta.data$CARTEx_630 <- Z(scores)
CART.combined.CD8pos@meta.data$CARTEx_630i <- integerize(CART.combined.CD8pos@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(CART.combined.CD8pos))
expr <- t(as.matrix(GetAssayData(CART.combined.CD8pos))[match(common, rownames(as.matrix(GetAssayData(CART.combined.CD8pos)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
CART.combined.CD8pos@meta.data$CARTEx_200 <- Z(scores)
CART.combined.CD8pos@meta.data$CARTEx_200i <- integerize(CART.combined.CD8pos@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(CART.combined.CD8pos))
expr <- t(as.matrix(GetAssayData(CART.combined.CD8pos))[match(common, rownames(as.matrix(GetAssayData(CART.combined.CD8pos)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
CART.combined.CD8pos@meta.data$CARTEx_84 <- Z(scores)
CART.combined.CD8pos@meta.data$CARTEx_84i <- integerize(CART.combined.CD8pos@meta.data$CARTEx_84)


####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

CART.combined.CD8pos <- AddModuleScore(CART.combined.CD8pos, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
CART.combined.CD8pos@meta.data$Activation <- scale(CART.combined.CD8pos@meta.data$State1)
CART.combined.CD8pos@meta.data$Anergy <- scale(CART.combined.CD8pos@meta.data$State2)
CART.combined.CD8pos@meta.data$Stemness <- scale(CART.combined.CD8pos@meta.data$State3)
CART.combined.CD8pos@meta.data$Senescence <- scale(CART.combined.CD8pos@meta.data$State4)

CART.combined.CD8pos@meta.data$Activationi <- integerize(CART.combined.CD8pos@meta.data$Activation)
CART.combined.CD8pos@meta.data$Anergyi <- integerize(CART.combined.CD8pos@meta.data$Anergy)
CART.combined.CD8pos@meta.data$Stemnessi <- integerize(CART.combined.CD8pos@meta.data$Stemness)
CART.combined.CD8pos@meta.data$Senescencei <- integerize(CART.combined.CD8pos@meta.data$Senescence)

CART.combined.CD8pos@meta.data$State1 <- NULL
CART.combined.CD8pos@meta.data$State2 <- NULL
CART.combined.CD8pos@meta.data$State3 <- NULL
CART.combined.CD8pos@meta.data$State4 <- NULL

saveRDS(CART.combined.CD8pos, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(CART.combined.CD8pos)







