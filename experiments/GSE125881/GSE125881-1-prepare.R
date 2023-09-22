# prepare seurat object
# Vandon Duong

experiment = 'GSE125881'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(stringr)

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


counts_data <- read.csv(file = "./data/GSE125881_raw.expMatrix.csv", header = TRUE, row.names = 1)
CART.combined <- CreateSeuratObject(counts = counts_data, project = "Kinetics", min.cells = 3, min.features = 200)

CART.combined[["percent.MT"]] <- PercentageFeatureSet(CART.combined, pattern = "^MT-")

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
CART.combined.CD8pos <- subset(CART.combined,cells=pos_ids)
CART.combined.CD8neg <- subset(CART.combined,cells=neg_ids)

#- Add meta.data
CART.combined.CD8pos@meta.data$SampleID <- sapply(rownames(CART.combined.CD8pos@meta.data), function(x) str_split(x, '[.]')[[1]][2])
CART.combined.CD8pos@meta.data$Patient <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$SampleID,
                                                          from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                          to = c('CLL-1', 'CLL-1', 'CLL-1', 'CLL-1', 'NHL-9', 'NHL-9', 'NHL-9', 'CLL-2', 'CLL-2', 'CLL-2', 'CLL-2', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-9'))
CART.combined.CD8pos@meta.data$Patient <- factor(CART.combined.CD8pos@meta.data$Patient, levels = c('CLL-1', 'CLL-2', 'NHL-9', 'NHL-10'))
CART.combined.CD8pos@meta.data$TimePoint <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$SampleID,
                                                            from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                            to = c('IP', 'd21', 'd38', 'd112', 'd12', 'd29', 'd102', 'IP', 'd12', 'd29', 'd83', 'IP', 'd12', 'd28', 'd89', 'IP'))
CART.combined.CD8pos@meta.data$TimePoint <- factor(CART.combined.CD8pos@meta.data$TimePoint, levels = c('IP', 'd12', 'd21', 'd28', 'd29', 'd38', 'd83', 'd89', 'd102', 'd112'))
CART.combined.CD8pos@meta.data$Group <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$SampleID,
                                                        from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                        to = c('IP', 'ExpansionPeak', 'Contraction', 'Late', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP'))
CART.combined.CD8pos@meta.data$Group <- factor(CART.combined.CD8pos@meta.data$Group, levels = c('IP', 'ExpansionPeak', 'Contraction', 'Late'))
CART.combined.CD8pos@meta.data$Disease <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$SampleID,
                                                          from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                          to = c('CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'NHL', 'NHL'))

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
rm(CART.combined)
rm(CART.combined.CD8neg)
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

rm(CART.combined)
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














