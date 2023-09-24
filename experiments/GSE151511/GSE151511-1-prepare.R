# prepare seurat object
# Vandon Duong

experiment = 'GSE151511'

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

CART.combined <- merge(pt01_obj, y = c(pt02_obj, pt03_obj, pt04_obj, pt05_obj, pt06_obj, pt07_obj, pt08_obj, pt09_obj,
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

CART.combined[["percent.MT"]] <- PercentageFeatureSet(CART.combined, pattern = "^MT-")
head(CART.combined)

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
CART.combined.CD8pos <- subset(CART.combined,cells=pos_ids)
CART.combined.CD8neg <- subset(CART.combined,cells=neg_ids)

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

CART.combined.CD8pos@meta.data$CARTEx_630_repcount <- rowSums(expr > 0)
CART.combined.CD8pos@meta.data$CARTEx_630_reppercent <- CART.combined.CD8pos@meta.data$CARTEx_630_repcount / 630

scatterplot_CARTEx_630_representation <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_630', feature2 = 'CARTEx_630_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_630_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_630_representation', sep = ''))

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(CART.combined.CD8pos))
expr <- t(as.matrix(GetAssayData(CART.combined.CD8pos))[match(common, rownames(as.matrix(GetAssayData(CART.combined.CD8pos)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
CART.combined.CD8pos@meta.data$CARTEx_200 <- Z(scores)
CART.combined.CD8pos@meta.data$CARTEx_200i <- integerize(CART.combined.CD8pos@meta.data$CARTEx_200)

CART.combined.CD8pos@meta.data$CARTEx_200_repcount <- rowSums(expr > 0)
CART.combined.CD8pos@meta.data$CARTEx_200_reppercent <- CART.combined.CD8pos@meta.data$CARTEx_200_repcount / 200

scatterplot_CARTEx_200_representation <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_200', feature2 = 'CARTEx_200_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_200_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_200_representation', sep = ''))

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(CART.combined.CD8pos))
expr <- t(as.matrix(GetAssayData(CART.combined.CD8pos))[match(common, rownames(as.matrix(GetAssayData(CART.combined.CD8pos)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
CART.combined.CD8pos@meta.data$CARTEx_84 <- Z(scores)
CART.combined.CD8pos@meta.data$CARTEx_84i <- integerize(CART.combined.CD8pos@meta.data$CARTEx_84)

CART.combined.CD8pos@meta.data$CARTEx_84_repcount <- rowSums(expr > 0)
CART.combined.CD8pos@meta.data$CARTEx_84_reppercent <- CART.combined.CD8pos@meta.data$CARTEx_84_repcount / 84

scatterplot_CARTEx_84_representation <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'CARTEx_84_reppercent', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_84_representation, paste('./plots/', experiment, '_scatterplot_CARTEx_84_representation', sep = ''))

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




#- Add meta.data (should have added earlier)
CART.combined.CD8pos@meta.data$ClinicalResponse <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$orig.ident,
                                                                   from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                                   to = c('CR', 'PD', 'PD', 'PD', 'CR', 'NE', 'CR', 'CR', 'CR', 'CR', 'PD', 'CR', 'PD', 'CR', 'PD', 'CR', 'PD', 'PD', 'PD', 'PR', 'PD', 'PD', 'PD', 'PD'))
CART.combined.CD8pos@meta.data$ClinicalResponse <- factor(CART.combined.CD8pos@meta.data$ClinicalResponse, levels = c('CR', 'PR', 'PD', 'NE'))
CART.combined.CD8pos@meta.data$EMR <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$orig.ident,
                                                      from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                      to = c('GOOD', 'NE', 'BAD', 'NE', 'GOOD', 'BAD', 'GOOD', 'GOOD', 'NE', 'GOOD', 'BAD', 'NE', 'BAD', 'NE', 'NE', 'GOOD', 'GOOD', 'BAD', 'BAD', 'BAD', 'GOOD', 'BAD', 'BAD', 'NE'))
CART.combined.CD8pos@meta.data$EMR <- factor(CART.combined.CD8pos@meta.data$EMR, levels = c('GOOD', 'BAD', 'NE'))

CART.combined.CD8pos@meta.data$Responder <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$orig.ident,
                                                            from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                            to = c('Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Responder', 'NE', 'Responder', 'Responder', 'Responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Responder', 'Non-responder', 'Non-responder', 'Non-responder', 'Non-responder'))
CART.combined.CD8pos@meta.data$Responder <- factor(CART.combined.CD8pos@meta.data$Responder, levels = c('Responder', 'Non-responder', 'NE'))

CART.combined.CD8pos@meta.data$CRS <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$orig.ident,
                                                      from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                      to = c(2, 2, 4, 3, 1, 2, 1, 2, 2, 2, 2, 3, 1, 2, 2, 2, 1, 1, 1, 1, 3, 1, 2, 2))
CART.combined.CD8pos@meta.data$CRS <- factor(CART.combined.CD8pos@meta.data$CRS, levels = c(1, 2, 3, 4))

CART.combined.CD8pos@meta.data$ICANS <- plyr::mapvalues(x = CART.combined.CD8pos@meta.data$orig.ident,
                                                      from = c('pt01', 'pt02', 'pt03', 'pt04', 'pt05', 'pt06', 'pt07', 'pt08', 'pt09', 'pt10', 'pt11', 'pt12', 'pt13', 'pt14', 'pt15', 'pt16', 'pt17', 'pt18', 'pt19', 'pt20', 'pt21', 'pt22', 'pt23', 'pt24'),
                                                      to = c(3, 0, 4, 3, 0, 4, 3, 0, 3, 0, 3, 3, 0, 2, 3, 1, 0, 2, 2, 0, 3, 3, 3, 0))
CART.combined.CD8pos@meta.data$ICANS <- factor(CART.combined.CD8pos@meta.data$ICANS, levels = c(0, 1, 2, 3, 4))




# Novelty score
CART.combined.CD8pos$log10GenesPerUMI <- log10(CART.combined.CD8pos$nFeature_RNA) / log10(CART.combined.CD8pos$nCount_RNA)



saveRDS(CART.combined.CD8pos, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(CART.combined.CD8pos)
























