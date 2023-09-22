# prepare seurat object
# Vandon Duong

experiment = 'GSE120575'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)


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


read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

express <- read.csv("./data/2019-04-29-GSE120575-dat.csv", row.names=1, header=TRUE)
metadata <- read.csv("./data/GSE120575-patient-data-formatted.csv", row.names=1, header=TRUE)
texpress <- read.tcsv("./data/2019-04-29-GSE120575-dat.csv", row.names=1, header=TRUE)
modmetadata <- read.csv("./data/GSE120575-patient-data-formatted-modified.csv", row.names=1, header=TRUE)
# write.csv(texpress, "./data/2019-04-29-GSE120575-dat-transpose.csv")

CART.combined <- CreateSeuratObject(counts = texpress, meta.data = modmetadata, min.cells = 3, min.features = 200)
CART.combined[["characteristics_response"]] <- modmetadata$characteristics_response
CART.combined[["characteristics_therapy"]] <- modmetadata$characteristics_therapy
CART.combined$characteristics_therapy <- factor(CART.combined$characteristics_therapy, levels = c("anti-CTLA4","anti-PD1","anti-CTLA4+PD1"))

sum(is.na(CART.combined[["characteristics_response"]]))

# assemble metadata
CART.combined@meta.data$Timepoint <- gsub("_.*$", "", modmetadata$characteristics_patinet_ID_Prebaseline_Post_ontreatment)
CART.combined$Timepoint <- factor(CART.combined$Timepoint, levels = c("Pre","Post"))
CART.combined$Timepoint_Response <- paste(CART.combined$Timepoint, CART.combined$characteristics_response, sep = "_")
CART.combined$Timepoint_Response <- factor(CART.combined$Timepoint_Response, levels = c("Pre_Non-responder","Pre_Responder","Post_Non-responder","Post_Responder"))

CART.combined[["percent.mt"]] <- PercentageFeatureSet(CART.combined, pattern = "^MT-")

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = CART.combined, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
CART.combined.CD8pos <- subset(CART.combined,cells=pos_ids)
CART.combined.CD8neg <- subset(CART.combined,cells=neg_ids)

# Quality filter
CART.combined.CD8pos <- subset(CART.combined.CD8pos, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# extract genes
all.genes <- rownames(CART.combined.CD8pos)
write.csv(all.genes, paste('./data/', experiment, '_allgenes.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = CART.combined.CD8pos, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'characteristics_response', ncol=3)
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

# examine other signatures
NK_like = rownames(read.csv("../../signatures/NK-like-dysfunction.csv", row.names=1, header = TRUE))
Wherry_Tex = rownames(read.csv("../../signatures/Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", row.names = 1, header = TRUE))
BBD_Tex = rownames(read.csv("../../signatures/Selli_2023_Blood_TBBDex.csv", row.names = 1, header = TRUE))
PD1_Tex = rownames(read.csv("../../signatures/Cai_2020_Pathology_PD1_Tex.csv", row.names = 1, header = TRUE))

CART.combined.CD8pos <- AddModuleScore(CART.combined.CD8pos, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

CART.combined.CD8pos@meta.data$NKlike_Tex <- scale(CART.combined.CD8pos@meta.data$Signature1)
CART.combined.CD8pos@meta.data$LCMV_Tex <- scale(CART.combined.CD8pos@meta.data$Signature2)
CART.combined.CD8pos@meta.data$BBD_Tex <- scale(CART.combined.CD8pos@meta.data$Signature3)
CART.combined.CD8pos@meta.data$PD1_Tex <- scale(CART.combined.CD8pos@meta.data$Signature4)

CART.combined.CD8pos@meta.data$NKlike_Texi <- integerize(CART.combined.CD8pos@meta.data$NKlike_Tex)
CART.combined.CD8pos@meta.data$LCMV_Texi <- integerize(CART.combined.CD8pos@meta.data$LCMV_Tex)
CART.combined.CD8pos@meta.data$BBD_Texi <- integerize(CART.combined.CD8pos@meta.data$BBD_Tex)
CART.combined.CD8pos@meta.data$PD1_Texi <- integerize(CART.combined.CD8pos@meta.data$PD1_Tex)

CART.combined.CD8pos@meta.data$Signature1 <- NULL
CART.combined.CD8pos@meta.data$Signature2 <- NULL
CART.combined.CD8pos@meta.data$Signature3 <- NULL
CART.combined.CD8pos@meta.data$Signature4 <- NULL


saveRDS(CART.combined.CD8pos, file = paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))

head(CART.combined.CD8pos)




















