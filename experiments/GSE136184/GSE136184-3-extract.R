# extract newborn data
# Vandon Duong

set.seed(123)
experiment = 'GSE136184'
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_annotated_cross.rds', sep = ''))
# expt.obj <- SetIdent(expt.obj, value = "AgeGroup")
head(expt.obj)

expt.obj@meta.data$extract.ident <- paste(expt.obj@meta.data$AgeGroup, expt.obj@meta.data$monaco, sep = ", ")


# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_630 <- Z(scores)
expt.obj@meta.data$CARTEx_630i <- integerize(expt.obj@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_200 <- Z(scores)
expt.obj@meta.data$CARTEx_200i <- integerize(expt.obj@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_84 <- Z(scores)
expt.obj@meta.data$CARTEx_84i <- integerize(expt.obj@meta.data$CARTEx_84)


activation_sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), row.names = 1, header = TRUE))

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

ridgeplt_CARTEx_84 <- RidgePlot(expt.obj, features = c("CARTEx_84"), group.by = "AgeGroup2") + xlim(c(NA, 4))
generate_figs(ridgeplt_CARTEx_84, paste('./plots/', experiment, '_extract_ridgeplt_CARTEx_84', sep = ''), c(8,4))


cellIDs_newborn <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Newborn'])
cellIDs_under30 <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Under 30'])
cellIDs_elderly <- Cells(expt.obj[, expt.obj[['AgeGroup2']] == 'Elderly'])

cellIDs_naiveCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Naive CD8 T cells'])
cellIDs_terminalCD8T <- Cells(expt.obj[, expt.obj[['monaco']] == 'Terminal effector CD8 T cells'])

expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% intersect(union(cellIDs_newborn, cellIDs_under30), cellIDs_naiveCD8T))]
expt.obj.elderly_terminal <- expt.obj[,(colnames(expt.obj) %in% intersect(cellIDs_elderly, cellIDs_terminalCD8T))]

# expt.obj.young_naive@meta.data$extract.ident <- "Young, naive CD8 T cells"
# expt.obj.elderly_terminal@meta.data$extract.ident <- "Elderly, terminal effector CD8 T cells"


SelectSortCells <- function(atlas, module, direction, counts){
  md <- atlas@meta.data
  if (direction == 'high'){cells <- rownames(md[order(-md[module]), ][1:counts,])}
  if (direction == 'low'){cells <- rownames(md[order(md[module]), ][1:counts,])}
  return(cells)
}

cellIDs_YN_low <- SelectSortCells(expt.obj.young_naive, 'CARTEx_84', 'low', 1500)
cellIDs_ET_high <- SelectSortCells(expt.obj.elderly_terminal, 'CARTEx_84', 'high', 1500)

# expt.obj.young_naive <- expt.obj[,(colnames(expt.obj) %in% cellIDs_YN_low)]
# expt.obj.elderly_terminal <- expt.obj[,(colnames(expt.obj) %in% cellIDs_ET_high)]


# rm(expt.obj)

# expt.obj <- merge(expt.obj.young_naive, y = expt.obj.elderly_terminal, project = "aging_extract")

expt.obj <- expt.obj[,(colnames(expt.obj) %in% union(cellIDs_YN_low, cellIDs_ET_high))]

# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# extract genes
all.genes <- rownames(expt.obj)
write.csv(all.genes, paste('./data/', experiment, '_allgenes_aging_extract.csv', sep = ''))

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.mt"), group.by = 'orig.ident', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_prepare_vlnplot_quality_aging_extract', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_prepare_varplt_labeled_aging_extract', sep = ''))


expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

Idents(expt.obj) <- "extract.ident"
expt.obj <- RenameIdents(expt.obj, `Old, Terminal effector CD8 T cells` = "Old Terminal", `Young, Naive CD8 T cells` = "Young Naive")
table(Idents(expt.obj))
expt.obj[['extract.ident']] <- Idents(expt.obj)
table(expt.obj[['extract.ident']])

saveRDS(expt.obj, file = paste('./data/', experiment, '_aging_extract.rds', sep = ''))

# expt.obj <- readRDS(paste('./data/', experiment, '_aging_extract.rds', sep = ''))
















# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical <- VlnPlot(expt.obj, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "AgeGroup2", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical, paste('./plots/', experiment, '_vlnplot_naiveT_canonical', sep = ''))

####################################################################################################
######################################## Extract newborn data ######################################
####################################################################################################


expt.obj <- SetIdent(expt.obj, value = "AgeGroup2")
expt.obj <- subset(expt.obj, idents = c("Newborn"))


# Quality filter
expt.obj <- subset(expt.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Examine features
vlnplot_quality <- VlnPlot(object = expt.obj, features = c('nFeature_RNA','nCount_RNA', "percent.MT"), group.by = 'AgeGroup2', ncol=3)
generate_figs(vlnplot_quality, paste('./plots/', experiment, '_vlnplot_quality_newborn', sep = ''))

# Variable features and initial UMAP analysis
expt.obj <- FindVariableFeatures(expt.obj, selection.method = "vst", nfeatures = 2000)
top10_CART <- head(VariableFeatures(expt.obj), 10)
varplt <- VariableFeaturePlot(expt.obj)
varplt_labeled <- LabelPoints(plot = varplt, points = top10_CART, repel = TRUE)
generate_figs(varplt_labeled, paste('./plots/', experiment, '_varplt_labeled_newborn', sep = ''))

expt.obj <- ScaleData(expt.obj, features = all.genes)
expt.obj <- RunPCA(expt.obj, features = VariableFeatures(object = expt.obj), npcs = 40)

expt.obj <- FindNeighbors(expt.obj, dims = 1:10)
expt.obj <- FindClusters(expt.obj, resolution = 0.5)
expt.obj <- RunUMAP(expt.obj, dims = 1:10)

saveRDS(expt.obj, file = paste('./data/', experiment, '_newborn.rds', sep = ''))


umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters_newborn', sep = ''))

# Naive T cells have high expression of CD62L and CCR7, but low or no expression of IL2RA, CD44, CD69, and PTPRC
vlnplot_naiveT_canonical_newborn <- VlnPlot(expt.obj, c("CD62L", "CCR7", "IL2RA", "CD44", "CD69", "PTPRC"), group.by = "seurat_clusters", pt.size = 0)
generate_figs(vlnplot_naiveT_canonical_newborn, paste('./plots/', experiment, '_vlnplot_naiveT_canonical_newborn', sep = ''))







