# integrate seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE136874'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)

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

# https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_rpca.html

expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))
aging.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE136184/data/GSE136184_annotated_cross.rds")
ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_allT_annotated.rds")
# ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_cellcycle.rds") # T cell reference ?
# ref.obj <- SetIdent(ref.obj, value = "celltype.l1")
# ref.obj <- subset(ref.obj, idents = c("CD4 T", "CD8 T", "other T"))

# add experiment identifier
expt.obj@meta.data[['identifier']] <- experiment
aging.obj@meta.data[['identifier']] <- "GSE136184"
ref.obj@meta.data[['identifier']] <- "GSE164378"

expt.obj@meta.data[['identifier2']] <- expt.obj@meta.data$orig.ident
aging.obj@meta.data[['identifier2']] <- aging.obj@meta.data$AgeGroup2
ref.obj@meta.data[['identifier2']] <- ref.obj@meta.data$monaco

expt.obj@meta.data[["split.ident"]] <- "Query"
aging.obj@meta.data[["split.ident"]] <- "Query"
ref.obj@meta.data[["split.ident"]] <- "Reference"

integration.list <- c(experiment = expt.obj, aging = aging.obj, reference = ref.obj)

# normalize and identify variable features for each dataset independently
integration.list <- lapply(X = integration.list, FUN = function(x) {
  x <- UpdateSeuratObject(x)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

integration.features <- SelectIntegrationFeatures(object.list = integration.list)

integration.list <- lapply(X = integration.list, FUN = function(x) {
  x <- ScaleData(x, features = integration.features)
  x <- RunPCA(x, features = integration.features)
})

# integration.list <- PrepSCTIntegration(object.list = integration.list, anchor.features = integration.features)
integration.anchors <- FindIntegrationAnchors(object.list = integration.list, anchor.features = integration.features, reduction = "rpca", reference = c(3)) # k.anchor = 100, k.filter = NA, dims = 1:10
integration.obj <- IntegrateData(anchorset = integration.anchors)
DefaultAssay(integration.obj) <- "integrated"

# Run the standard workflow for visualization and clustering
integration.obj <- ScaleData(integration.obj)
integration.obj <- RunPCA(integration.obj, npcs = 30)
integration.obj <- RunUMAP(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindNeighbors(integration.obj, reduction = "pca", dims = 1:30)
integration.obj <- FindClusters(integration.obj, resolution = 0.5)

head(integration.obj)

saveRDS(integration.obj, file = paste('./data/', experiment, '_integrated.rds', sep = ''))

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

integrated_umap_identifier <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier2", split.by = "identifier", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier, paste('./plots/', experiment, '_integrated_umap_identifier', sep = ''), c(12, 5))

integrated_umap_identifier_expts <- DimPlot(integration.obj, reduction = "umap", group.by = "identifier", split.by = "split.ident", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_identifier_expts, paste('./plots/', experiment, '_integrated_umap_identifier_expts', sep = ''), c(12, 5))

integrated_umap_seurat_clusters <- DimPlot(integration.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_seurat_clusters, paste('./plots/', experiment, '_integrated_umap_seurat_clusters', sep = ''))

integrated_umap_azimuth <- DimPlot(integration.obj, reduction = "umap", group.by = "azimuth", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_azimuth, paste('./plots/', experiment, '_integrated_umap_azimuth', sep = ''))

integrated_umap_monaco <- DimPlot(integration.obj, reduction = "umap", group.by = "monaco", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_monaco, paste('./plots/', experiment, '_integrated_umap_monaco', sep = ''))

integrated_umap_dice <- DimPlot(integration.obj, reduction = "umap", group.by = "dice", shuffle = TRUE, seed = 123, raster = FALSE)
generate_figs(integrated_umap_dice, paste('./plots/', experiment, '_integrated_umap_dice', sep = ''))

# BAR CHARTS of CELL TYPES
md <- integration.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_integrated_barplot_dice_seurat_clusters', sep = ''))

####################################################################################################
####################################### Examine query datasets #####################################
####################################################################################################

integration.obj <- readRDS(paste('./data/', experiment, '_integrated.rds', sep = ''))

query.list <- SplitObject(integration.obj, split.by = "split.ident")
query.obj <- query.list$Query
rm(aging.obj, expt.obj, ref.obj, integration.obj)

md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = query.obj, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = query.obj, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = query.obj, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
query.obj <- subset(query.obj, cells=pos_ids)

# CD8 T cell detected in any reference
query.obj$CD8Tref_1 <- with(query.obj, ifelse((query.obj@meta.data$azimuth == "CD8 Naive") | (query.obj@meta.data$azimuth == "CD8 TEM") | (query.obj@meta.data$azimuth == "CD8 Proliferating") | 
                                              (query.obj@meta.data$monaco == "Naive CD8 T cells") | (query.obj@meta.data$monaco == "Effector CD8 T cells") | (query.obj@meta.data$monaco == "Central Memory CD8 T cells") | (query.obj@meta.data$monaco == "Terminal effector CD8 T cells") |
                                              (query.obj@meta.data$dice == "T cells, CD8+, naive") | (query.obj@meta.data$dice == "T cells, CD8+, naive, stimulated") , 
                                            'CD8 T cell', 'Other'))

dplyr::count(query.obj@meta.data, CD8Tref_1, sort = TRUE)

umap_predicted_CD8Tref_1 <- DimPlot(query.obj, reduction = "umap", group.by = "CD8Tref_1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_CD8Tref_1, paste('./plots/', experiment, '_query_umap_predicted_CD8Tref_1', sep = ''))

query.obj <- SetIdent(query.obj, value = "CD8Tref_1")
query.obj <- subset(query.obj, idents = c("CD8 T cell"))

head(query.obj)


####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_630 <- Z(scores)
query.obj@meta.data$CARTEx_630i <- integerize(query.obj@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_200 <- Z(scores)
query.obj@meta.data$CARTEx_200i <- integerize(query.obj@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_84 <- Z(scores)
query.obj@meta.data$CARTEx_84i <- integerize(query.obj@meta.data$CARTEx_84)


####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation_sig <- rownames(read.csv("../../signatures/panther-activation.csv", header = TRUE, row.names = 1))
anergy_sig <- rownames(read.csv("../../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", header = TRUE, row.names = 1))
stemness_sig <- rownames(read.csv("../../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
senescence_sig <- rownames(read.csv("../../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

query.obj <- AddModuleScore(query.obj, features = list(activation_sig, anergy_sig, stemness_sig, senescence_sig), name="State", search = TRUE)

# z score normalization
query.obj@meta.data$Activation <- scale(query.obj@meta.data$State1)
query.obj@meta.data$Anergy <- scale(query.obj@meta.data$State2)
query.obj@meta.data$Stemness <- scale(query.obj@meta.data$State3)
query.obj@meta.data$Senescence <- scale(query.obj@meta.data$State4)

query.obj@meta.data$Activationi <- integerize(query.obj@meta.data$Activation)
query.obj@meta.data$Anergyi <- integerize(query.obj@meta.data$Anergy)
query.obj@meta.data$Stemnessi <- integerize(query.obj@meta.data$Stemness)
query.obj@meta.data$Senescencei <- integerize(query.obj@meta.data$Senescence)

query.obj@meta.data$State1 <- NULL
query.obj@meta.data$State2 <- NULL
query.obj@meta.data$State3 <- NULL
query.obj@meta.data$State4 <- NULL

saveRDS(query.obj, file = paste('./data/', experiment, '_query_scored.rds', sep = ''))

query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

head(query.obj)

# CARTEx violin plot
query.obj <- SetIdent(query.obj, value = "identifier2")
split.ident.order = c("CD19", "GD2", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly")
Idents(query.obj) <- factor(Idents(query.obj), levels = split.ident.order)

query.obj$identifier2 <- factor(query.obj$identifier2, levels = c("CD19", "GD2", "Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
vlnplot_CARTEx_84 <- VlnPlot(query.obj, features = c("CARTEx_84"), group.by = 'identifier2', y.max = 6, pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 5)
generate_figs(vlnplot_CARTEx_84, paste('./plots/', experiment, '_query_vlnplot_CARTEx_84', sep = ''))

# BAR CHARTS of CELL TYPES
md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("azimuth", "monaco", "dice")]

md_temp <- md[, .N, by = c('azimuth', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_identifier2, paste('./plots/', experiment, '_query_barplot_azimuth_identifier2', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_identifier2, paste('./plots/', experiment, '_query_barplot_monaco_identifier2', sep = ''))

md_temp <- md[, .N, by = c('dice', 'identifier2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_identifier2 <- ggplot(md_temp, aes(x = identifier2, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_identifier2, paste('./plots/', experiment, '_query_barplot_dice_identifier2', sep = ''))






