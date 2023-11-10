# prepare seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE136184'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR)
library(scuttle)
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

expt.obj <- readRDS(paste('./data/GSE136184_seu_obj2.Rds', sep = ''))
expt.obj[["percent.mt"]] <- PercentageFeatureSet(expt.obj, pattern = "^MT-")

# Filter for CD8A +  cells, removing CD4+
CD4_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD4",]
CD8A_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8A",]
CD8B_expression <- GetAssayData(object = expt.obj, assay = "RNA", slot = "data")["CD8B",]
pos_ids <- names(which(CD8A_expression > 0 & CD8B_expression > 0 & CD4_expression == 0))
neg_ids <- names(which(CD8A_expression == 0 & CD8B_expression == 0 & CD4_expression > 0))
expt.obj <- subset(expt.obj,cells=pos_ids)

#- Add meta.data
expt.obj@meta.data$Sex <- plyr::mapvalues(x = expt.obj@meta.data$Code,
                                                          from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                          to = c('M', 'F', 'M', 'F', 'F', 'M', 'F', 'M', 'F', 'F', 'M', 'M', 'F', 'M', 'M', 'F', 'F', 'M', 'M', 'F', 'F', 'M', 'F', 'M'))
expt.obj@meta.data$Sex <- factor(expt.obj@meta.data$Sex, levels = c('M', 'F'))

expt.obj@meta.data <- mutate(expt.obj@meta.data, Age = case_when(
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

expt.obj@meta.data <- mutate(expt.obj@meta.data, AgeGroup = case_when(
  Age <= 30 ~ 'Young',
  Age <= 69 ~ 'Middle',
  Age > 69 ~ 'Old'
  ))
  
expt.obj@meta.data$AgeGroup <- factor(expt.obj@meta.data$AgeGroup, levels = c('Young', 'Middle', 'Old'))

expt.obj@meta.data <- mutate(expt.obj@meta.data, AgeGroup2 = case_when(
  Age == 0 ~ 'Newborn',
  Age < 30 ~ 'Under 30',
  Age < 50 ~ 'Under 50',
  Age < 70 ~ 'Under 70',
  Age >= 70 ~ 'Elderly'
  ))

expt.obj@meta.data$AgeGroup2 <- factor(expt.obj@meta.data$AgeGroup2, levels = c('Newborn', 'Under 30', 'Under 50', 'Under 70', 'Elderly'))

expt.obj@meta.data$Cohort <- plyr::mapvalues(x = expt.obj@meta.data$Code,
                                                      from = c('sc_d1', 'sc_d2', 'sc_d3', 'sc_d4', 'sc_d5', 'sc_d6', 'sc_d7', 'sc_d8', 'sc_d9', 'sc_d10', 'sc_d11', 'sc_d12', 'sc_d13', 'sc_d14', 'sc_d15', 'sc_d16', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
                                                      to = c('Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Cross-sectional', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal', 'Longitudinal'))
expt.obj@meta.data$Cohort <- factor(expt.obj@meta.data$Cohort, levels = c('Cross-sectional', 'Longitudinal'))

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


umap_seurat_clusters <- DimPlot(expt.obj, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_age_group <- DimPlot(expt.obj, reduction = "umap", group.by = "AgeGroup", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group, paste('./plots/', experiment, '_umap_age_group', sep = ''))

umap_age_group_2 <- DimPlot(expt.obj, reduction = "umap", group.by = "AgeGroup2", shuffle = TRUE, seed = 123)
generate_figs(umap_age_group_2, paste('./plots/', experiment, '_umap_age_group_2', sep = ''))

umap_sex <- DimPlot(expt.obj, reduction = "umap", group.by = "Sex", shuffle = TRUE, seed = 123)
generate_figs(umap_sex, paste('./plots/', experiment, '_umap_sex', sep = ''))

umap_visit <- DimPlot(expt.obj, reduction = "umap", group.by = "visit", shuffle = TRUE, seed = 123)
generate_figs(umap_visit, paste('./plots/', experiment, '_umap_visit', sep = ''))

umap_code <- DimPlot(expt.obj, reduction = "umap", group.by = "Code", shuffle = TRUE, seed = 123)
generate_figs(umap_code, paste('./plots/', experiment, '_umap_code', sep = ''))

umap_cohort <- DimPlot(expt.obj, reduction = "umap", group.by = "Cohort", shuffle = TRUE, seed = 123)
generate_figs(umap_cohort, paste('./plots/', experiment, '_umap_cohort', sep = ''))



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
azimuth.obj <- subset(azimuth.obj, idents = c("CD8 T", "other T"))
# downsample
azimuth.obj <- SetIdent(azimuth.obj, value = "celltype.l2")
azimuth.obj <- subset(azimuth.obj, downsample = 100)
azimuth_assay <- LayerData(azimuth.obj)

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html#subsetting-summarizedexperiment
monaco.obj <- MonacoImmuneData(ensembl=F)
monaco.obj.md <- monaco.obj@colData %>% as.data.table
monaco.obj.md[, .N, by = c("label.main", "label.fine")]
monaco.index <- monaco.obj$label.main %in% c('CD8+ T cells', 'T cells')
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

for (cell_ref in c("azimuth", "monaco", "dice")){
  for (age_group in unique(md$AgeGroup2)) {
    print(cell_ref)
    print(age_group)
    print(md[AgeGroup2 == age_group, .N, by = cell_ref])
    md_temp <- md[AgeGroup2 == age_group, .N, by = cell_ref]
    setkey(md_temp[, cell_ref := cell_ref], cell_ref)
    setorder(md_temp, cols = cell_ref)
    md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
    md_temp$label <- paste(md_temp$cell_ref, " (", md_temp$percent,"%)", sep = "")
    md_temp$cols <- hcl.colors(n = nrow(md_temp), palette = "Temps")
    barplot_temp <- ggplot(data= md_temp, aes(x=cell_ref, y=percent)) + geom_bar(stat="identity", fill = md_temp$cols)
    generate_figs(barplot_temp, paste('./plots/', experiment, '_barplot_', cell_ref, '_', age_group, sep = ''))
  }
}

md_temp <- md[, .N, by = c('azimuth', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_seurat_clusters, paste('./plots/', experiment, '_barplot_azimuth_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_seurat_clusters, paste('./plots/', experiment, '_barplot_monaco_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('dice', 'seurat_clusters')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_seurat_clusters <- ggplot(md_temp, aes(x = seurat_clusters, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_seurat_clusters, paste('./plots/', experiment, '_barplot_dice_seurat_clusters', sep = ''))

md_temp <- md[, .N, by = c('azimuth', 'AgeGroup2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_azimuth_age_group <- ggplot(md_temp, aes(x = AgeGroup2, y = N, fill = azimuth)) + geom_col(position = "fill")
generate_figs(barplot_azimuth_age_group, paste('./plots/', experiment, '_barplot_azimuth_age_group', sep = ''))

md_temp <- md[, .N, by = c('monaco', 'AgeGroup2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_monaco_age_group <- ggplot(md_temp, aes(x = AgeGroup2, y = N, fill = monaco)) + geom_col(position = "fill")
generate_figs(barplot_monaco_age_group, paste('./plots/', experiment, '_barplot_monaco_age_group', sep = ''))

md_temp <- md[, .N, by = c('dice', 'AgeGroup2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_dice_age_group <- ggplot(md_temp, aes(x = AgeGroup2, y = N, fill = dice)) + geom_col(position = "fill")
generate_figs(barplot_dice_age_group, paste('./plots/', experiment, '_barplot_dice_age_group', sep = ''))

md_temp <- md[, .N, by = c('Phase', 'AgeGroup2')]
md_temp$percent <- round(100*md_temp$N / sum(md_temp$N), digits = 1)
barplot_phase_age_group <- ggplot(md_temp, aes(x = AgeGroup2, y = N, fill = Phase)) + geom_col(position = "fill")
generate_figs(barplot_phase_age_group, paste('./plots/', experiment, '_barplot_phase_age_group', sep = ''))




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

# expt.obj <- SetIdent(expt.obj, value = 'CD8Tref_1')
# expt.obj <- subset(expt.obj, idents = 'CD8 T cell')

# saveRDS(expt.obj, file = paste('./data/', experiment, '_CD8.rds', sep = ''))


expt.obj <- SetIdent(expt.obj, value = 'Cohort')
expt.list <- SplitObject(expt.obj, split.by = "Cohort")

saveRDS(expt.list$`Cross-sectional`, file = paste('./data/', experiment, '_annotated_cross.rds', sep = ''))
saveRDS(expt.list$`Longitudinal`, file = paste('./data/', experiment, '_annotated_long.rds', sep = ''))



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

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv("../../weights/cartex-200-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_200_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_200_weights[match(common, rownames(cartex_200_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_200 <- Z(scores)
expt.obj@meta.data$CARTEx_200i <- integerize(expt.obj@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_84_weights), rownames(expt.obj))
expr <- t(as.matrix(GetAssayData(expt.obj))[match(common, rownames(as.matrix(GetAssayData(expt.obj)))),])
weights <- cartex_84_weights[match(common, rownames(cartex_84_weights)),]
scores <- expr %*% as.matrix(weights)
expt.obj@meta.data$CARTEx_84 <- Z(scores)
expt.obj@meta.data$CARTEx_84i <- integerize(expt.obj@meta.data$CARTEx_84)


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

saveRDS(expt.obj, file = paste('./data/', experiment, '_scored.rds', sep = ''))

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

head(expt.obj)


# UMAP of scores
umap_CARTEx_630i <- DimPlot(expt.obj, group.by = "CARTEx_630i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_630i, paste('./plots/', experiment, '_umap_CARTEx_630i', sep = ''))

umap_CARTEx_200i <- DimPlot(expt.obj, group.by = "CARTEx_200i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_200i, paste('./plots/', experiment, '_umap_CARTEx_200i', sep = ''))

umap_CARTEx_84i <- DimPlot(expt.obj, group.by = "CARTEx_84i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_84i, paste('./plots/', experiment, '_umap_CARTEx_84i', sep = ''))

umap_sig_activationi <- DimPlot(expt.obj, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_umap_sig_activationi', sep = ''))

umap_sig_anergyi <- DimPlot(expt.obj, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_umap_sig_anergyi', sep = ''))

umap_sig_stemnessi <- DimPlot(expt.obj, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_umap_sig_stemnessi', sep = ''))

umap_sig_senescencei <- DimPlot(expt.obj, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_umap_sig_senescencei', sep = ''))

# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(expt.obj, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(expt.obj, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(expt.obj, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(expt.obj, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(expt.obj, features = c("Senescence"), order = TRUE) + fix.sc

generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_umap_sig_activation', sep = ''))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_umap_sig_anergy', sep = ''))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_umap_sig_stemness', sep = ''))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_umap_sig_senescence', sep = ''))



# quality control 

scatter_CARTEx_630_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_630", feature2 = "CARTEx_630_countsproportion")
scatter_CARTEx_200_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_200", feature2 = "CARTEx_200_countsproportion")
scatter_CARTEx_84_countsproportion <- FeatureScatter(expt.obj, feature1 = "CARTEx_84", feature2 = "CARTEx_84_countsproportion")

generate_figs(scatter_CARTEx_630_countsproportion, paste('./plots/', experiment, '_scatter_CARTEx_630_countsproportion', sep = ''))
generate_figs(scatter_CARTEx_200_countsproportion, paste('./plots/', experiment, '_scatter_CARTEx_200_countsproportion', sep = ''))
generate_figs(scatter_CARTEx_84_countsproportion, paste('./plots/', experiment, '_scatter_CARTEx_84_countsproportion', sep = ''))




expt.obj <- SetIdent(expt.obj, value = "monaco")
levels(expt.obj) <- c('Naive CD8 T cells', 'Effector memory CD8 T cells', 'Central memory CD8 T cells', 'Terminal effector CD8 T cells')

vlnplot_CARTEx_84_monaco <- VlnPlot(expt.obj, features = c("CARTEx_84"), group.by = 'monaco', pt.size = 0) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_monaco, paste('./plots/', experiment, '_vlnplot_CARTEx_84_monaco', sep = ''))






