#!/usr/bin/env Rscript

### Script name: GSE160160-4-explore.R
### Description: explore seurat object for GSE160160
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
###################################### Gene expression heatmaps ####################################
####################################################################################################

# expt.obj <- readRDS(paste('./data/', experiment, '.rds', sep = ''))
expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))

allgenes_exCARTEx <- rownames(expt.obj)[!(rownames(expt.obj) %in% rownames(cartex_630_weights))]

expt.obj.downsample <- subset(expt.obj, downsample = 100)

heatmap_allgenes_exposure <- DoHeatmap(expt.obj.downsample, group.by = "exposure")
heatmap_allgenes_exCARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = allgenes_exCARTEx, group.by = "exposure")
heatmap_CARTEx630_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "exposure")
heatmap_CARTEx200_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "exposure")
heatmap_CARTEx84_exposure <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_84_weights), group.by = "exposure")
heatmap_TSG_exposure <- DoHeatmap(expt.obj.downsample, features = TSG, group.by = "exposure")

generate_figs(heatmap_allgenes_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exposure'), c(15, 12))
generate_figs(heatmap_allgenes_exCARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_allgenes_exCARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx630_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx200_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_exposure'), c(15, 12))
generate_figs(heatmap_CARTEx84_exposure, paste0('./plots/', experiment, '_explore_heatmap_CARTEx84_exposure'), c(15, 12))
generate_figs(heatmap_TSG_exposure, paste0('./plots/', experiment, '_explore_heatmap_TSG_exposure'), c(15, 12))

heatmap_CARTEx630_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_630_weights), group.by = "monaco")
heatmap_CARTEx200_monaco <- DoHeatmap(expt.obj.downsample, features = rownames(cartex_200_weights), group.by = "monaco")
generate_figs(heatmap_CARTEx630_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx630_monaco'), c(15, 12))
generate_figs(heatmap_CARTEx200_monaco, paste0('./plots/', experiment, '_explore_heatmap_CARTEx200_monaco'), c(15, 12))




DimPlot(expt.obj, group.by = 'exposure')
FeaturePlot(expt.obj, features = c('HNF1A'))
FeaturePlot(expt.obj, features = c('HNF1B'))

# https://www.mcponline.org/article/S1535-9476(20)30681-2/fulltext
SYK_substrates <- FeaturePlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'))
generate_figs(SYK_substrates, paste0('./plots/', experiment, '_explore_SYK_substrates'), c(15, 15))

vlnplot_SYK_substrates <- VlnPlot(expt.obj, features = c('SYK', 'ARHGDIA', 'BTK', 'DBNL', 'DCTN2', 'FKBP8', 'GCET2', 'HCA59', 'HCLS1', 'HIRIP3', 'MAP1B', 'NFYA', 'PAG1', 'PIH1D1', 'PLCG2', 'PTPN11', 'RTN4', 'SF3A3', 'SFPQ', 'SMAP', 'STIP1', 'TBCB', 'TOMM34', 'TUBB', 'TUBB4B', 'UBA1'), group.by = 'exposure')
generate_figs(vlnplot_SYK_substrates, paste0('./plots/', experiment, '_explore_vlnplot_SYK_substrates'), c(15, 15))


# Compare CAE to day0

cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

de_genes <- FindMarkers(expt.obj, ident.1 = "CAE", ident.2 = "day0", group.by = "exposure", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAEvday0 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         selectLab = rownames(signif), drawConnectors = TRUE,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_CAEvday0, paste('./plots/', experiment, '_explore_volcano_CAEvday0', sep = ''), c(10, 8))


FeaturePlot(expt.obj, features = c('CCL3', 'CCL4', 'PPARG', 'METRNL', 'IFNG', 'IL2RA', 'PMCH', 'EPAS1', 'RDH10'))


# Compare CAE exhausted to CAE not exhausted
# 3, 1, 8, 0 vs 6, 9, 7

de_genes <- FindMarkers(expt.obj, ident.1 = c(0,1,3,8), ident.2 = c(6,7,9), group.by = "seurat_clusters", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_CAE2groups <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         selectLab = rownames(signif), drawConnectors = TRUE,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_CAE2groups, paste('./plots/', experiment, '_explore_volcano_CAE2groups', sep = ''), c(10, 8))

de_genes_sorted_CAE2groups <- arrange(filter(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2), avg_log2FC)
write.csv(de_genes_sorted_CAE2groups, "data/de_genes_sorted_CAE2groups.csv")


FeaturePlot(expt.obj, features = rownames(filter(de_genes_sorted_CAE2groups, avg_log2FC < 0)))





