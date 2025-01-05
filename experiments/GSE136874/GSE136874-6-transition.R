
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))


# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

ggplot(md, aes(x = dpt, y = CAR, colour = CAR)) +
  geom_quasirandom(groupOnX = FALSE)

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))


transitplot_CAR_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_84', sep = ''), c(3.5,2))

transitplot_CAR_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_200', sep = ''), c(3.5,2))

transitplot_CAR_CARTEx_630 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_630)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_CARTEx_630, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_630', sep = ''), c(3.5,2))

transitplot_CAR_activation <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Activation)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_activation, paste('./plots/', experiment, '_transition_transitplot_CAR_activation', sep = ''), c(3.5,2))

transitplot_CAR_anergy <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Anergy)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_anergy, paste('./plots/', experiment, '_transition_transitplot_CAR_anergy', sep = ''), c(3.5,2))

transitplot_CAR_senescence <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Senescence)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_senescence, paste('./plots/', experiment, '_transition_transitplot_CAR_senescence', sep = ''), c(3.5,2))

transitplot_CAR_stemness <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Stemness)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_stemness, paste('./plots/', experiment, '_transition_transitplot_CAR_stemness', sep = ''), c(3.5,2))

transitplot_CAR_LCMV_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = LCMV_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_LCMV_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_LCMV_Tex', sep = ''), c(3.5,2))

transitplot_CAR_NKlike_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = NKlike_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_NKlike_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_NKlike_Tex', sep = ''), c(3.5,2))

transitplot_CAR_BBD_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = BBD_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_BBD_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_BBD_Tex', sep = ''), c(3.5,2))

transitplot_CAR_PD1_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = PD1_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 1500)
generate_figs(transitplot_CAR_PD1_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_PD1_Tex', sep = ''), c(3.5,2))


###



expt.obj$UMAP1 <- expt.obj[['umap']]@cell.embeddings[,1]
test <- ggplot(md, aes(x = UMAP1, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc




umap_DC1rank <- FeaturePlot(expt.obj, feature = 'DC1rank', order = FALSE, pt.size = 0.1) + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c("yellow", "orange", "red", "violetred", "darkorchid"))
generate_figs(umap_DC1rank, paste('./plots/', experiment, '_prepare_umap_DC1rank', sep = ''), c(2.9, 2))


### 





### exploded volcano
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)


# Compare GD2 vs CD19
de_genes <- FindMarkers(expt.obj, ident.1 = "GD2", ident.2 = "CD19", group.by = "orig.ident", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_GD2vCD19 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 20, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_GD2vCD19, paste('./plots/', experiment, '_transition_volcano_GD2vCD19', sep = ''), c(6, 5))



signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_200_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_200_weights), "CARTEx")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_GD2vCD19_200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 20, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_GD2vCD19_200, paste('./plots/', experiment, '_transition_volcano_GD2vCD19_200', sep = ''), c(6, 5))








# examine cells with DC1rank > 1200
table(expt.obj$DC1rank > 1200)

expt.obj@meta.data$DC1rank_g1200 <- expt.obj$DC1rank > 1200

de_genes <- FindMarkers(expt.obj, ident.1 = "TRUE", ident.2 = "FALSE", group.by = "DC1rank_g1200", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_DC1rank_g1200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                              selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                              shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_DC1rank_g1200, paste('./plots/', experiment, '_transition_volcano_DC1rank_g1200', sep = ''), c(6, 5))


# AND gate
table(expt.obj$DC1rank > 1200 & expt.obj$CAR == "GD2")

expt.obj@meta.data$DC1rank_g1200_aGD2 <- expt.obj$DC1rank > 1200 & expt.obj$CAR == "GD2"

de_genes <- FindMarkers(expt.obj, ident.1 = "TRUE", ident.2 = "FALSE", group.by = "DC1rank_g1200_aGD2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% rownames(cartex_630_weights),]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, rownames(cartex_630_weights), "C5")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_DC1rank_g1200_aGD2 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                                   pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                                   selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                                   shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                                   xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_DC1rank_g1200_aGD2, paste('./plots/', experiment, '_transition_volcano_DC1rank_g1200_aGD2', sep = ''), c(6, 5))





# examine CARTEx C2 genes
cartex_C2 <- rownames(read.csv(paste(PATH_CONSTRUCTION, "./data/cartex-cluster-2.csv", sep = ''), header = TRUE, row.names = 1))

# Compare GD2 vs CD19
de_genes <- FindMarkers(expt.obj, ident.1 = "GD2", ident.2 = "CD19", group.by = "orig.ident", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 0.5)
signif <- signif[rownames(signif) %in% cartex_C2,]

# create custom key-value pairs for CARTEx genes
keyvals <- CustomKeyValPairsVolcanoPlot(de_genes, cartex_C2, "C2")

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_GD2vCD19_C2genes <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         pCutoff = 10e-6, FCcutoff = 0.5, title = NULL, subtitle = NULL,
                                         selectLab = rownames(signif), drawConnectors = TRUE, typeConnectors = 'closed', endsConnectors = 'last', directionConnectors = 'both', colConnectors = 'black', max.overlaps = 15, 
                                         shapeCustom = keyvals$shape, colAlpha = 0.75, pointSize = keyvals$ptsize,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) + theme_classic() + theme(legend.position = "top", legend.title=element_blank()) # + coord_flip()

generate_figs(plot_volcano_GD2vCD19_C2genes, paste('./plots/', experiment, '_transition_volcano_GD2vCD19_C2genes', sep = ''), c(6, 5))





