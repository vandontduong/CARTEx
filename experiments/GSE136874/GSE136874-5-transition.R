
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
dm <- readRDS(paste('./data/', experiment, '_dmap.rds', sep = ''))
dpt <- DPT(dm)

expt.obj$dpt <- dpt$dpt

# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

ggplot(md, aes(x = dpt, y = CAR, colour = CAR)) +
  geom_quasirandom(groupOnX = FALSE)

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
transitplot_CAR_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_84', sep = ''), c(8, 5))

transitplot_CAR_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_200', sep = ''), c(8, 5))

transitplot_CAR_CARTEx_630 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_630)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_CARTEx_630, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_630', sep = ''), c(8, 5))

transitplot_CAR_activation <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Activation)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_activation, paste('./plots/', experiment, '_transition_transitplot_CAR_activation', sep = ''), c(8, 5))

transitplot_CAR_anergy <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Anergy)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_anergy, paste('./plots/', experiment, '_transition_transitplot_CAR_anergy', sep = ''), c(8, 5))

transitplot_CAR_senescence <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Senescence)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_senescence, paste('./plots/', experiment, '_transition_transitplot_CAR_senescence', sep = ''), c(8, 5))

transitplot_CAR_stemness <- ggplot(md, aes(x = DC1rank, y = CAR, colour = Stemness)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_stemness, paste('./plots/', experiment, '_transition_transitplot_CAR_stemness', sep = ''), c(8, 5))

transitplot_CAR_LCMV_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = LCMV_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_LCMV_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_LCMV_Tex', sep = ''), c(8, 5))

transitplot_CAR_NKlike_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = NKlike_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_NKlike_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_NKlike_Tex', sep = ''), c(8, 5))

transitplot_CAR_BBD_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = BBD_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_BBD_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_BBD_Tex', sep = ''), c(8, 5))

transitplot_CAR_PD1_Tex <- ggplot(md, aes(x = DC1rank, y = CAR, colour = PD1_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_PD1_Tex, paste('./plots/', experiment, '_transition_transitplot_CAR_PD1_Tex', sep = ''), c(8, 5))


###



expt.obj$UMAP1 <- expt.obj[['umap']]@cell.embeddings[,1]
test <- ggplot(md, aes(x = UMAP1, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc




### exploded volcano

cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)


# Compare GD2 vs CD19
de_genes <- FindMarkers(expt.obj, ident.1 = "GD2", ident.2 = "CD19", group.by = "orig.ident", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2)
signif <- signif[rownames(signif) %in% rownames(cartex_84_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_GD2vCD19 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         selectLab = rownames(signif), drawConnectors = TRUE,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_GD2vCD19, paste('./plots/', experiment, '_transition_volcano_GD2vCD19', sep = ''), c(10, 8))


# examine cells with DC1rank > 1200
table(expt.obj$DC1rank > 1200)

expt.obj@meta.data$DC1rank_g1200 <- expt.obj$DC1rank > 1200

de_genes <- FindMarkers(expt.obj, ident.1 = "TRUE", ident.2 = "FALSE", group.by = "DC1rank_g1200", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2)
signif <- signif[rownames(signif) %in% rownames(cartex_84_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_DC1rank_g1200 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                         selectLab = rownames(signif), drawConnectors = TRUE,
                                         xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_DC1rank_g1200, paste('./plots/', experiment, '_transition_volcano_DC1rank_g1200', sep = ''), c(10, 8))


# AND gate
table(expt.obj$DC1rank > 1200 & expt.obj$CAR == "GD2")

expt.obj@meta.data$DC1rank_g1200_aGD2 <- expt.obj$DC1rank > 1200 & expt.obj$CAR == "GD2"

de_genes <- FindMarkers(expt.obj, ident.1 = "TRUE", ident.2 = "FALSE", group.by = "DC1rank_g1200_aGD2", min.pct = 0.25)
log2fc_lim <- min(ceiling(max(abs(de_genes$avg_log2FC[which(!is.infinite(de_genes$avg_log2FC))]))), 10)
head(de_genes)
signif <- subset(de_genes, p_val < 10e-6 & abs(avg_log2FC) > 2)
signif <- signif[rownames(signif) %in% rownames(cartex_84_weights),]

# change 'log2FoldChange' to 'avg_log2FC' and 'pvalue' to 'p_val'
plot_volcano_DC1rank_g1200_aGD2 <- EnhancedVolcano(de_genes, lab = rownames(de_genes), x = 'avg_log2FC', y = 'p_val', 
                                              selectLab = rownames(signif), drawConnectors = TRUE,
                                              xlim = c(-log2fc_lim, log2fc_lim), labSize = 4.0) # + coord_flip()

generate_figs(plot_volcano_DC1rank_g1200_aGD2, paste('./plots/', experiment, '_transition_volcano_DC1rank_g1200_aGD2', sep = ''), c(10, 8))



