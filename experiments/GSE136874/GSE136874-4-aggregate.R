#!/usr/bin/env Rscript

### Script name: GSE136874-4-aggregate.R
### Description: aggregate seurat object for GSE136874
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/sc-mode.html#pseudo-bulk-aggregation
# https://rdrr.io/bioc/SingleR/man/aggregateReference.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html


query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

query.obj@meta.data$pblabels <- PseudoBulkLabels(query.obj, 5)
md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("identifier2", "monaco")]
md[, .N, by = c("identifier2", "pblabels")]
md[, .N, by = c("identifier2", "monaco", "pblabels")]

unique(query.obj@meta.data$identifier2)
query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'monaco', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'monaco'), return.seurat = TRUE)


####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
query.obj.agg@meta.data$CARTEx_630 <- Z(SignatureScore(query.obj.agg, cartex_630_weights))
query.obj.agg@meta.data$CARTEx_630i <- integerize(query.obj.agg@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
query.obj.agg@meta.data$CARTEx_200 <- Z(SignatureScore(query.obj.agg, cartex_200_weights))
query.obj.agg@meta.data$CARTEx_200i <- integerize(query.obj.agg@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
query.obj.agg@meta.data$CARTEx_84 <- Z(SignatureScore(query.obj.agg, cartex_84_weights))
query.obj.agg@meta.data$CARTEx_84i <- integerize(query.obj.agg@meta.data$CARTEx_84)

####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = TRUE, row.names = 1))
senescence.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = TRUE, row.names = 1))

query.obj.agg <- AddModuleScore(query.obj.agg, features = list(activation.sig, anergy.sig, stemness.sig, senescence.sig), name="State", search = TRUE)

# z score normalization
query.obj.agg@meta.data$Activation <- scale(query.obj.agg@meta.data$State1)
query.obj.agg@meta.data$Anergy <- scale(query.obj.agg@meta.data$State2)
query.obj.agg@meta.data$Stemness <- scale(query.obj.agg@meta.data$State3)
query.obj.agg@meta.data$Senescence <- scale(query.obj.agg@meta.data$State4)

query.obj.agg@meta.data$Activationi <- integerize(query.obj.agg@meta.data$Activation)
query.obj.agg@meta.data$Anergyi <- integerize(query.obj.agg@meta.data$Anergy)
query.obj.agg@meta.data$Stemnessi <- integerize(query.obj.agg@meta.data$Stemness)
query.obj.agg@meta.data$Senescencei <- integerize(query.obj.agg@meta.data$Senescence)

query.obj.agg@meta.data$Activationi <- factor(query.obj.agg@meta.data$Activationi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
query.obj.agg@meta.data$Anergyi <- factor(query.obj.agg@meta.data$Anergyi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
query.obj.agg@meta.data$Stemnessi <- factor(query.obj.agg@meta.data$Stemnessi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
query.obj.agg@meta.data$Senescencei <- factor(query.obj.agg@meta.data$Senescencei, levels = c(4,3,2,1,0,-1,-2,-3,-4))

query.obj.agg@meta.data$State1 <- NULL
query.obj.agg@meta.data$State2 <- NULL
query.obj.agg@meta.data$State3 <- NULL
query.obj.agg@meta.data$State4 <- NULL

# examine other signatures
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))

query.obj.agg <- AddModuleScore(query.obj.agg, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

query.obj.agg@meta.data$NKlike_Tex <- scale(query.obj.agg@meta.data$Signature1)
query.obj.agg@meta.data$LCMV_Tex <- scale(query.obj.agg@meta.data$Signature2)
query.obj.agg@meta.data$BBD_Tex <- scale(query.obj.agg@meta.data$Signature3)
query.obj.agg@meta.data$PD1_Tex <- scale(query.obj.agg@meta.data$Signature4)

query.obj.agg@meta.data$NKlike_Texi <- integerize(query.obj.agg@meta.data$NKlike_Tex)
query.obj.agg@meta.data$LCMV_Texi <- integerize(query.obj.agg@meta.data$LCMV_Tex)
query.obj.agg@meta.data$BBD_Texi <- integerize(query.obj.agg@meta.data$BBD_Tex)
query.obj.agg@meta.data$PD1_Texi <- integerize(query.obj.agg@meta.data$PD1_Tex)

query.obj.agg@meta.data$Signature1 <- NULL
query.obj.agg@meta.data$Signature2 <- NULL
query.obj.agg@meta.data$Signature3 <- NULL
query.obj.agg@meta.data$Signature4 <- NULL


table(query.obj.agg$identifier2)

query.obj.agg$identifier2 <- factor(query.obj.agg$identifier2, levels = c("CD19", "GD2", "YoungNaive", "OldTerminal"))

saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition


md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier2")]

glimpse(md)

md_mean_values <- md %>% group_by(identifier2) %>% summarise(avg = mean(CARTEx_84), stdev = sd(CARTEx_84))

aggplot_qk_CARTEx_84 <- md_mean_values %>% ggplot(aes(identifier2, avg)) +
  geom_col(aes(fill = identifier2), color = "black", width = 0.85) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#22292F", width = 0.1)
generate_figs(aggplot_qk_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_qk_CARTEx_84', sep = ''), c(6,5))


aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier2, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 0.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(5.5,5))

aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.6) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(5.5,5))

aggplot_CARTEx_630 <- md %>% ggplot(aes(identifier2, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.6) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_630, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630', sep = ''), c(5.5,5))



aggplot_activation <- md %>% ggplot(aes(identifier2, Activation)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_activation, paste('./plots/', experiment, '_query_agg_aggplot_activation', sep = ''), c(5.5,5))

aggplot_anergy <- md %>% ggplot(aes(identifier2, Anergy)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_anergy, paste('./plots/', experiment, '_query_agg_aggplot_anergy', sep = ''), c(5.5,5))

aggplot_senescence <- md %>% ggplot(aes(identifier2, Senescence)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_senescence, paste('./plots/', experiment, '_query_agg_aggplot_senescence', sep = ''), c(5.5,5))


aggplot_stemness <- md %>% ggplot(aes(identifier2, Stemness)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_stemness, paste('./plots/', experiment, '_query_agg_aggplot_stemness', sep = ''), c(5.5,5))


aggplot_NKlike_Tex <- md %>% ggplot(aes(identifier2, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.6) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_NKlike_Tex, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex', sep = ''), c(5.5,5))

aggplot_LCMV_Tex <- md %>% ggplot(aes(identifier2, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.6) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_LCMV_Tex, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex', sep = ''), c(5.5,5))

aggplot_BBD_Tex <- md %>% ggplot(aes(identifier2, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.3) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.6) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_BBD_Tex, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex', sep = ''), c(5.5,5))

aggplot_PD1_Tex <- md %>% ggplot(aes(identifier2, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 1.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 1.8) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-2, 2)
generate_figs(aggplot_PD1_Tex, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex', sep = ''), c(5.5,5))



# compare PD-1 levels
md$PDCD1 <- query.obj.agg@assays$RNA$data['PDCD1',]

aggplot_PDCD1 <- md %>% ggplot(aes(identifier2, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 0.15) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 0.18) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 0.21) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_query_agg_aggplot_PDCD1', sep = ''), c(6,5))


