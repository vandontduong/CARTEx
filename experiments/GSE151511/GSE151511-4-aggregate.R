#!/usr/bin/env Rscript

### Script name: GSE151511-4-aggregate.R
### Description: aggregate seurat object for GSE151511
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE151511'
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

Idents(query.obj.agg) <- "identifier2"
query.obj.agg <- subset(query.obj.agg, idents = "NE", invert = TRUE)

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


query.obj.agg$identifier2 <- factor(query.obj.agg$identifier2, levels = c("NR", "R", "NE", "YoungNaive", "OldTerminal"))

saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition
md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier2")]

glimpse(md)


aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier2, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(5.5,5))

aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(5.5,5))

aggplot_CARTEx_630 <- md %>% ggplot(aes(identifier2, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_630, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630', sep = ''), c(5.5,5))


aggplot_LCMV_Tex <- md %>% ggplot(aes(identifier2, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_LCMV_Tex, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex', sep = ''), c(5.5,5))

aggplot_NKlike_Tex <- md %>% ggplot(aes(identifier2, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_NKlike_Tex, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex', sep = ''), c(5.5,5))

aggplot_BBD_Tex <- md %>% ggplot(aes(identifier2, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_BBD_Tex, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex', sep = ''), c(5.5,5))

aggplot_PD1_Tex <- md %>% ggplot(aes(identifier2, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_PD1_Tex, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex', sep = ''), c(5.5,5))




####################################################################################################
#################################### Repeat for clinical response ##################################
####################################################################################################


unique(query.obj@meta.data$identifier3_clinical)

query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier3_clinical', 'pblabels'), return.seurat = TRUE)

Idents(query.obj.agg) <- "identifier3_clinical"
query.obj.agg <- subset(query.obj.agg, idents = "NE", invert = TRUE)


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


query.obj.agg$identifier3_clinical <- factor(query.obj.agg$identifier3_clinical, levels = c("CR", "PR", "PD", "YoungNaive", "OldTerminal"))


saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored_clinical.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored_clinical.rds', sep = ''))

head(query.obj.agg)

md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier3_clinical")]

glimpse(md)

aggplot_CARTEx_84_clinical <- md %>% ggplot(aes(identifier3_clinical, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_clinical), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "slateblue", "firebrick", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','PD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('PD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','YoungNaive')), label = "p.signif", label.y = 1.4) +
  xlab("CAR T cells and controls") + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_84_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84_clinical', sep = ''), c(5.5,5))

aggplot_CARTEx_200_clinical <- md %>% ggplot(aes(identifier3_clinical, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_clinical), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "slateblue", "firebrick", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','PD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('PD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','YoungNaive')), label = "p.signif", label.y = 1.4) +
  xlab("CAR T cells and controls") + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_200_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200_clinical', sep = ''), c(5.5,5))

aggplot_CARTEx_630_clinical <- md %>% ggplot(aes(identifier3_clinical, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_clinical), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "slateblue", "firebrick", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','PD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('PD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','YoungNaive')), label = "p.signif", label.y = 1.4) +
  xlab("CAR T cells and controls") + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_630_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630_clinical', sep = ''), c(5.5,5))





####################################################################################################
########################################### Repeat for EMR #########################################
####################################################################################################



unique(query.obj@meta.data$identifier3_EMR)

query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier3_EMR', 'pblabels'), return.seurat = TRUE)

Idents(query.obj.agg) <- "identifier3_EMR"
query.obj.agg <- subset(query.obj.agg, idents = "NE", invert = TRUE)


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


query.obj.agg$identifier3_EMR <- factor(query.obj.agg$identifier3_EMR, levels = c("BAD", "GOOD", "YoungNaive", "OldTerminal"))


saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored_EMR.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored_EMR.rds', sep = ''))

head(query.obj.agg)

md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier3_EMR")]

glimpse(md)

aggplot_CARTEx_84_EMR <- md %>% ggplot(aes(identifier3_EMR, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_84_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84_EMR', sep = ''), c(5.5,5))

aggplot_CARTEx_200_EMR <- md %>% ggplot(aes(identifier3_EMR, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_200_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200_EMR', sep = ''), c(5.5,5))

aggplot_CARTEx_630_EMR <- md %>% ggplot(aes(identifier3_EMR, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_630_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630_EMR', sep = ''), c(5.5,5))

aggplot_LCMV_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_LCMV_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex_EMR', sep = ''), c(5.5,5))

aggplot_NKlike_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_NKlike_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex_EMR', sep = ''), c(5.5,5))

aggplot_BBD_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_BBD_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex_EMR', sep = ''), c(5.5,5))

aggplot_PD1_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_PD1_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex_EMR', sep = ''), c(5.5,5))








