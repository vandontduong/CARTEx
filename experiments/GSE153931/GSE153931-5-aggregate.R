#!/usr/bin/env Rscript

### Script name: GSE125881-5-aggregate.R
### Description: aggregate seurat object for GSE125881
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE153931'
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
md[, .N, by = c("identifier5", "monaco")]
md[, .N, by = c("identifier5", "pblabels")]
md[, .N, by = c("identifier5", "monaco", "pblabels")]

unique(query.obj@meta.data$identifier5)
query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier5', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier5', 'monaco', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier5', 'monaco'), return.seurat = TRUE)

Idents(query.obj.agg) <- "identifier5"
# query.obj.agg <- subset(query.obj.agg, idents = "void", invert = TRUE)
# query.obj.agg$identifier4 <- factor(query.obj.agg$identifier4, levels = c('Mild', 'Moderate', 'Severe', 'YoungNaive', 'OldTerminal'))

query.obj.agg$identifier5 <- factor(query.obj.agg$identifier5, levels = c('Mild', 'Severe', 'YoungNaive', 'OldTerminal'))


####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)



# query.obj.agg$identifier5 <- factor(query.obj.agg$identifier5, levels = c('Mild', 'Moderate', 'Severe', 'void', 'YoungNaive', 'OldTerminal'))

saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition


md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier5")]

glimpse(md)


aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier5, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(3,3))


aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier5, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(3,3))

aggplot_CARTEx_630 <- md %>% ggplot(aes(identifier5, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_CARTEx_630, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630', sep = ''), c(3,3))





aggplot_activation <- md %>% ggplot(aes(identifier5, Activation)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("Activation") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_activation, paste('./plots/', experiment, '_query_agg_aggplot_activation', sep = ''), c(3,3))

aggplot_anergy <- md %>% ggplot(aes(identifier5, Anergy)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("Anergy") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_anergy, paste('./plots/', experiment, '_query_agg_aggplot_anergy', sep = ''), c(3,3))

aggplot_senescence <- md %>% ggplot(aes(identifier5, Senescence)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("Senescence") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_senescence, paste('./plots/', experiment, '_query_agg_aggplot_senescence', sep = ''), c(3,3))

aggplot_stemness <- md %>% ggplot(aes(identifier5, Stemness)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("Stemness") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_stemness, paste('./plots/', experiment, '_query_agg_aggplot_stemness', sep = ''), c(3,3))




aggplot_NKlike_Tex <- md %>% ggplot(aes(identifier5, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_NKlike_Tex, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex', sep = ''), c(3,3))


aggplot_LCMV_Tex <- md %>% ggplot(aes(identifier5, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_LCMV_Tex, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex', sep = ''), c(3,3))


aggplot_BBD_Tex <- md %>% ggplot(aes(identifier5, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_BBD_Tex, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex', sep = ''), c(3,3))


aggplot_PD1_Tex <- md %>% ggplot(aes(identifier5, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c('cadetblue', 'indianred', "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','Severe')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','YoungNaive')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Mild','OldTerminal')), label = "p.signif", label.y = 2.2) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("M", "S", "YN", "OT"))
generate_figs(aggplot_PD1_Tex, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex', sep = ''), c(3,3))






















