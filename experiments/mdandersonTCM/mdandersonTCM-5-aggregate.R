#!/usr/bin/env Rscript

### Script name: mdandersonTCM-5-aggregate.R
### Description: aggregate seurat object for mdandersonTCM
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'mdandersonTCM'
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


####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)


query.obj.agg$identifier2 <- factor(query.obj.agg$identifier2, levels = c("Healthy donor", "Primary tumor tissue", "Metastatic tumor tissue", "Uninvolved normal tissue", "YoungNaive", "OldTerminal"))

saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition


md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier2")]

glimpse(md)
table(md$identifier2)

aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "pink", "salmon", "skyblue","royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','Primary tumor tissue')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','YoungNaive')), label = "p.signif", label.y = 2.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','OldTerminal')), label = "p.signif", label.y = 2.6) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Metastatic tumor tissue','Primary tumor tissue')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Metastatic tumor tissue','Uninvolved normal tissue')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Primary tumor tissue','OldTerminal')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Uninvolved normal tissue','OldTerminal')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("H", "P", "M", "U", "YN", "OT"))
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(3.5,3))



aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier2, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "pink", "salmon", "skyblue","royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','Primary tumor tissue')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','YoungNaive')), label = "p.signif", label.y = 2.2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Healthy donor','OldTerminal')), label = "p.signif", label.y = 2.6) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Metastatic tumor tissue','Primary tumor tissue')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Metastatic tumor tissue','Uninvolved normal tissue')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Primary tumor tissue','OldTerminal')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Uninvolved normal tissue','OldTerminal')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-2, 3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("H", "P", "M", "U", "YN", "OT"))
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(3.5,3))











