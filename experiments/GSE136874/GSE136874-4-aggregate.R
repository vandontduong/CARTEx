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
  xlab("CAR T cells and controls")
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(6,5))


# compare PD-1 levels
md$PDCD1 <- query.obj.agg@assays$RNA$data['PDCD1',]

aggplot_PDCD1 <- md %>% ggplot(aes(identifier2, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("dodgerblue", "indianred", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','GD2')), label = "p.signif", label.y = 0.15) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CD19','YoungNaive')), label = "p.signif", label.y = 0.18) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GD2','OldTerminal')), label = "p.signif", label.y = 0.21) +
  xlab("CAR T cells and controls")
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_query_agg_aggplot_PDCD1', sep = ''), c(6,5))


