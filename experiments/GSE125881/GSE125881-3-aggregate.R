### Script name: GSE125881-3-aggregate.R
### Description: aggregate seurat object for GSE125881
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/sc-mode.html#pseudo-bulk-aggregation
# https://rdrr.io/bioc/SingleR/man/aggregateReference.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

# https://github.com/satijalab/seurat/issues/6007



# pseudo_samples <- 3
# num_cells <- length(Cells(expt.obj))
# labels <- sample.int(pseudo_samples, num_cells, replace = TRUE)

query.obj@meta.data$pblabels <- PseudoBulkLabels(query.obj, 5)
md <- query.obj@meta.data %>% as.data.table
md[, .N, by = c("identifier2", "monaco")]
md[, .N, by = c("identifier2", "pblabels")]
md[, .N, by = c("identifier2", "monaco", "pblabels")]

unique(query.obj@meta.data$identifier2)
query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'monaco', 'pblabels'), return.seurat = TRUE)
# query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier2', 'monaco'), return.seurat = TRUE)

# filter by cell type?
# use identifier instead of Group when aggregating query data


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


saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)


# CARTEx violin plot

query.obj.agg$identifier2 <- factor(query.obj.agg$identifier2, levels = c("IP", "ExpansionPeak", "Contraction", "Late", "YoungNaive", "OldTerminal"))
vlnplot_CARTEx_84 <- VlnPlot(query.obj.agg, features = c("CARTEx_84"), group.by = 'identifier2', y.max = 6, pt.size = 2) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','ExpansionPeak')), label = "p.signif", label.y = 4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Late')), label = "p.signif", label.y = 4.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','YoungNaive')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('ExpansionPeak','YoungNaive')), label = "p.signif", label.y = 5.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Late','YoungNaive')), label = "p.signif", label.y = 3.5)
# generate_figs(vlnplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_vlnplot_CARTEx_84', sep = ''))


md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier2")]

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition

glimpse(md)

md_mean_values <- md %>% group_by(identifier2) %>% summarise(avg = mean(CARTEx_84), stdev = sd(CARTEx_84))

md_mean_values %>% 
  ggplot(aes(identifier2, avg)) +
  geom_col(aes(fill = identifier2), color = "black", width = 0.85) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#22292F", width = 0.1)

md %>% ggplot(aes(identifier2, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','ExpansionPeak')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','Late')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('IP','YoungNaive')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('ExpansionPeak','YoungNaive')), label = "p.signif", label.y = 2.5) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Late','YoungNaive')), label = "p.signif", label.y = 1)





