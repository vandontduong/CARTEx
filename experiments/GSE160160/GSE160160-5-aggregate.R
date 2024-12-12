#!/usr/bin/env Rscript

### Script name: GSE160160-4-aggregate.R
### Description: aggregate seurat object for GSE160160
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
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
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)


query.obj.agg$identifier2 <- factor(query.obj.agg$identifier2, levels = c("day0", "day20", "YoungNaive", "OldTerminal"))

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
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#22292F", width = 0.1) +
  scale_fill_manual(values = c("skyblue", "cadetblue", "royalblue", "orchid"))
generate_figs(aggplot_qk_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_qk_CARTEx_84', sep = ''), c(3,3))


aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier2, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(3,3))


aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(3,3))


aggplot_CARTEx_630 <- md %>% ggplot(aes(identifier2, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_CARTEx_630, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630', sep = ''), c(3,3))




aggplot_activation <- md %>% ggplot(aes(identifier2, Activation)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_activation, paste('./plots/', experiment, '_query_agg_aggplot_activation', sep = ''), c(3,3))

aggplot_anergy <- md %>% ggplot(aes(identifier2, Anergy)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_anergy, paste('./plots/', experiment, '_query_agg_aggplot_anergy', sep = ''), c(3,3))

aggplot_senescence <- md %>% ggplot(aes(identifier2, Senescence)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_senescence, paste('./plots/', experiment, '_query_agg_aggplot_senescence', sep = ''), c(3,3))


aggplot_stemness <- md %>% ggplot(aes(identifier2, Stemness)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point()
generate_figs(aggplot_stemness, paste('./plots/', experiment, '_query_agg_aggplot_stemness', sep = ''), c(3,3))


aggplot_NKlike_Tex <- md %>% ggplot(aes(identifier2, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_NKlike_Tex, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex', sep = ''), c(3,3))

aggplot_LCMV_Tex <- md %>% ggplot(aes(identifier2, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_LCMV_Tex, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex', sep = ''), c(3,3))

aggplot_BBD_Tex <- md %>% ggplot(aes(identifier2, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_BBD_Tex, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex', sep = ''), c(3,3))

aggplot_PD1_Tex <- md %>% ggplot(aes(identifier2, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day20','OldTerminal')), label = "p.signif", label.y = 2) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','YoungNaive')), label = "p.signif", label.y = 2.5) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-1.5,3) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20", "YN", "OT"))
generate_figs(aggplot_PD1_Tex, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex', sep = ''), c(3,3))




# compare checkpoint levels

md$PDCD1 <- query.obj.agg@assays$RNA$data['PDCD1',] # PD-1
md$HAVCR2 <- query.obj.agg@assays$RNA$data['HAVCR2',] # TIM3
md$LAG3 <- query.obj.agg@assays$RNA$data['LAG3',]
md$CTLA4 <- query.obj.agg@assays$RNA$data['CTLA4',]
md$NT5E <- query.obj.agg@assays$RNA$data['NT5E',] # CD73

aggplot_PDCD1 <- md %>% ggplot(aes(identifier2, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point() + ylim(0, 1)
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_query_agg_aggplot_PDCD1', sep = ''), c(6,5))

aggplot_HAVCR2 <- md %>% ggplot(aes(identifier2, HAVCR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point() + ylim(0, 1)
generate_figs(aggplot_HAVCR2, paste('./plots/', experiment, '_query_agg_aggplot_HAVCR2', sep = ''), c(6,5))

aggplot_LAG3 <- md %>% ggplot(aes(identifier2, LAG3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point() + ylim(0, 1)
generate_figs(aggplot_LAG3, paste('./plots/', experiment, '_query_agg_aggplot_LAG3', sep = ''), c(6,5))

aggplot_CTLA4 <- md %>% ggplot(aes(identifier2, CTLA4)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point() + ylim(0, 1)
generate_figs(aggplot_CTLA4, paste('./plots/', experiment, '_query_agg_aggplot_CTLA4', sep = ''), c(6,5))

aggplot_NT5E <- md %>% ggplot(aes(identifier2, NT5E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2)) +
  scale_fill_manual(values=c("skyblue", "cadetblue", "royalblue", "orchid")) +
  xlab("CAR T cells and controls") + geom_point() + ylim(0, 1)
generate_figs(aggplot_NT5E, paste('./plots/', experiment, '_query_agg_aggplot_NT5E', sep = ''), c(6,5))









