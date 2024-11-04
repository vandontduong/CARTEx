#!/usr/bin/env Rscript

### Script name: GSE120575-4-aggregate.R
### Description: aggregate seurat object for GSE120575
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE120575'
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


####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)


unique(query.obj.agg$identifier5)

query.obj.agg$identifier5 <- factor(query.obj.agg$identifier5, levels = c("Pre-NR", "Post-NR", "Pre-R", "Post-R", "YoungNaive", "OldTerminal"))

saveRDS(query.obj.agg, file = paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

# query.obj.agg <- readRDS(paste('./data/', experiment, '_query_agg_scored.rds', sep = ''))

head(query.obj.agg)

# https://ggplot2tutor.com/tutorials/barchart_simple
# https://stackoverflow.com/questions/65675688/ggplot-filling-color-based-on-condition


md <- query.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("identifier5")]

glimpse(md)

md_mean_values <- md %>% group_by(identifier5) %>% summarise(avg = mean(CARTEx_84), stdev = sd(CARTEx_84))

aggplot_qk_CARTEx_84 <- md_mean_values %>% ggplot(aes(identifier5, avg)) +
  geom_col(aes(fill = identifier5), color = "black", width = 0.85) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#22292F", width = 0.1)
generate_figs(aggplot_qk_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_qk_CARTEx_84', sep = ''), c(6,5))


aggplot_CARTEx_84 <- md %>% ggplot(aes(identifier5, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Post-NR')), label = "p.signif", label.y = 0.75) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','Post-R')), label = "p.signif", label.y = 0.75) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.25) +
  xlab("T cells and controls") + geom_point()
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(6,5))



# exclude post- (baseline only)
Idents(query.obj.agg) <- "identifier5"
query.obj.agg <- subset(query.obj.agg, idents = c("Post-NR", "Post-R"), invert = TRUE)


# re-normalize / z-score
query.obj.agg@meta.data$CARTEx_84 <- scale(query.obj.agg@meta.data$CARTEx_84)
query.obj.agg@meta.data$CARTEx_200 <- scale(query.obj.agg@meta.data$CARTEx_200)
query.obj.agg@meta.data$CARTEx_630 <- scale(query.obj.agg@meta.data$CARTEx_630)

query.obj.agg@meta.data$Activation <- scale(query.obj.agg@meta.data$Activation)
query.obj.agg@meta.data$Anergy <- scale(query.obj.agg@meta.data$Anergy)
query.obj.agg@meta.data$Stemness <- scale(query.obj.agg@meta.data$Stemness)
query.obj.agg@meta.data$Senescence <- scale(query.obj.agg@meta.data$Senescence)

query.obj.agg@meta.data$NKlike_Tex <- scale(query.obj.agg@meta.data$NKlike_Tex)
query.obj.agg@meta.data$LCMV_Tex <- scale(query.obj.agg@meta.data$LCMV_Tex)
query.obj.agg@meta.data$BBD_Tex <- scale(query.obj.agg@meta.data$BBD_Tex)
query.obj.agg@meta.data$PD1_Tex <- scale(query.obj.agg@meta.data$PD1_Tex)


md <- query.obj.agg@meta.data %>% as.data.table


md_mean_values <- md %>% group_by(identifier5) %>% summarise(avg = mean(CARTEx_84), stdev = sd(CARTEx_84))

aggplot_qk_CARTEx_84_baseline <- md_mean_values %>% ggplot(aes(identifier5, avg)) +
  geom_col(aes(fill = identifier5), color = "black", width = 0.85) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  geom_errorbar(aes(ymin = avg - stdev, ymax = avg + stdev), color = "#22292F", width = 0.1)
generate_figs(aggplot_qk_CARTEx_84_baseline, paste('./plots/', experiment, '_query_agg_aggplot_qk_CARTEx_84_baseline', sep = ''), c(6,5))


aggplot_CARTEx_84_baseline <- md %>% ggplot(aes(identifier5, CARTEx_84)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_84_baseline, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84_baseline', sep = ''), c(3,3))

aggplot_CARTEx_200_baseline <- md %>% ggplot(aes(identifier5, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_200_baseline, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200_baseline', sep = ''), c(3,3))

aggplot_CARTEx_630_baseline <- md %>% ggplot(aes(identifier5, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_630_baseline, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630_baseline', sep = ''), c(3,3))



aggplot_activation_baseline <- md %>% ggplot(aes(identifier5, Activation)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  xlab("T cells and controls") + geom_point() + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_activation_baseline, paste('./plots/', experiment, '_query_agg_aggplot_activation_baseline', sep = ''), c(3,3))

aggplot_anergy_baseline <- md %>% ggplot(aes(identifier5, Anergy)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  xlab("T cells and controls") + geom_point() + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_anergy_baseline, paste('./plots/', experiment, '_query_agg_aggplot_anergy_baseline', sep = ''), c(3,3))

aggplot_senescence_baseline <- md %>% ggplot(aes(identifier5, Senescence)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  xlab("T cells and controls") + geom_point() + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_senescence_baseline, paste('./plots/', experiment, '_query_agg_aggplot_senescence_baseline', sep = ''), c(3,3))

aggplot_stemness_baseline <- md %>% ggplot(aes(identifier5, Stemness)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  xlab("T cells and controls") + geom_point() + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_stemness_baseline, paste('./plots/', experiment, '_query_agg_aggplot_stemness_baseline', sep = ''), c(3,3))



aggplot_NKlike_Tex_baseline <- md %>% ggplot(aes(identifier5, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_NKlike_Tex_baseline, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex_baseline', sep = ''), c(3,3))

aggplot_LCMV_Tex_baseline <- md %>% ggplot(aes(identifier5, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_LCMV_Tex_baseline, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex_baseline', sep = ''), c(3,3))

aggplot_BBD_Tex_baseline <- md %>% ggplot(aes(identifier5, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_BBD_Tex_baseline, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex_baseline', sep = ''), c(3,3))

aggplot_PD1_Tex_baseline <- md %>% ggplot(aes(identifier5, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier5), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','Pre-R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-R','OldTerminal')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_PD1_Tex_baseline, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex_baseline', sep = ''), c(3,3))











