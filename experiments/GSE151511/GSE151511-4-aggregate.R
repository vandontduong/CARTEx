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
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)

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
  ylab("CARTEx 84") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_84, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84', sep = ''), c(3,3))

aggplot_CARTEx_200 <- md %>% ggplot(aes(identifier2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(3,3))

aggplot_CARTEx_630 <- md %>% ggplot(aes(identifier2, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_CARTEx_630, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630', sep = ''), c(3,3))


aggplot_LCMV_Tex <- md %>% ggplot(aes(identifier2, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_LCMV_Tex, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex', sep = ''), c(3,3))

aggplot_NKlike_Tex <- md %>% ggplot(aes(identifier2, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_NKlike_Tex, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex', sep = ''), c(3,3))

aggplot_BBD_Tex <- md %>% ggplot(aes(identifier2, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_BBD_Tex, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex', sep = ''), c(3,3))

aggplot_PD1_Tex <- md %>% ggplot(aes(identifier2, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("PD1") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_PD1_Tex, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex', sep = ''), c(3,3))



aggplot_activation <- md %>% ggplot(aes(identifier2, Activation)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("Activation") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_activation, paste('./plots/', experiment, '_query_agg_aggplot_activation', sep = ''), c(3,3))

aggplot_anergy <- md %>% ggplot(aes(identifier2, Anergy)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("Anergy") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_anergy, paste('./plots/', experiment, '_query_agg_aggplot_anergy', sep = ''), c(3,3))


aggplot_senescence <- md %>% ggplot(aes(identifier2, Senescence)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("Senescence") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_senescence, paste('./plots/', experiment, '_query_agg_aggplot_senescence', sep = ''), c(3,3))



aggplot_stemness <- md %>% ggplot(aes(identifier2, Stemness)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','R')), label = "p.signif", label.y = 1.4) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('R','YoungNaive')), label = "p.signif", label.y = 1.1) +
  ylab("Stemness") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("NR", "R", "YN", "OT"))
generate_figs(aggplot_stemness, paste('./plots/', experiment, '_query_agg_aggplot_stemness', sep = ''), c(3,3))








####################################################################################################
#################################### Repeat for clinical response ##################################
####################################################################################################


unique(query.obj@meta.data$identifier3_clinical)

query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier3_clinical', 'pblabels'), return.seurat = TRUE)

Idents(query.obj.agg) <- "identifier3_clinical"
query.obj.agg <- subset(query.obj.agg, idents = "NE", invert = TRUE)


####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)

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
generate_figs(aggplot_CARTEx_84_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84_clinical', sep = ''), c(3.5,3))

aggplot_CARTEx_200_clinical <- md %>% ggplot(aes(identifier3_clinical, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_clinical), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "slateblue", "firebrick", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','PD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('PD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','YoungNaive')), label = "p.signif", label.y = 1.4) +
  xlab("CAR T cells and controls") + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_200_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200_clinical', sep = ''), c(3.5,3))

aggplot_CARTEx_630_clinical <- md %>% ggplot(aes(identifier3_clinical, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_clinical), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("seagreen", "slateblue", "firebrick", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','PD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('PD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('CR','YoungNaive')), label = "p.signif", label.y = 1.4) +
  xlab("CAR T cells and controls") + geom_point() + ylim(-2, 2)
generate_figs(aggplot_CARTEx_630_clinical, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630_clinical', sep = ''), c(3.5,3))





####################################################################################################
########################################### Repeat for EMR #########################################
####################################################################################################



unique(query.obj@meta.data$identifier3_EMR)

query.obj.agg <- AggregateExpression(query.obj, group.by = c('identifier3_EMR', 'pblabels'), return.seurat = TRUE)

Idents(query.obj.agg) <- "identifier3_EMR"
query.obj.agg <- subset(query.obj.agg, idents = "NE", invert = TRUE)


####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

query.obj.agg <- ScoreSubroutine(query.obj.agg)

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
generate_figs(aggplot_CARTEx_84_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_84_EMR', sep = ''), c(3,3))

aggplot_CARTEx_200_EMR <- md %>% ggplot(aes(identifier3_EMR, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_200_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200_EMR', sep = ''), c(3,3))

aggplot_CARTEx_630_EMR <- md %>% ggplot(aes(identifier3_EMR, CARTEx_630)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("CARTEx 630") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_CARTEx_630_EMR, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_630_EMR', sep = ''), c(3,3))

aggplot_LCMV_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, LCMV_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("LCMV") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_LCMV_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_LCMV_Tex_EMR', sep = ''), c(3,3))

aggplot_NKlike_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, NKlike_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("NK-like") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_NKlike_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_NKlike_Tex_EMR', sep = ''), c(3,3))

aggplot_BBD_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, BBD_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_BBD_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_BBD_Tex_EMR', sep = ''), c(3,3))

aggplot_PD1_Tex_EMR <- md %>% ggplot(aes(identifier3_EMR, PD1_Tex)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = identifier3_EMR), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(values = c("firebrick", "seagreen", "royalblue", "orchid")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','BAD')), label = "p.signif", label.y = 1.1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('BAD','OldTerminal')), label = "p.signif", label.y = 1.7) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('GOOD','YoungNaive')), label = "p.signif", label.y = 1.4) +
  ylab("BBD") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16)) 
generate_figs(aggplot_PD1_Tex_EMR, paste('./plots/', experiment, '_query_agg_aggplot_PD1_Tex_EMR', sep = ''), c(3,3))








