#!/usr/bin/env Rscript

### Script name: GSE160160-7-specific-genes.R
### Description: explore seurat object for GSE160160
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39











# aggregate for CARTEx scores (BAR CHART)


expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('exposure', 'pblabels'), return.seurat = TRUE)
expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg$exposure <- factor(expt.obj.agg$exposure, levels = c("day0", "day20"))

md <- expt.obj.agg@meta.data %>% as.data.table
glimpse(md)



# examine canonical markers
md$PDCD1 <- expt.obj.agg@assays$RNA$data['PDCD1',] # PD-1
md$HAVCR2 <- expt.obj.agg@assays$RNA$data['HAVCR2',] # TIM3
md$LAG3 <- expt.obj.agg@assays$RNA$data['LAG3',]
md$CTLA4 <- expt.obj.agg@assays$RNA$data['CTLA4',]
md$NT5E <- expt.obj.agg@assays$RNA$data['NT5E',] # CD73
md$ENTPD1 <- expt.obj.agg@assays$RNA$data['ENTPD1',] # CD39

# other genes
md$METRNL <- expt.obj.agg@assays$RNA$data['METRNL',]
md$SNX9 <- expt.obj.agg@assays$RNA$data['SNX9',]
md$CXCL10 <- expt.obj.agg@assays$RNA$data['CXCL10',]
md$XKRX <- expt.obj.agg@assays$RNA$data['XKRX',]
md$ATP9A <- expt.obj.agg@assays$RNA$data['ATP9A',]
md$TRIB1 <- expt.obj.agg@assays$RNA$data['TRIB1',]

# CRISPR screen genes (METRNL and CXCL10 already captured)
md$RAB32 <- expt.obj.agg@assays$RNA$data['RAB32',]
md$EBI3 <- expt.obj.agg@assays$RNA$data['RAB32',]
md$ING3 <- expt.obj.agg@assays$RNA$data['ING3',]
md$ITPRIPL2 <- expt.obj.agg@assays$RNA$data['ITPRIPL2',]
md$MTURN <- expt.obj.agg@assays$RNA$data['MTURN',]
md$CTTN <- expt.obj.agg@assays$RNA$data['CTTN',]
md$KRT3 <- expt.obj.agg@assays$RNA$data['KRT3',]
md$CYP1B1 <- expt.obj.agg@assays$RNA$data['CYP1B1',]
md$IL2RA <- expt.obj.agg@assays$RNA$data['IL2RA',]
md$DUSP5 <- expt.obj.agg@assays$RNA$data['DUSP5',]
md$MT1E <- expt.obj.agg@assays$RNA$data['MT1E',]
md$RGS16 <- expt.obj.agg@assays$RNA$data['RGS16',]
md$AHI1 <- expt.obj.agg@assays$RNA$data['AHI1',]
md$EGR2 <- expt.obj.agg@assays$RNA$data['EGR2',]
md$LMCD1 <- expt.obj.agg@assays$RNA$data['LMCD1',]
md$RNF19B <- expt.obj.agg@assays$RNA$data['RNF19B',]




aggplot_PDCD1 <- md %>% ggplot(aes(exposure, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("PDCD1") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_explore_agg_aggplot_PDCD1', sep = ''), c(2,2))

aggplot_HAVCR2 <- md %>% ggplot(aes(exposure, HAVCR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("HAVCR2") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_HAVCR2, paste('./plots/', experiment, '_explore_agg_aggplot_HAVCR2', sep = ''), c(2,2))

aggplot_LAG3 <- md %>% ggplot(aes(exposure, LAG3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("LAG3") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_LAG3, paste('./plots/', experiment, '_explore_agg_aggplot_LAG3', sep = ''), c(2,2))

aggplot_CTLA4 <- md %>% ggplot(aes(exposure, CTLA4)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("CTLA4") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_CTLA4, paste('./plots/', experiment, '_explore_agg_aggplot_CTLA4', sep = ''), c(2,2))

aggplot_NT5E <- md %>% ggplot(aes(exposure, NT5E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("NT5E") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_NT5E, paste('./plots/', experiment, '_explore_agg_aggplot_NT5E', sep = ''), c(2,2))

aggplot_ENTPD1 <- md %>% ggplot(aes(exposure, ENTPD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("ENTPD1") + xlab(NULL) + geom_point() + ylim(-0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_ENTPD1, paste('./plots/', experiment, '_explore_agg_aggplot_ENTPD1', sep = ''), c(2,2))




aggplot_all_canonicals <- (aggplot_PDCD1 | aggplot_HAVCR2 | aggplot_LAG3 | aggplot_CTLA4 | aggplot_NT5E | aggplot_ENTPD1) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_all_canonicals, paste('./plots/', experiment, '_explore_agg_aggplot_all_canonicals', sep = ''), c(6,4))






aggplot_METRNL <- md %>% ggplot(aes(exposure, METRNL)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("METRNL") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_METRNL, paste('./plots/', experiment, '_explore_agg_aggplot_METRNL', sep = ''), c(2, 2))

aggplot_SNX9 <- md %>% ggplot(aes(exposure, SNX9)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("SNX9") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_SNX9, paste('./plots/', experiment, '_explore_agg_aggplot_SNX9', sep = ''), c(2, 2))

aggplot_CXCL10 <- md %>% ggplot(aes(exposure, CXCL10)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("CXCL10") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_CXCL10, paste('./plots/', experiment, '_explore_agg_aggplot_CXCL10', sep = ''), c(2, 2))

aggplot_XKRX <- md %>% ggplot(aes(exposure, XKRX)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("XKRX") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_XKRX, paste('./plots/', experiment, '_explore_agg_aggplot_XKRX', sep = ''), c(2, 2))

aggplot_ATP9A <- md %>% ggplot(aes(exposure, ATP9A)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("ATP9A") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_ATP9A, paste('./plots/', experiment, '_explore_agg_aggplot_ATP9A', sep = ''), c(2, 2))

aggplot_TRIB1 <- md %>% ggplot(aes(exposure, TRIB1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("TRIB1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_TRIB1, paste('./plots/', experiment, '_explore_agg_aggplot_TRIB1', sep = ''), c(2, 2))

aggplot_RAB32 <- md %>% ggplot(aes(exposure, RAB32)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("RAB32") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_RAB32, paste('./plots/', experiment, '_explore_agg_aggplot_RAB32', sep = ''), c(2, 2))

aggplot_EBI3 <- md %>% ggplot(aes(exposure, EBI3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("EBI3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_EBI3, paste('./plots/', experiment, '_explore_agg_aggplot_EBI3', sep = ''), c(2, 2))

aggplot_ING3 <- md %>% ggplot(aes(exposure, ING3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("ING3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_ING3, paste('./plots/', experiment, '_explore_agg_aggplot_ING3', sep = ''), c(2, 2))

aggplot_ITPRIPL2 <- md %>% ggplot(aes(exposure, ITPRIPL2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("ITPRIPL2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_ITPRIPL2, paste('./plots/', experiment, '_explore_agg_aggplot_ITPRIPL2', sep = ''), c(2, 2))

aggplot_MTURN <- md %>% ggplot(aes(exposure, MTURN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("MTURN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_MTURN, paste('./plots/', experiment, '_explore_agg_aggplot_MTURN', sep = ''), c(2, 2))

aggplot_CTTN <- md %>% ggplot(aes(exposure, CTTN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("CTTN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_CTTN, paste('./plots/', experiment, '_explore_agg_aggplot_CTTN', sep = ''), c(2, 2))

aggplot_KRT3 <- md %>% ggplot(aes(exposure, KRT3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("KRT3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_KRT3, paste('./plots/', experiment, '_explore_agg_aggplot_KRT3', sep = ''), c(2, 2))

aggplot_CYP1B1 <- md %>% ggplot(aes(exposure, CYP1B1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("CYP1B1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_CYP1B1, paste('./plots/', experiment, '_explore_agg_aggplot_CYP1B1', sep = ''), c(2, 2))

aggplot_IL2RA <- md %>% ggplot(aes(exposure, IL2RA)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("IL2RA") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_IL2RA, paste('./plots/', experiment, '_explore_agg_aggplot_IL2RA', sep = ''), c(2, 2))

aggplot_DUSP5 <- md %>% ggplot(aes(exposure, DUSP5)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("DUSP5") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_DUSP5, paste('./plots/', experiment, '_explore_agg_aggplot_DUSP5', sep = ''), c(2, 2))

aggplot_MT1E <- md %>% ggplot(aes(exposure, MT1E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("MT1E") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_MT1E, paste('./plots/', experiment, '_explore_agg_aggplot_MT1E', sep = ''), c(2, 2))

aggplot_RGS16 <- md %>% ggplot(aes(exposure, RGS16)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("RGS16") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_RGS16, paste('./plots/', experiment, '_explore_agg_aggplot_RGS16', sep = ''), c(2, 2))

aggplot_AHI1 <- md %>% ggplot(aes(exposure, AHI1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("AHI1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_AHI1, paste('./plots/', experiment, '_explore_agg_aggplot_AHI1', sep = ''), c(2, 2))

aggplot_EGR2 <- md %>% ggplot(aes(exposure, EGR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("EGR2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_EGR2, paste('./plots/', experiment, '_explore_agg_aggplot_EGR2', sep = ''), c(2, 2))

aggplot_LMCD1 <- md %>% ggplot(aes(exposure, LMCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("LMCD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_LMCD1, paste('./plots/', experiment, '_explore_agg_aggplot_LMCD1', sep = ''), c(2, 2))

aggplot_RNF19B <- md %>% ggplot(aes(exposure, RNF19B)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure), color = "black") +
  scale_fill_manual(values = c("skyblue", "cadetblue")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("RNF19B") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("D0", "D20"))
generate_figs(aggplot_RNF19B, paste('./plots/', experiment, '_explore_agg_aggplot_RNF19B', sep = ''), c(2, 2))




aggplot_other_genes <- (aggplot_METRNL | aggplot_SNX9 | aggplot_CXCL10 | aggplot_XKRX | aggplot_ATP9A | aggplot_TRIB1) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_other_genes, paste('./plots/', experiment, '_explore_agg_aggplot_other_genes', sep = ''), c(6,4))




aggplot_CRISPR_screens <- (aggplot_METRNL | aggplot_CXCL10 | aggplot_IL2RA | aggplot_RAB32 | aggplot_EBI3 | aggplot_ING3 | aggplot_ITPRIPL2 | aggplot_MTURN | aggplot_CTTN |
                             aggplot_KRT3 | aggplot_CYP1B1 | aggplot_DUSP5 | aggplot_MT1E | aggplot_RGS16 | aggplot_AHI1 | aggplot_EGR2 | aggplot_LMCD1 | aggplot_RNF19B) +
  plot_layout(ncol = 6, nrow = 3, guides = "collect")

generate_figs(aggplot_CRISPR_screens, paste('./plots/', experiment, '_explore_agg_aggplot_CRISPR_screens', sep = ''), c(12,6))













# aggregate by subtype as well

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('exposure', 'monaco', 'pblabels'), return.seurat = TRUE)


expt.obj.agg <- ScoreSubroutine(expt.obj.agg)

expt.obj.agg$exposure <- factor(expt.obj.agg$exposure, levels = c("day0", "day20"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


md <- expt.obj.agg@meta.data %>% as.data.table



md[, subtype := fcase(
  monaco == "Naive CD8 T cells", "naive",
  monaco == "Central memory CD8 T cells", "central_memory",
  monaco == "Effector memory CD8 T cells", "effector_memory",
  monaco == "Terminal effector CD8 T cells", "terminal_effector",
  default = NA_character_  # Optional: handle unexpected values
)]

md[, exposure_subtype := paste(exposure, subtype, sep = "_")]
md[, exposure_subtype := factor(exposure_subtype, levels = c("day0_naive", "day0_central_memory", "day0_effector_memory", "day0_terminal_effector",
                                                             "day20_naive", "day20_central_memory", "day20_effector_memory", "day20_terminal_effector"))]



head(md)


md$PDCD1 <- expt.obj.agg@assays$RNA$data['PDCD1',] # PD-1
md$HAVCR2 <- expt.obj.agg@assays$RNA$data['HAVCR2',] # TIM3
md$LAG3 <- expt.obj.agg@assays$RNA$data['LAG3',]
md$CTLA4 <- expt.obj.agg@assays$RNA$data['CTLA4',]
md$NT5E <- expt.obj.agg@assays$RNA$data['NT5E',] # CD73
md$ENTPD1 <- expt.obj.agg@assays$RNA$data['ENTPD1',] # CD39


# other genes
md$METRNL <- expt.obj.agg@assays$RNA$data['METRNL',]
md$SNX9 <- expt.obj.agg@assays$RNA$data['SNX9',]
md$CXCL10 <- expt.obj.agg@assays$RNA$data['CXCL10',]
md$XKRX <- expt.obj.agg@assays$RNA$data['XKRX',]
md$ATP9A <- expt.obj.agg@assays$RNA$data['ATP9A',]
md$TRIB1 <- expt.obj.agg@assays$RNA$data['TRIB1',]

# CRISPR screen genes (METRNL and CXCL10 already captured)
md$RAB32 <- expt.obj.agg@assays$RNA$data['RAB32',]
md$EBI3 <- expt.obj.agg@assays$RNA$data['RAB32',]
md$ING3 <- expt.obj.agg@assays$RNA$data['ING3',]
md$ITPRIPL2 <- expt.obj.agg@assays$RNA$data['ITPRIPL2',]
md$MTURN <- expt.obj.agg@assays$RNA$data['MTURN',]
md$CTTN <- expt.obj.agg@assays$RNA$data['CTTN',]
md$KRT3 <- expt.obj.agg@assays$RNA$data['KRT3',]
md$CYP1B1 <- expt.obj.agg@assays$RNA$data['CYP1B1',]
md$IL2RA <- expt.obj.agg@assays$RNA$data['IL2RA',]
md$DUSP5 <- expt.obj.agg@assays$RNA$data['DUSP5',]
md$MT1E <- expt.obj.agg@assays$RNA$data['MT1E',]
md$RGS16 <- expt.obj.agg@assays$RNA$data['RGS16',]
md$AHI1 <- expt.obj.agg@assays$RNA$data['AHI1',]
md$EGR2 <- expt.obj.agg@assays$RNA$data['EGR2',]
md$LMCD1 <- expt.obj.agg@assays$RNA$data['LMCD1',]
md$RNF19B <- expt.obj.agg@assays$RNA$data['RNF19B',]


aggplot_subtype_PDCD1 <- md %>% ggplot(aes(exposure_subtype, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("PDCD1") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_PDCD1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_PDCD1', sep = ''), c(2,2))

aggplot_subtype_HAVCR2 <- md %>% ggplot(aes(exposure_subtype, HAVCR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("HAVCR2") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_HAVCR2, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_HAVCR2', sep = ''), c(2,2))

aggplot_subtype_LAG3 <- md %>% ggplot(aes(exposure_subtype, LAG3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("LAG3") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_LAG3, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_LAG3', sep = ''), c(2,2))

aggplot_subtype_CTLA4 <- md %>% ggplot(aes(exposure_subtype, CTLA4)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("CTLA4") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_CTLA4, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_CTLA4', sep = ''), c(2,2))

aggplot_subtype_NT5E <- md %>% ggplot(aes(exposure_subtype, NT5E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("NT5E") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_NT5E, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_NT5E', sep = ''), c(2,2))

aggplot_subtype_ENTPD1 <- md %>% ggplot(aes(exposure_subtype, ENTPD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") + # geom_hline(yintercept=0) +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('day0','day20')), label = "p.signif", label.y = 1.8) +
  ylab("ENTPD1") + xlab(NULL) + geom_point() + ylim(0,2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) + annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_ENTPD1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_ENTPD1', sep = ''), c(2,2))


aggplot_subtype_all_canonicals <- (aggplot_subtype_PDCD1 | aggplot_subtype_HAVCR2 | aggplot_subtype_LAG3 | aggplot_subtype_CTLA4 | aggplot_subtype_NT5E | aggplot_subtype_ENTPD1) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_subtype_all_canonicals, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_all_canonicals', sep = ''), c(10,4))



aggplot_subtype_METRNL <- md %>% ggplot(aes(exposure_subtype, METRNL)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("METRNL") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_METRNL, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_METRNL', sep = ''), c(2, 2))

aggplot_subtype_SNX9 <- md %>% ggplot(aes(exposure_subtype, SNX9)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("SNX9") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_SNX9, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_SNX9', sep = ''), c(2, 2))

aggplot_subtype_CXCL10 <- md %>% ggplot(aes(exposure_subtype, CXCL10)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("CXCL10") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_CXCL10, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_CXCL10', sep = ''), c(2, 2))

aggplot_subtype_XKRX <- md %>% ggplot(aes(exposure_subtype, XKRX)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("XKRX") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_XKRX, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_XKRX', sep = ''), c(2, 2))

aggplot_subtype_ATP9A <- md %>% ggplot(aes(exposure_subtype, ATP9A)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("ATP9A") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_ATP9A, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_ATP9A', sep = ''), c(2, 2))

aggplot_subtype_TRIB1 <- md %>% ggplot(aes(exposure_subtype, TRIB1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("TRIB1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_TRIB1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_TRIB1', sep = ''), c(2, 2))

aggplot_subtype_RAB32 <- md %>% ggplot(aes(exposure_subtype, RAB32)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("RAB32") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_RAB32, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_RAB32', sep = ''), c(2, 2))

aggplot_subtype_EBI3 <- md %>% ggplot(aes(exposure_subtype, EBI3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("EBI3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_EBI3, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_EBI3', sep = ''), c(2, 2))

aggplot_subtype_ING3 <- md %>% ggplot(aes(exposure_subtype, ING3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("ING3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_ING3, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_ING3', sep = ''), c(2, 2))

aggplot_subtype_ITPRIPL2 <- md %>% ggplot(aes(exposure_subtype, ITPRIPL2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("ITPRIPL2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_ITPRIPL2, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_ITPRIPL2', sep = ''), c(2, 2))

aggplot_subtype_MTURN <- md %>% ggplot(aes(exposure_subtype, MTURN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("MTURN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_MTURN, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_MTURN', sep = ''), c(2, 2))

aggplot_subtype_CTTN <- md %>% ggplot(aes(exposure_subtype, CTTN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("CTTN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_CTTN, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_CTTN', sep = ''), c(2, 2))

aggplot_subtype_KRT3 <- md %>% ggplot(aes(exposure_subtype, KRT3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("KRT3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_KRT3, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_KRT3', sep = ''), c(2, 2))

aggplot_subtype_CYP1B1 <- md %>% ggplot(aes(exposure_subtype, CYP1B1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("CYP1B1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_CYP1B1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_CYP1B1', sep = ''), c(2, 2))

aggplot_subtype_IL2RA <- md %>% ggplot(aes(exposure_subtype, IL2RA)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("IL2RA") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_IL2RA, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_IL2RA', sep = ''), c(2, 2))

aggplot_subtype_DUSP5 <- md %>% ggplot(aes(exposure_subtype, DUSP5)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("DUSP5") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_DUSP5, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_DUSP5', sep = ''), c(2, 2))

aggplot_subtype_MT1E <- md %>% ggplot(aes(exposure_subtype, MT1E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("MT1E") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_MT1E, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_MT1E', sep = ''), c(2, 2))

aggplot_subtype_RGS16 <- md %>% ggplot(aes(exposure_subtype, RGS16)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("RGS16") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_RGS16, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_RGS16', sep = ''), c(2, 2))

aggplot_subtype_AHI1 <- md %>% ggplot(aes(exposure_subtype, AHI1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("AHI1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_AHI1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_AHI1', sep = ''), c(2, 2))

aggplot_subtype_EGR2 <- md %>% ggplot(aes(exposure_subtype, EGR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("EGR2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_EGR2, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_EGR2', sep = ''), c(2, 2))

aggplot_subtype_LMCD1 <- md %>% ggplot(aes(exposure_subtype, LMCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("LMCD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_LMCD1, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_LMCD1', sep = ''), c(2, 2))

aggplot_subtype_RNF19B <- md %>% ggplot(aes(exposure_subtype, RNF19B)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = exposure_subtype), color = "black") +
  scale_fill_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3', 'deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  ylab("RNF19B") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "CM", "EM", "TE", "N", "CM", "EM", "TE")) +
  annotate("text", x = c(2.5, 6.5), y = 2, label = c("Day 0", "Day 20"), size = 4)
generate_figs(aggplot_subtype_RNF19B, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_RNF19B', sep = ''), c(2, 2))




aggplot_subtype_other_genes <- (aggplot_subtype_METRNL | aggplot_subtype_SNX9 | aggplot_subtype_CXCL10 | aggplot_subtype_XKRX | aggplot_subtype_ATP9A | aggplot_subtype_TRIB1) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_subtype_other_genes, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_other_genes', sep = ''), c(10,4))




aggplot_subtype_CRISPR_screens <- (aggplot_subtype_METRNL | aggplot_subtype_CXCL10 | aggplot_subtype_IL2RA | aggplot_subtype_RAB32 | aggplot_subtype_EBI3 | aggplot_subtype_ING3 | aggplot_subtype_ITPRIPL2 | aggplot_subtype_MTURN | aggplot_subtype_CTTN |
                             aggplot_subtype_KRT3 | aggplot_subtype_CYP1B1 | aggplot_subtype_DUSP5 | aggplot_subtype_MT1E | aggplot_subtype_RGS16 | aggplot_subtype_AHI1 | aggplot_subtype_EGR2 | aggplot_subtype_LMCD1 | aggplot_subtype_RNF19B) +
  plot_layout(ncol = 6, nrow = 3, guides = "collect")

generate_figs(aggplot_subtype_CRISPR_screens, paste('./plots/', experiment, '_explore_agg_aggplot_subtype_CRISPR_screens', sep = ''), c(20,6))












