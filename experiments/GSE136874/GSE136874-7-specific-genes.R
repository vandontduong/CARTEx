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
experiment = 'GSE136874'
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

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('CAR', 'pblabels'), return.seurat = TRUE)
expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg$CAR <- factor(expt.obj.agg$CAR, levels = c("CD19", "GD2"))

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
# md$CXCL10 <- expt.obj.agg@assays$RNA$data['CXCL10',]
# md$XKRX <- expt.obj.agg@assays$RNA$data['XKRX',]
# md$ATP9A <- expt.obj.agg@assays$RNA$data['ATP9A',]
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


# PDCD1, HAVCR2, LAG3, CTLA4, NT5E, ENTPD1, METRNL, SNX9, CXCL10, XKRX, ATP9A, TRIB1, RAB32, EBI3, ING3, ITPRIPL2, MTURN, CTTN, KRT3, CYP1B1, IL2RA, DUSP5, MT1E, RGS16, AHI1, EGR2, LMCD1, RNF19B



# Plot for PDCD1
aggplot_PDCD1 <- md %>% ggplot(aes(CAR, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("PDCD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_explore_agg_aggplot_PDCD1', sep = ''), c(2, 2))

# Plot for HAVCR2
aggplot_HAVCR2 <- md %>% ggplot(aes(CAR, HAVCR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("HAVCR2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_HAVCR2, paste('./plots/', experiment, '_explore_agg_aggplot_HAVCR2', sep = ''), c(2, 2))

# Plot for LAG3
aggplot_LAG3 <- md %>% ggplot(aes(CAR, LAG3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("LAG3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_LAG3, paste('./plots/', experiment, '_explore_agg_aggplot_LAG3', sep = ''), c(2, 2))

# Plot for CTLA4
aggplot_CTLA4 <- md %>% ggplot(aes(CAR, CTLA4)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("CTLA4") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_CTLA4, paste('./plots/', experiment, '_explore_agg_aggplot_CTLA4', sep = ''), c(2, 2))

# Plot for NT5E
aggplot_NT5E <- md %>% ggplot(aes(CAR, NT5E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("NT5E") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_NT5E, paste('./plots/', experiment, '_explore_agg_aggplot_NT5E', sep = ''), c(2, 2))

# Plot for ENTPD1
aggplot_ENTPD1 <- md %>% ggplot(aes(CAR, ENTPD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("ENTPD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_ENTPD1, paste('./plots/', experiment, '_explore_agg_aggplot_ENTPD1', sep = ''), c(2, 2))

# Plot for METRNL
aggplot_METRNL <- md %>% ggplot(aes(CAR, METRNL)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("METRNL") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_METRNL, paste('./plots/', experiment, '_explore_agg_aggplot_METRNL', sep = ''), c(2, 2))

# Plot for SNX9
aggplot_SNX9 <- md %>% ggplot(aes(CAR, SNX9)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("SNX9") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_SNX9, paste('./plots/', experiment, '_explore_agg_aggplot_SNX9', sep = ''), c(2, 2))

# Plot for TRIB1
aggplot_TRIB1 <- md %>% ggplot(aes(CAR, TRIB1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("TRIB1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_TRIB1, paste('./plots/', experiment, '_explore_agg_aggplot_TRIB1', sep = ''), c(2, 2))

# Plot for RAB32
aggplot_RAB32 <- md %>% ggplot(aes(CAR, RAB32)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("RAB32") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_RAB32, paste('./plots/', experiment, '_explore_agg_aggplot_RAB32', sep = ''), c(2, 2))

# Plot for EBI3
aggplot_EBI3 <- md %>% ggplot(aes(CAR, EBI3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("EBI3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_EBI3, paste('./plots/', experiment, '_explore_agg_aggplot_EBI3', sep = ''), c(2, 2))

# Plot for ING3
aggplot_ING3 <- md %>% ggplot(aes(CAR, ING3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("ING3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_ING3, paste('./plots/', experiment, '_explore_agg_aggplot_ING3', sep = ''), c(2, 2))

# Plot for ITPRIPL2
aggplot_ITPRIPL2 <- md %>% ggplot(aes(CAR, ITPRIPL2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("ITPRIPL2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_ITPRIPL2, paste('./plots/', experiment, '_explore_agg_aggplot_ITPRIPL2', sep = ''), c(2, 2))

# Plot for MTURN
aggplot_MTURN <- md %>% ggplot(aes(CAR, MTURN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("MTURN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_MTURN, paste('./plots/', experiment, '_explore_agg_aggplot_MTURN', sep = ''), c(2, 2))

# Plot for CTTN
aggplot_CTTN <- md %>% ggplot(aes(CAR, CTTN)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("CTTN") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_CTTN, paste('./plots/', experiment, '_explore_agg_aggplot_CTTN', sep = ''), c(2, 2))

# Plot for KRT3
aggplot_KRT3 <- md %>% ggplot(aes(CAR, KRT3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("KRT3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_KRT3, paste('./plots/', experiment, '_explore_agg_aggplot_KRT3', sep = ''), c(2, 2))

# Plot for CYP1B1
aggplot_CYP1B1 <- md %>% ggplot(aes(CAR, CYP1B1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("CYP1B1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_CYP1B1, paste('./plots/', experiment, '_explore_agg_aggplot_CYP1B1', sep = ''), c(2, 2))

# Plot for IL2RA
aggplot_IL2RA <- md %>% ggplot(aes(CAR, IL2RA)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("IL2RA") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_IL2RA, paste('./plots/', experiment, '_explore_agg_aggplot_IL2RA', sep = ''), c(2, 2))

# Plot for DUSP5
aggplot_DUSP5 <- md %>% ggplot(aes(CAR, DUSP5)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("DUSP5") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_DUSP5, paste('./plots/', experiment, '_explore_agg_aggplot_DUSP5', sep = ''), c(2, 2))

# Plot for MT1E
aggplot_MT1E <- md %>% ggplot(aes(CAR, MT1E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("MT1E") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_MT1E, paste('./plots/', experiment, '_explore_agg_aggplot_MT1E', sep = ''), c(2, 2))

# Plot for RGS16
aggplot_RGS16 <- md %>% ggplot(aes(CAR, RGS16)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("RGS16") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_RGS16, paste('./plots/', experiment, '_explore_agg_aggplot_RGS16', sep = ''), c(2, 2))

# Plot for AHI1
aggplot_AHI1 <- md %>% ggplot(aes(CAR, AHI1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("AHI1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_AHI1, paste('./plots/', experiment, '_explore_agg_aggplot_AHI1', sep = ''), c(2, 2))

# Plot for EGR2
aggplot_EGR2 <- md %>% ggplot(aes(CAR, EGR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("EGR2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_EGR2, paste('./plots/', experiment, '_explore_agg_aggplot_EGR2', sep = ''), c(2, 2))

# Plot for LMCD1
aggplot_LMCD1 <- md %>% ggplot(aes(CAR, LMCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("LMCD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_LMCD1, paste('./plots/', experiment, '_explore_agg_aggplot_LMCD1', sep = ''), c(2, 2))

# Plot for RNF19B
aggplot_RNF19B <- md %>% ggplot(aes(CAR, RNF19B)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = CAR), color = "black") +
  scale_fill_manual(values = c('dodgerblue', 'indianred')) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("RNF19B") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("CD19", "GD2"))
generate_figs(aggplot_RNF19B, paste('./plots/', experiment, '_explore_agg_aggplot_RNF19B', sep = ''), c(2, 2))




aggplot_all_canonicals <- (aggplot_PDCD1 | aggplot_HAVCR2 | aggplot_LAG3 | aggplot_CTLA4 | aggplot_NT5E | aggplot_ENTPD1) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_all_canonicals, paste('./plots/', experiment, '_explore_agg_aggplot_all_canonicals', sep = ''), c(6,4))




aggplot_other_genes <- (aggplot_METRNL | aggplot_SNX9 | aggplot_TRIB1) +
  plot_layout(ncol = 3, nrow = 1, guides = "collect")

generate_figs(aggplot_other_genes, paste('./plots/', experiment, '_explore_agg_aggplot_other_genes', sep = ''), c(6,2))




aggplot_CRISPR_screens <- (aggplot_METRNL | aggplot_IL2RA | aggplot_RAB32 | aggplot_EBI3 | aggplot_ING3 | aggplot_ITPRIPL2 | aggplot_MTURN | aggplot_CTTN |
                             aggplot_KRT3 | aggplot_CYP1B1 | aggplot_DUSP5 | aggplot_MT1E | aggplot_RGS16 | aggplot_AHI1 | aggplot_EGR2 | aggplot_LMCD1 | aggplot_RNF19B) +
  plot_layout(ncol = 6, nrow = 3, guides = "collect")

generate_figs(aggplot_CRISPR_screens, paste('./plots/', experiment, '_explore_agg_aggplot_CRISPR_screens', sep = ''), c(12,6))











