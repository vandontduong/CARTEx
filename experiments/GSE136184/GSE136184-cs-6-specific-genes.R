#!/usr/bin/env Rscript

### Script name: GSE136184-cs-6-specific-genes.R
### Description: explore seurat object for GSE136184
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136184'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))


####################################################################################################
#################################### Examine canonical exhaustion ##################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_cs_scored.rds', sep = ''))

# PDCD1 corresponds to PD-1
# HAVCR2 corresponds to TIM-3
# LAG3
# CTLA4
# TIGIT
# ENTPD1 corresponds to CD39





# aggregate for CARTEx scores (BAR CHART)


expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AgeGroup2', 'pblabels'), return.seurat = TRUE)
expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg$AgeGroup2 <- factor(expt.obj.agg$AgeGroup2, levels = c("Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))

md <- expt.obj.agg@meta.data %>% as.data.table
glimpse(md)



# examine canonical markers
md$PDCD1 <- expt.obj.agg@assays$RNA$data['PDCD1',] # PD-1
md$HAVCR2 <- expt.obj.agg@assays$RNA$data['HAVCR2',] # TIM3
md$LAG3 <- expt.obj.agg@assays$RNA$data['LAG3',]
md$CTLA4 <- expt.obj.agg@assays$RNA$data['CTLA4',]
md$TIGIT <- expt.obj.agg@assays$RNA$data['TIGIT',]
md$NT5E <- expt.obj.agg@assays$RNA$data['NT5E',] # CD73
# md$ENTPD1 <- expt.obj.agg@assays$RNA$data['ENTPD1',] # CD39
# md$TNFSRF9 <- expt.obj.agg@assays$RNA$data['TNFSRF9',]


# other genes
md$METRNL <- expt.obj.agg@assays$RNA$data['METRNL',]
md$SNX9 <- expt.obj.agg@assays$RNA$data['SNX9',]
# md$CXCL10 <- expt.obj.agg@assays$RNA$data['CXCL10',]
# md$XKRX <- expt.obj.agg@assays$RNA$data['XKRX',]
# md$ATP9A <- expt.obj.agg@assays$RNA$data['ATP9A',]
# md$TRIB1 <- expt.obj.agg@assays$RNA$data['TRIB1',]




umap_age_group_2_cols <- colorRampPalette(c("lightblue","orange", "orangered","violet"))(length(unique(expt.obj@meta.data$AgeGroup2)))


# Plot for PDCD1
aggplot_PDCD1 <- md %>% ggplot(aes(AgeGroup2, PDCD1)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("PDCD1") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_PDCD1, paste('./plots/', experiment, '_explore_agg_aggplot_PDCD1', sep = ''), c(2, 2))

# Plot for HAVCR2
aggplot_HAVCR2 <- md %>% ggplot(aes(AgeGroup2, HAVCR2)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("HAVCR2") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_HAVCR2, paste('./plots/', experiment, '_explore_agg_aggplot_HAVCR2', sep = ''), c(2, 2))

# Plot for LAG3
aggplot_LAG3 <- md %>% ggplot(aes(AgeGroup2, LAG3)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("LAG3") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_LAG3, paste('./plots/', experiment, '_explore_agg_aggplot_LAG3', sep = ''), c(2, 2))

# Plot for CTLA4
aggplot_CTLA4 <- md %>% ggplot(aes(AgeGroup2, CTLA4)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("CTLA4") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_CTLA4, paste('./plots/', experiment, '_explore_agg_aggplot_CTLA4', sep = ''), c(2, 2))

# Plot for TIGIT
aggplot_TIGIT <- md %>% ggplot(aes(AgeGroup2, TIGIT)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("TIGIT") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_TIGIT, paste('./plots/', experiment, '_explore_agg_aggplot_TIGIT', sep = ''), c(2, 2))


# Plot for NT5E
aggplot_NT5E <- md %>% ggplot(aes(AgeGroup2, NT5E)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") +
  scale_fill_manual(values = umap_age_group_2_cols) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c("CD19", "GD2")), label = "p.signif", label.y = 1.8) +
  ylab("NT5E") + xlab(NULL) + geom_point() + ylim(0, 2) + theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_NT5E, paste('./plots/', experiment, '_explore_agg_aggplot_NT5E', sep = ''), c(2, 2))




aggplot_all_canonicals <- (aggplot_PDCD1 | aggplot_HAVCR2 | aggplot_LAG3 | aggplot_CTLA4 | aggplot_TIGIT | aggplot_NT5E) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect")

generate_figs(aggplot_all_canonicals, paste('./plots/', experiment, '_explore_agg_aggplot_all_canonicals', sep = ''), c(9,4))


