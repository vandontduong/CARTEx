#!/usr/bin/env Rscript

### Script name: GSE136184-4-aggregate.R
### Description: aggregate seurat object for GSE136184 (extract)
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
##################################### Load single-cell dataset #####################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_cs_scored.rds', sep = ''))

# examine differentiation

vlnplot_CARTEx_age_group_monaco_split <- VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AgeGroup2', split.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + ylab("CARTEx 200") + ylim(-2,4)

generate_figs(vlnplot_CARTEx_age_group_monaco_split, paste('./plots/', experiment, '_cs_prepare_vlnplot_CARTEx_age_group_monaco_split', sep = ''), c(12,5)) 

custom_labels <- c("Naive", "Central memory", "Effector memory", "Terminal effector")
vlnplot_CARTEx_age_group_monaco <- plot_grid(VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'AgeGroup2', pt.size = 0, cols = colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$AgeGroup2)))) +
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4),
                                             VlnPlot(expt.obj, features = 'CARTEx_200', group.by = 'monaco', pt.size = 0, cols = c('deepskyblue','seagreen','darkgoldenrod','plum3')) + 
                                               theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                                               ylab("CARTEx 200") + ylim(-2,4) + scale_x_discrete(labels = custom_labels))

generate_figs(vlnplot_CARTEx_age_group_monaco, paste('./plots/', experiment, '_cs_prepare_vlnplot_CARTEx_age_group_monaco', sep = ''), c(12,5)) 







####################################################################################################
######################################## Pseudobulk the data #######################################
####################################################################################################

table(expt.obj@meta.data$AgeGroup2, expt.obj@meta.data$monaco)
# select single terminal effector CD8 T cell Newborn
id_newborn_teff <- Cells(expt.obj[, (expt.obj[['AgeGroup2']] == 'Newborn') & (expt.obj[['monaco']] == 'Terminal effector CD8 T cells')])
expt.obj <- subset(expt.obj, cells = setdiff(Cells(expt.obj), id_newborn_teff))

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 5)
md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup2", "Age")]
md[, .N, by = c("AgeGroup2", "Age", "Code")]
md[, .N, by = c("AgeGroup2", "Age", "Code", "pblabels")]


swarmplot_CARTEx_age_group_monaco <- ggplot(md, aes(x = AgeGroup2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + ylim(-2,4) +
  labs(x = "Age Group", y = "CARTEx 200", color = "Cell Type") +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position="none")
generate_figs(swarmplot_CARTEx_age_group_monaco, paste('./plots/', experiment, '_cs_prepare_swarmplot_CARTEx_age_group_monaco', sep = ''), c(6,5)) 


expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('AgeGroup2', 'monaco', 'pblabels'), return.seurat = TRUE)

expt.obj.agg.v2 <- AggregateExpression(expt.obj, group.by = c('AgeGroup2', 'pblabels'), return.seurat = TRUE)




####################################################################################################
######################################### Signature scoring ########################################
####################################################################################################

expt.obj.agg <- ScoreSubroutine(expt.obj.agg)
expt.obj.agg.v2 <- ScoreSubroutine(expt.obj.agg.v2)


expt.obj.agg$AgeGroup2 <- factor(expt.obj.agg$AgeGroup2, levels = c("Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))

expt.obj.agg.v2$AgeGroup2 <- factor(expt.obj.agg.v2$AgeGroup2, levels = c("Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))

saveRDS(expt.obj.agg, file = paste('./data/', experiment, '_cs_agg_scored.rds', sep = ''))
saveRDS(expt.obj.agg.v2, file = paste('./data/', experiment, '_cs_agg_scored_v2.rds', sep = ''))

# expt.obj.agg <- readRDS(paste('./data/', experiment, '_cs_agg_scored.rds', sep = ''))
# expt.obj.agg.v2 <- readRDS(paste('./data/', experiment, '_cs_agg_scored_v2.rds', sep = ''))

head(expt.obj.agg)
head(expt.obj.agg.v2)







md <- expt.obj.agg@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup2")]

glimpse(md)


aggplot_CARTEx_200_age_group <- md %>% ggplot(aes(AgeGroup2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2)) +
  scale_fill_manual(values = colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$AgeGroup2)))) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-3, 3) + theme(legend.position="none")


aggplot_CARTEx_200_monaco <- md %>% ggplot(aes(monaco, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = monaco)) +
  scale_fill_manual(values = c('deepskyblue','seagreen','darkgoldenrod','plum3')) +
  ylab("CARTEx 200") + xlab(NULL) + geom_point() + ylim(-3, 3) + theme(legend.position="none")

aggplot_CARTEx_200_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = CARTEx_200, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_CARTEx_200_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_CARTEx_200_age_group_monaco_split', sep = ''), c(6,5)) 


aggplot_activation_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = Activation, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_activation_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_activation_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_anergy_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = Anergy, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_anergy_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_anergy_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_senescence_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = Senescence, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_senescence_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_senescence_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_stemness_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = Stemness, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_stemness_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_stemness_age_group_monaco_split', sep = ''), c(6,5)) 


aggplot_NKlike_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = NKlike_Tex, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_NKlike_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_NKlike_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_LCMV_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = LCMV_Tex, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_LCMV_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_LCMV_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_BBD_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = BBD_Tex, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_BBD_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_BBD_age_group_monaco_split', sep = ''), c(6,5)) 

aggplot_PD1_age_group_monaco_split <- md %>% ggplot(aes(x = AgeGroup2, y = PD1_Tex, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE) + ylim(-2,2) +
  scale_color_manual(values = c('deepskyblue', 'seagreen', 'darkgoldenrod', 'plum3')) +
  theme_bw() + theme(axis.title.x = element_blank())
generate_figs(aggplot_PD1_age_group_monaco_split, paste('./plots/', experiment, '_cs_aggplot_PD1_age_group_monaco_split', sep = ''), c(6,5)) 



####







md <- expt.obj.agg.v2@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup2")]

glimpse(md)


umap_age_group_2_cols <- colorRampPalette(c("lightblue","orange", "orangered","violet"))(length(unique(expt.obj@meta.data$AgeGroup2)))
aggplot_CARTEx_200 <- md %>% ggplot(aes(AgeGroup2, CARTEx_200)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = AgeGroup2), color = "black") + geom_hline(yintercept=0) +
  scale_fill_manual(labels=c("N", "U30", "U50", "U70", "E"),values=umap_age_group_2_cols) + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Newborn','Under 30')), label = "p.signif", label.y = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Newborn','Under 50')), label = "p.signif", label.y = 1) +
  # stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre-NR','OldTerminal')), label = "p.signif", label.y = 1.7) +
  ylab("CARTEx") + xlab(NULL) + geom_point() + ylim(-2, 2) + theme_classic() + theme(legend.position="none", text=element_text(size=16, color = "black")) +
  scale_x_discrete(labels = c("N", "U30", "U50", "U70", "E"))
generate_figs(aggplot_CARTEx_200, paste('./plots/', experiment, '_query_agg_aggplot_CARTEx_200', sep = ''), c(3,3))





