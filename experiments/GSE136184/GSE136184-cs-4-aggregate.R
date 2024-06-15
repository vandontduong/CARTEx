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






####################################################################################################
########################################### CARTEx scoring #########################################
####################################################################################################

# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_630 <- Z(SignatureScore(expt.obj.agg, cartex_630_weights))
expt.obj.agg@meta.data$CARTEx_630i <- integerize(expt.obj.agg@meta.data$CARTEx_630)

# CARTEx with weights // 200 genes
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_200 <- Z(SignatureScore(expt.obj.agg, cartex_200_weights))
expt.obj.agg@meta.data$CARTEx_200i <- integerize(expt.obj.agg@meta.data$CARTEx_200)

# CARTEx with weights // 84 genes
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
expt.obj.agg@meta.data$CARTEx_84 <- Z(SignatureScore(expt.obj.agg, cartex_84_weights))
expt.obj.agg@meta.data$CARTEx_84i <- integerize(expt.obj.agg@meta.data$CARTEx_84)

####################################################################################################
########################################### Module scoring #########################################
####################################################################################################

activation.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = TRUE, row.names = 1))
senescence.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = TRUE, row.names = 1))

expt.obj.agg <- AddModuleScore(expt.obj.agg, features = list(activation.sig, anergy.sig, stemness.sig, senescence.sig), name="State", search = TRUE)

# z score normalization
expt.obj.agg@meta.data$Activation <- scale(expt.obj.agg@meta.data$State1)
expt.obj.agg@meta.data$Anergy <- scale(expt.obj.agg@meta.data$State2)
expt.obj.agg@meta.data$Stemness <- scale(expt.obj.agg@meta.data$State3)
expt.obj.agg@meta.data$Senescence <- scale(expt.obj.agg@meta.data$State4)

expt.obj.agg@meta.data$Activationi <- integerize(expt.obj.agg@meta.data$Activation)
expt.obj.agg@meta.data$Anergyi <- integerize(expt.obj.agg@meta.data$Anergy)
expt.obj.agg@meta.data$Stemnessi <- integerize(expt.obj.agg@meta.data$Stemness)
expt.obj.agg@meta.data$Senescencei <- integerize(expt.obj.agg@meta.data$Senescence)

expt.obj.agg@meta.data$Activationi <- factor(expt.obj.agg@meta.data$Activationi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Anergyi <- factor(expt.obj.agg@meta.data$Anergyi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Stemnessi <- factor(expt.obj.agg@meta.data$Stemnessi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
expt.obj.agg@meta.data$Senescencei <- factor(expt.obj.agg@meta.data$Senescencei, levels = c(4,3,2,1,0,-1,-2,-3,-4))

expt.obj.agg@meta.data$State1 <- NULL
expt.obj.agg@meta.data$State2 <- NULL
expt.obj.agg@meta.data$State3 <- NULL
expt.obj.agg@meta.data$State4 <- NULL

# examine other signatures
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))

expt.obj.agg <- AddModuleScore(expt.obj.agg, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex), name="Signature", search = TRUE)

expt.obj.agg@meta.data$NKlike_Tex <- scale(expt.obj.agg@meta.data$Signature1)
expt.obj.agg@meta.data$LCMV_Tex <- scale(expt.obj.agg@meta.data$Signature2)
expt.obj.agg@meta.data$BBD_Tex <- scale(expt.obj.agg@meta.data$Signature3)
expt.obj.agg@meta.data$PD1_Tex <- scale(expt.obj.agg@meta.data$Signature4)

expt.obj.agg@meta.data$NKlike_Texi <- integerize(expt.obj.agg@meta.data$NKlike_Tex)
expt.obj.agg@meta.data$LCMV_Texi <- integerize(expt.obj.agg@meta.data$LCMV_Tex)
expt.obj.agg@meta.data$BBD_Texi <- integerize(expt.obj.agg@meta.data$BBD_Tex)
expt.obj.agg@meta.data$PD1_Texi <- integerize(expt.obj.agg@meta.data$PD1_Tex)

expt.obj.agg@meta.data$Signature1 <- NULL
expt.obj.agg@meta.data$Signature2 <- NULL
expt.obj.agg@meta.data$Signature3 <- NULL
expt.obj.agg@meta.data$Signature4 <- NULL


expt.obj.agg$AgeGroup2 <- factor(expt.obj.agg$AgeGroup2, levels = c("Newborn", "Under 30", "Under 50", "Under 70", "Elderly"))
expt.obj.agg$monaco <- factor(expt.obj.agg$monaco, levels = c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells"))


saveRDS(expt.obj.agg, file = paste('./data/', experiment, '_cs_agg_scored.rds', sep = ''))

# expt.obj.agg <- readRDS(paste('./data/', experiment, '_cs_agg_scored.rds', sep = ''))

head(expt.obj.agg)



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






