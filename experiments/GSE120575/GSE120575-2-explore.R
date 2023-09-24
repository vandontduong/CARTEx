# explore seurat object
# Vandon Duong

experiment = 'GSE120575'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(gtools)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

CART.combined.CD8pos <- readRDS(paste('./data/', experiment, '_CD8pos_scored.rds', sep = ''))
CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "characteristics_response")

head(CART.combined.CD8pos)


####################################################################################################
########################################## UMAP generation #########################################
####################################################################################################

umap_seurat_clusters <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, seed = 123)
generate_figs(umap_seurat_clusters, paste('./plots/', experiment, '_umap_seurat_clusters', sep = ''))

umap_patient <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, seed = 123)
generate_figs(umap_patient, paste('./plots/', experiment, '_umap_patient', sep = ''))

umap_response <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "characteristics_response", shuffle = TRUE, seed = 123)
generate_figs(umap_response, paste('./plots/', experiment, '_umap_response', sep = ''))

umap_therapy <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "characteristics_therapy", shuffle = TRUE, seed = 123)
generate_figs(umap_therapy, paste('./plots/', experiment, '_umap_therapy', sep = ''))

umap_timepoint <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Timepoint", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint, paste('./plots/', experiment, '_umap_timepoint', sep = ''))

umap_timepoint_response <- DimPlot(CART.combined.CD8pos, reduction = "umap", group.by = "Timepoint_Response", shuffle = TRUE, seed = 123)
generate_figs(umap_timepoint_response, paste('./plots/', experiment, '_umap_timepoint_response', sep = ''))


# ANALYZE CELL CYCLE DATA FROM SEURAT OBJECT
md <- CART.combined.CD8pos@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "Phase")]
pie_phase <- md[, .N, by = c("Phase")]
setorder(pie_phase, cols = "Phase")
pie_phase$percent <- round(100*pie_phase$N / sum(pie_phase$N), digits = 1)
pie_phase$label <- paste(pie_phase$Phase, " (", pie_phase$percent,"%)", sep = "")
pie_phase$cols <- hcl.colors(n = nrow(pie_phase), palette = "Temps")
barplot_phase <- ggplot(data= pie_phase, aes(x=Phase, y=percent)) + geom_bar(stat="identity", fill = pie_phase$cols)
generate_figs(barplot_phase, paste('./plots/', experiment, '_barplot_phase', sep = ''))
umap_phase <- DimPlot(CART.combined.CD8pos, group.by = "Phase", cols = pie_phase$cols)
generate_figs(umap_phase, paste('./plots/', experiment, '_umap_phase', sep = ''))


# UMAP of scores
umap_CARTEx_630i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_630i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_630i, paste('./plots/', experiment, '_umap_CARTEx_630i', sep = ''))

umap_CARTEx_200i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_200i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_200i, paste('./plots/', experiment, '_umap_CARTEx_200i', sep = ''))

umap_CARTEx_84i <- DimPlot(CART.combined.CD8pos, group.by = "CARTEx_84i", shuffle = TRUE, seed = 123)
generate_figs(umap_CARTEx_84i, paste('./plots/', experiment, '_umap_CARTEx_84i', sep = ''))

umap_sig_activationi <- DimPlot(CART.combined.CD8pos, group.by = "Activationi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_activationi, paste('./plots/', experiment, '_umap_sig_activationi', sep = ''))

umap_sig_anergyi <- DimPlot(CART.combined.CD8pos, group.by = "Anergyi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_anergyi, paste('./plots/', experiment, '_umap_sig_anergyi', sep = ''))

umap_sig_stemnessi <- DimPlot(CART.combined.CD8pos, group.by = "Stemnessi", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_stemnessi, paste('./plots/', experiment, '_umap_sig_stemnessi', sep = ''))

umap_sig_senescencei <- DimPlot(CART.combined.CD8pos, group.by = "Senescencei", shuffle = TRUE, seed = 123)
generate_figs(umap_sig_senescencei, paste('./plots/', experiment, '_umap_sig_senescencei', sep = ''))


# better way to generate UMAPs for scored cells

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
umap_CARTEx_84 <- FeaturePlot(CART.combined.CD8pos, features = c("CARTEx_84"), order = TRUE) + fix.sc
umap_sig_activation <- FeaturePlot(CART.combined.CD8pos, features = c("Activation"), order = TRUE) + fix.sc
umap_sig_anergy <- FeaturePlot(CART.combined.CD8pos, features = c("Anergy"), order = TRUE) + fix.sc
umap_sig_stemness <- FeaturePlot(CART.combined.CD8pos, features = c("Stemness"), order = TRUE) + fix.sc
umap_sig_senescence <- FeaturePlot(CART.combined.CD8pos, features = c("Senescence"), order = TRUE) + fix.sc


generate_figs(umap_CARTEx_84, paste('./plots/', experiment, '_umap_CARTEx_84', sep = ''))
generate_figs(umap_sig_activation, paste('./plots/', experiment, '_umap_sig_activation', sep = ''))
generate_figs(umap_sig_anergy, paste('./plots/', experiment, '_umap_sig_anergy', sep = ''))
generate_figs(umap_sig_stemness, paste('./plots/', experiment, '_umap_sig_stemness', sep = ''))
generate_figs(umap_sig_senescence, paste('./plots/', experiment, '_umap_sig_senescence', sep = ''))

####################################################################################################
############################################ More plots ###########################################
####################################################################################################

scatterplot_CARTEx_activation <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Activation', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_anergy <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Anergy', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_stemness <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Stemness', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_senescence <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'Senescence', shuffle = TRUE, seed = 123)
scatterplot_CARTEx_allstates <- CombinePlots(plots = list(scatterplot_CARTEx_activation, scatterplot_CARTEx_anergy, scatterplot_CARTEx_stemness, scatterplot_CARTEx_senescence))

generate_figs(scatterplot_CARTEx_activation, paste('./plots/', experiment, '_scatterplot_CARTEx_activation', sep = ''))
generate_figs(scatterplot_CARTEx_anergy, paste('./plots/', experiment, '_scatterplot_CARTEx_anergy', sep = ''))
generate_figs(scatterplot_CARTEx_stemness, paste('./plots/', experiment, '_scatterplot_CARTEx_stemness', sep = ''))
generate_figs(scatterplot_CARTEx_senescence, paste('./plots/', experiment, '_scatterplot_CARTEx_senescence', sep = ''))
generate_figs(scatterplot_CARTEx_allstates, paste('./plots/', experiment, '_scatterplot_CARTEx_allstates', sep = ''), c(7,5))


CARTEx_630 <- rownames(read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1))
CARTEx_84 <- rownames(read.csv("../../weights/cartex-84-weights.csv", header = TRUE, row.names = 1))
TSG <- rownames(read.csv("../../signatures/tumor-suppressor-genes.csv", row.names = 1, header = TRUE))


heatmap_allgenes <- DoHeatmap(CART.combined.CD8pos, group.by = "characteristics_response", group.bar = TRUE, angle = 0)
generate_figs(heatmap_allgenes, paste('./plots/', experiment, '_heatmap_allgenes', sep = ''), c(7,10))

heatmap_CARTEx <- DoHeatmap(CART.combined.CD8pos, features = CARTEx_84, group.by = "characteristics_response", group.bar = TRUE, angle = 0)
generate_figs(heatmap_CARTEx, paste('./plots/', experiment, '_heatmap_CARTEx', sep = ''), c(7,10))

heatmap_TSG <- DoHeatmap(CART.combined.CD8pos, features = intersect(TSG, CARTEx_630), group.by = "characteristics_response", group.bar = TRUE, angle = 0)
generate_figs(heatmap_TSG, paste('./plots/', experiment, '_heatmap_TSG', sep = ''), c(7,10))



# Correlations of canonical exhaustion markers

canonical_exhaustion <- c('PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'TCF7', 'TOX', 'GZMB', 'IFNG', 'IL10', 'PRDM1', 'TNFRSF9')

dotplot_canonical_exhaustion <- DotPlot(CART.combined.CD8pos, features = canonical_exhaustion) + RotatedAxis()
generate_figs(dotplot_canonical_exhaustion, paste('./plots/', experiment, '_dotplot_canonical_exhaustion', sep = ''), c(10,4))

scatterplot_CARTEx_PDCD1 <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'PDCD1', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_PDCD1, paste('./plots/', experiment, '_scatterplot_CARTEx_PDCD1', sep = ''))

scatterplot_CARTEx_TCF7 <- FeatureScatter(CART.combined.CD8pos, feature1 = 'CARTEx_84', feature2 = 'TCF7', shuffle = TRUE, seed = 123)
generate_figs(scatterplot_CARTEx_TCF7, paste('./plots/', experiment, '_scatterplot_CARTEx_TCF7', sep = ''))

featureplot_canonical_exhaustion <- FeaturePlot(CART.combined.CD8pos, features = c('CARTEx_84', canonical_exhaustion), order = TRUE, ncol = 3)
generate_figs(featureplot_canonical_exhaustion, paste('./plots/', experiment, '_featureplot_canonical_exhaustion', sep = ''))


# Tumor suppressor genes

dotplot_TSG <- DotPlot(CART.combined.CD8pos, features = intersect(TSG, CARTEx_630)) + RotatedAxis()
generate_figs(dotplot_TSG, paste('./plots/', experiment, '_dotplot_TSG', sep = ''), c(18,4))


# CARTEx violin plot

vlnplot_CARTEx_84_response <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'characteristics_response', y.max = 6) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Responder','Non-responder')), label = "p.signif", label.y = 5)
generate_figs(vlnplot_CARTEx_84_response, paste('./plots/', experiment, '_vlnplot_CARTEx_84_response', sep = ''))

vlnplot_CARTEx_84_therapy <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'characteristics_therapy', y.max = 6) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('anti-CTLA4', 'anti-PD1')), label = "p.signif", label.y = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('anti-CTLA4+PD1', 'anti-PD1')), label = "p.signif", label.y = 4.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('anti-CTLA4+PD1', 'anti-CTLA4')), label = "p.signif", label.y = 5.5)
generate_figs(vlnplot_CARTEx_84_therapy, paste('./plots/', experiment, '_vlnplot_CARTEx_84_therapy', sep = ''))

vlnplot_CARTEx_84_therapy_splitby_responder <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'characteristics_therapy', y.max = 5, split.by = 'characteristics_response', split.plot = TRUE) + 
  theme(legend.position = 'top') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_therapy_splitby_responder, paste('./plots/', experiment, '_vlnplot_CARTEx_84_therapy_splitby_responder', sep = ''))

vlnplot_CARTEx_84_responder_splitby_timepoint <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'characteristics_response', y.max = 5, split.by = 'Timepoint', split.plot = TRUE) + 
  theme(legend.position = 'top') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_responder_splitby_timepoint, paste('./plots/', experiment, '_vlnplot_CARTEx_84_responder_splitby_timepoint', sep = ''))

vlnplot_CARTEx_84_timepoint_splitby_responder <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'Timepoint', y.max = 5, split.by = 'characteristics_response', split.plot = TRUE) + 
  theme(legend.position = 'top') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_timepoint_splitby_responder, paste('./plots/', experiment, '_vlnplot_CARTEx_84_timepoint_splitby_responder', sep = ''))

vlnplot_CARTEx_84_therapy_splitby_timepoint <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'characteristics_therapy', y.max = 5, split.by = 'Timepoint', split.plot = TRUE) + 
  theme(legend.position = 'top') + geom_boxplot(width=0.2, color="black", alpha=0)
generate_figs(vlnplot_CARTEx_84_therapy_splitby_timepoint, paste('./plots/', experiment, '_vlnplot_CARTEx_84_therapy_splitby_timepoint', sep = ''))


vlnplot_NKlike_Tex_response <- VlnPlot(CART.combined.CD8pos, features = c("NKlike_Tex"), group.by = 'characteristics_response', y.max = 6) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Responder','Non-responder')), label = "p.signif", label.y = 5)
generate_figs(vlnplot_NKlike_Tex_response, paste('./plots/', experiment, '_vlnplot_NKlike_Tex_response', sep = ''))







# Statistical test for therapy x response
# http://sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

obj.therapy <- SplitObject(CART.combined.CD8pos, split.by = "characteristics_therapy")
CTLA4_R <- as.numeric(unlist(subset(obj.therapy$`anti-CTLA4`, subset = characteristics_response == "Responder")[['CARTEx_84']]))
CTLA4_NR <- as.numeric(unlist(subset(obj.therapy$`anti-CTLA4`, subset = characteristics_response == "Non-responder")[['CARTEx_84']]))
CTLA4_res <- wilcox.test(CTLA4_R, CTLA4_NR, alternative = "two.sided")
PD1_R <- as.numeric(unlist(subset(obj.therapy$`anti-PD1`, subset = characteristics_response == "Responder")[['CARTEx_84']]))
PD1_NR <- as.numeric(unlist(subset(obj.therapy$`anti-PD1`, subset = characteristics_response == "Non-responder")[['CARTEx_84']]))
PD1_res <- wilcox.test(PD1_R, PD1_NR, alternative = "two.sided")
CTLA4_PD1_R <- as.numeric(unlist(subset(obj.therapy$`anti-CTLA4+PD1`, subset = characteristics_response == "Responder")[['CARTEx_84']]))
CTLA4_PD1_NR <- as.numeric(unlist(subset(obj.therapy$`anti-CTLA4+PD1`, subset = characteristics_response == "Non-responder")[['CARTEx_84']]))
CTLA4_PD1_res <- wilcox.test(CTLA4_PD1_R, CTLA4_PD1_NR, alternative = "two.sided")

stars.pval(c(CTLA4_res$p.value, PD1_res$p.value, CTLA4_PD1_res$p.value))

# alternative
# https://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html

obj.response <- SplitObject(CART.combined.CD8pos, split.by = "characteristics_response")
R_pre <- as.numeric(unlist(subset(obj.response$'Responder', subset = Timepoint == 'Pre')[['CARTEx_84']]))
R_post <- as.numeric(unlist(subset(obj.response$'Responder', subset = Timepoint == 'Post')[['CARTEx_84']]))
R_res <- wilcox.test(R_pre, R_post, alternative = "two.sided")
NR_pre <- as.numeric(unlist(subset(obj.response$'Non-responder', subset = Timepoint == 'Pre')[['CARTEx_84']]))
NR_post <- as.numeric(unlist(subset(obj.response$'Non-responder', subset = Timepoint == 'Post')[['CARTEx_84']]))
NR_res <- wilcox.test(NR_pre, NR_post, alternative = "two.sided")

stars.pval(c(R_res$p.value, NR_res$p.value))

# examine just pre (baseline)
obj.timepoint <- SplitObject(CART.combined.CD8pos, split.by = "Timepoint")
vlnplot_CARTEx_84_response_baseline <- VlnPlot(obj.timepoint$Pre, features = c("CARTEx_84"), group.by = 'characteristics_response', y.max = 5) + 
  theme(legend.position = 'none') + geom_boxplot(width=0.2, color="black", alpha=0) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Non-responder', 'Responder')), label = "p.signif", label.y = 4.5)
generate_figs(vlnplot_CARTEx_84_response_baseline, paste('./plots/', experiment, '_vlnplot_CARTEx_84_response_baseline', sep = ''))


# plot CARTEx scores by patient, organized by responder/non-responder

# CART.combined.CD8pos <- SetIdent(CART.combined.CD8pos, value = "characteristics_response")
vlnplot_patients_R <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'orig.ident', idents = c('Responder'), y.max = 6) + 
  theme(legend.position = 'none')
vlnplot_patients_NR <- VlnPlot(CART.combined.CD8pos, features = c("CARTEx_84"), group.by = 'orig.ident', idents = c('Non-responder'), y.max = 6) + 
  theme(legend.position = 'none')

vlnplot_patients_response <- CombinePlots(plots = list(vlnplot_patients_R, vlnplot_patients_NR), legend = 'none', nrow = length(x = 1))
generate_figs(vlnplot_patients_response, paste('./plots/', experiment, '_vlnplot_patients_response', sep = ''), c(20,5))


# aggregate gene expression by patient
modmetadata <- read.csv("./data/GSE120575-patient-data-formatted-modified.csv", row.names=1, header=TRUE)
df_express_patients <- AggregateExpression(CART.combined.CD8pos, group.by = 'orig.ident')
head(df_express_patients$RNA[c('CCL4'),])

# aggregate CARTEx by patient
obj.patient <- SplitObject(CART.combined.CD8pos, split.by = "orig.ident")

library(stringr)
str_extract(rownames(modmetadata), "[^_]+")


# COLORFUL gradient

CARTEx_score <- CART.combined.CD8pos@meta.data$CARTEx_84
Response <- CART.combined.CD8pos@meta.data$characteristics_response
Therapy <- CART.combined.CD8pos@meta.data$characteristics_therapy
Timepoint <- CART.combined.CD8pos@meta.data$Timepoint

meta_trunc <- data.frame(CARTEx_score, Response, Therapy, Timepoint)

gradient_timepoint_response <- ggplot(meta_trunc, aes(x=factor(Timepoint, level = c("Pre","Post")),y=CARTEx_score)) + 
  geom_boxplot() + geom_jitter(aes(color = Response), alpha = 0.5) + xlab("Therapy") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Pre','Post')), label = "p.signif")

generate_figs(gradient_timepoint_response, paste('./plots/', experiment, '_gradient_timepoint_response', sep = ''))










