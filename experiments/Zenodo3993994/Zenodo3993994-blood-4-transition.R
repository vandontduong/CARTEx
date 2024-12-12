
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'Zenodo3993994'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_blood_scored.rds', sep = ''))


# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))


transitplot_disease_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = Disease, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 30000)
generate_figs(transitplot_disease_CARTEx_84, paste('./plots/', experiment, '_blood_transition_transitplot_disease_CARTEx_84', sep = ''), c(3.5,2))

transitplot_disease_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = Disease, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 30000) +
  scale_y_discrete(labels = c('P', 'C'))
generate_figs(transitplot_disease_CARTEx_200, paste('./plots/', experiment, '_blood_transition_transitplot_disease_CARTEx_200', sep = ''), c(3.8,2))

transitplot_disease_monaco <- ggplot(md, aes(x = DC1rank, y = Disease, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 30000) +
  # theme(legend.position="none") +
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  scale_y_discrete(labels = c('P', 'C'))
generate_figs(transitplot_disease_monaco, paste('./plots/', experiment, '_blood_transition_transitplot_disease_monaco', sep = ''), c(3.5,2))


umap_DC1rank <- FeaturePlot(expt.obj, feature = 'DC1rank', order = FALSE, pt.size = 0.1) + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c("yellow", "orange", "red", "violetred", "darkorchid"))
generate_figs(umap_DC1rank, paste('./plots/', experiment, '_blood_prepare_umap_DC1rank', sep = ''), c(3, 2))


# c("royalblue","palevioletred","darkorchid")









