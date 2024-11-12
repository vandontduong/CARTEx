



####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE146264'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))


# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))


transitplot_severity_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = Severity, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 4000) +
  scale_y_discrete(labels = c('C', 'M', 'S'))
generate_figs(transitplot_severity_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_severity_CARTEx_84', sep = ''), c(3.5,2))

transitplot_severity_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = Severity, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 4000) +
  scale_y_discrete(labels = c('C', 'M', 'S'))
generate_figs(transitplot_severity_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_severity_CARTEx_200', sep = ''), c(3.5,2))



transitplot_severity_monaco <- ggplot(md, aes(x = DC1rank, y = Severity, color = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic() + xlim(0, 4000) +
  # theme(legend.position="none") +
  scale_color_manual(labels=c("N", "CM", "EM", "TE"), values = c('Naive CD8 T cells' = 'deepskyblue', 'Central memory CD8 T cells' = 'seagreen', 'Effector memory CD8 T cells' = 'darkgoldenrod', 'Terminal effector CD8 T cells' = 'plum3')) +
  scale_y_discrete(labels = c('C', 'M', 'S'))
generate_figs(transitplot_severity_monaco, paste('./plots/', experiment, '_transition_transitplot_severity_monaco', sep = ''), c(3.5,2))


FeaturePlot(expt.obj, feature = 'DC1rank')










