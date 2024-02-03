

# https://github.com/kstreet13/slingshot/issues/43
# http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/Slingshot
# https://danhdtruong.com/Diffusion-Pseudotime/
# https://bioconductor.org/packages/release/bioc/vignettes/destiny/inst/doc/DPT.html
# https://rnabioco.github.io/cellar/previous/2019/docs/5_trajectories.html
# http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/diffusionMaps
# https://bioconductor.org/packages/release/bioc/vignettes/destiny/inst/doc/DPT.html
# https://bioconductor.org/packages/release/bioc/vignettes/destiny/inst/doc/Diffusion-Maps.html
# https://rnabioco.github.io/cellar/previous/2019/docs/3_norm_viz_clustering.html

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
dm <- readRDS(paste('./data/', experiment, '_dmap.rds', sep = ''))



md <- expt.obj@meta.data %>% as.data.table

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
transitplot_group_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = Group, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc
generate_figs(transitplot_group_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_group_CARTEx_84', sep = ''), c(8, 5))

transitplot_group_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = Group, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc
generate_figs(transitplot_group_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_group_CARTEx_200', sep = ''), c(8, 5))

transitplot_timepoint_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = TimePoint, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc
generate_figs(transitplot_timepoint_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_timepoint_CARTEx_200', sep = ''), c(8, 5))

transitplot_patient_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = Patient, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc
generate_figs(transitplot_patient_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_patient_CARTEx_200', sep = ''), c(8, 5))


# check how timepoint is distributed along diffusion map

timepoint.sc <- colorRampPalette(c("lightgrey","lightblue","mediumblue"))(length(unique(expt.obj@meta.data$TimePoint)))
transitplot_patient_timepoint <- ggplot(md, aes(x = DC1rank, y = Patient, colour = TimePoint)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + scale_color_manual(values = timepoint.sc)
generate_figs(transitplot_patient_timepoint, paste('./plots/', experiment, '_transition_transitplot_patient_timepoint', sep = ''), c(8, 5))

transitplot_patient_group <- ggplot(md, aes(x = DC1rank, y = Patient, colour = Group)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + scale_color_manual(values = c("red", "violetred", "violet", "purple"))
generate_figs(transitplot_patient_group, paste('./plots/', experiment, '_transition_transitplot_patient_group', sep = ''), c(8, 5))


transitplot_timepoint_patient <- ggplot(md, aes(x = DC1rank, y = TimePoint, colour = Patient)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1)
generate_figs(transitplot_timepoint_patient, paste('./plots/', experiment, '_transition_transitplot_timepoint_patient', sep = ''), c(8, 5))

transitplot_timepoint_phase <- ggplot(md, aes(x = DC1rank, y = TimePoint, colour = Phase)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + scale_color_manual(values = c("#089392", "#EAE29C", "#CF597E"))
generate_figs(transitplot_timepoint_phase, paste('./plots/', experiment, '_transition_transitplot_timepoint_phase', sep = ''), c(8, 5))

transitplot_timepoint_monaco <- ggplot(md, aes(x = DC1rank, y = TimePoint, colour = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + scale_color_manual(values = c("deepskyblue", "seagreen", "darkgoldenrod", "plum3"))
generate_figs(transitplot_timepoint_monaco, paste('./plots/', experiment, '_transition_transitplot_timepoint_monaco', sep = ''), c(9, 5))

