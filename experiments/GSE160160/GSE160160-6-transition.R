
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))

# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

ggplot(md, aes(x = dpt, y = exposure, colour = exposure)) +
  geom_quasirandom(groupOnX = FALSE)

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))

transitplot_exposure_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = exposure, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc + theme_classic() + scale_x_continuous(breaks = seq(0, 14000, by = 2000))
generate_figs(transitplot_exposure_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_exposure_CARTEx_84', sep = ''), c(10, 5))

transitplot_exposure_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = exposure, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc + theme_classic() + scale_x_continuous(breaks = seq(0, 14000, by = 2000))
generate_figs(transitplot_exposure_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_exposure_CARTEx_200', sep = ''), c(10, 5))

transitplot_exposure_LCMV_Tex <- ggplot(md, aes(x = DC1rank, y = exposure, colour = LCMV_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc + theme_classic() + scale_x_continuous(breaks = seq(0, 14000, by = 2000))
generate_figs(transitplot_exposure_LCMV_Tex, paste('./plots/', experiment, '_transition_transitplot_exposure_LCMV_Tex', sep = ''), c(10, 5))

transitplot_exposure_NKlike_Tex <- ggplot(md, aes(x = DC1rank, y = exposure, colour = NKlike_Tex)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc + theme_classic() + scale_x_continuous(breaks = seq(0, 14000, by = 2000))
generate_figs(transitplot_exposure_NKlike_Tex, paste('./plots/', experiment, '_transition_transitplot_exposure_NKlike_Tex', sep = ''), c(10, 5))

dmap_exposure <- DimPlot(expt.obj, reduction = "dm", group.by = "exposure", shuffle = TRUE, seed = 123, cols = c('skyblue', 'cadetblue'))



