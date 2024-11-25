
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE207935'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))


# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

ggplot(md, aes(x = dpt, y = AffstatStim, colour = AffstatStim)) +
  geom_quasirandom(groupOnX = FALSE)

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))


transitplot_affstatstim_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_affstatstim_CARTEx_84', sep = ''), c(10, 5))

transitplot_affstatstim_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_affstatstim_CARTEx_200', sep = ''), c(10, 5))

transitplot_affstatstim_CARTEx_630 <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = CARTEx_630)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_CARTEx_630, paste('./plots/', experiment, '_transition_transitplot_affstatstim_CARTEx_630', sep = ''), c(10, 5))

transitplot_affstatstim_activation <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Activation)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_activation, paste('./plots/', experiment, '_transition_transitplot_affstatstim_activation', sep = ''), c(10, 5))

transitplot_affstatstim_anergy <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Anergy)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_anergy, paste('./plots/', experiment, '_transition_transitplot_affstatstim_anergy', sep = ''), c(10, 5))

transitplot_affstatstim_senescence <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Senescence)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_senescence, paste('./plots/', experiment, '_transition_transitplot_affstatstim_senescence', sep = ''), c(10, 5))

transitplot_affstatstim_stemness <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Stemness)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_stemness, paste('./plots/', experiment, '_transition_transitplot_affstatstim_stemness', sep = ''), c(10, 5))

transitplot_affstatstim_LCMV_Tex <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = LCMV_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_LCMV_Tex, paste('./plots/', experiment, '_transition_transitplot_affstatstim_LCMV_Tex', sep = ''), c(10, 5))

transitplot_affstatstim_NKlike_Tex <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = NKlike_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_NKlike_Tex, paste('./plots/', experiment, '_transition_transitplot_affstatstim_NKlike_Tex', sep = ''), c(10, 5))

transitplot_affstatstim_BBD_Tex <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = BBD_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_BBD_Tex, paste('./plots/', experiment, '_transition_transitplot_affstatstim_BBD_Tex', sep = ''), c(10, 5))

transitplot_affstatstim_PD1_Tex <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = PD1_Tex)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_PD1_Tex, paste('./plots/', experiment, '_transition_transitplot_affstatstim_PD1_Tex', sep = ''), c(10, 5))

transitplot_affstatstim_TSR <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = TSR)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_TSR, paste('./plots/', experiment, '_transition_transitplot_affstatstim_TSR', sep = ''), c(10, 5))


transitplot_affstatstim_Tex_Term <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Tex_Term)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_Tex_Term, paste('./plots/', experiment, '_transition_transitplot_affstatstim_Tex_Term', sep = ''), c(10, 5))

transitplot_affstatstim_Tex_KLR <- ggplot(md, aes(x = DC1rank, y = AffstatStim, colour = Tex_KLR)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc + theme_classic()
generate_figs(transitplot_affstatstim_Tex_KLR, paste('./plots/', experiment, '_transition_transitplot_affstatstim_Tex_KLR', sep = ''), c(10, 5))



###










