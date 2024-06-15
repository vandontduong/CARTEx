
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
experiment = 'GSE136184'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
##################################### Load single-cell dataset #####################################
####################################################################################################

expt.obj <- readRDS(paste('./data/', experiment, '_cs_scored.rds', sep = ''))

md <- expt.obj@meta.data %>% as.data.table
md[, .N, by = c("AgeGroup2", "Age")]


fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
transitplot_group_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = AgeGroup2, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1) + fix.sc
generate_figs(transitplot_group_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_group_CARTEx_200', sep = ''), c(8, 5))

transitplot_group_monaco <- ggplot(md, aes(x = DC1rank, y = AgeGroup2, colour = monaco)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.1)
generate_figs(transitplot_group_monaco, paste('./plots/', experiment, '_transition_transitplot_group_monaco', sep = ''), c(8, 5))





