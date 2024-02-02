
####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE136874'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))



expt.obj <- readRDS(paste('./data/', experiment, '_scored.rds', sep = ''))
dm <- readRDS(paste('./data/', experiment, '_dmap.rds', sep = ''))
dpt <- DPT(dm)

df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], 
                 dptval = dpt$dpt, CAR = dpt$CAR)

expt.obj$dpt <- dpt$dpt

# https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html

md <- expt.obj@meta.data %>% as.data.table

ggplot(md, aes(x = dpt, y = CAR, colour = CAR)) +
  geom_quasirandom(groupOnX = FALSE)

fix.sc <- scale_color_gradientn(colours = c("blue","lightgrey","red"), limits = c(-4,4))
transitplot_CAR_CARTEx_84 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_CARTEx_84, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_84', sep = ''), c(8, 5))

transitplot_CAR_CARTEx_200 <- ggplot(md, aes(x = DC1rank, y = CAR, colour = CARTEx_200)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
generate_figs(transitplot_CAR_CARTEx_200, paste('./plots/', experiment, '_transition_transitplot_CAR_CARTEx_200', sep = ''), c(8, 5))


expt.obj$UMAP1 <- expt.obj[['umap']]@cell.embeddings[,1]
test <- ggplot(md, aes(x = UMAP1, y = CAR, colour = CARTEx_84)) +
  geom_quasirandom(groupOnX = FALSE) + fix.sc
