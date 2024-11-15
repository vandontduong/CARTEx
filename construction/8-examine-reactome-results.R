
setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")


cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

library(readxl)
library(stringr)
library(ggvenn)

reactome_CARTEx_200 <- read.csv(paste0('./data/', 'reactome_result_CARTEx_200.csv'), header = TRUE)
reactome_CARTEx_630 <- read.csv(paste0('./data/', 'reactome_result_CARTEx_630.csv'), header = TRUE)
reactome_NK_like <- read.csv(paste0('./data/', 'reactome_result_NK_like.csv'), header = TRUE)

reactome_venn_overall <- list(
  CARTEx_200 = reactome_CARTEx_200$Pathway.name,
  CARTEx_630 = reactome_CARTEx_630$Pathway.name,
  NK_like = reactome_NK_like$Pathway.name
)

ggvenn(reactome_venn_overall)




select_significant_pathways <- function(reactome_results, threshold){
  return(reactome_results[reactome_results$Entities.pValue < threshold, ]$Pathway.name)
}

reactome_CARTEx_200_sig_pathways <- select_significant_pathways(reactome_CARTEx_200, 0.05)
reactome_CARTEx_630_sig_pathways <- select_significant_pathways(reactome_CARTEx_630, 0.05)
reactome_NK_like_sig_pathways <- select_significant_pathways(reactome_NK_like, 0.05)


reactome_venn_sig <- list(
  CARTEx_200 = reactome_CARTEx_200_sig_pathways,
  CARTEx_630 = reactome_CARTEx_630_sig_pathways,
  NK_like = reactome_NK_like_sig_pathways
)

ggvenn(reactome_venn_sig)



# create barchart of significant pathways in reactome for CARTEx (200)

reactome_CARTEx_200_sig_paths <- reactome_CARTEx_200 %>% filter(Entities.pValue <= 0.05 & Entities.FDR <= 0.05)

# entities found - number of curated molecules in pathway
barplot_significant_pathways <- ggplot(reactome_CARTEx_200_sig_paths, aes(x = reorder(Pathway.name, X.Entities.found), y = X.Entities.found)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_classic() + ylab('# of CARTEx genes in pathway') + xlab(NULL)

generate_figs(barplot_significant_pathways, "./plots/barplot_significant_pathways", c(6,2))















