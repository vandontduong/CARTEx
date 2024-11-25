
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

ggvenn(reactome_venn_sig, stroke_size = 0.5, set_name_size = 3, text_size = 3)



# create barchart of significant pathways in reactome for CARTEx (200)

reactome_CARTEx_200_sig_paths <- reactome_CARTEx_200 %>% filter(Entities.pValue <= 0.05 & Entities.FDR <= 0.05)

# entities found - number of curated molecules in pathway
barplot_significant_pathways <- ggplot(reactome_CARTEx_200_sig_paths, aes(x = reorder(Pathway.name, X.Entities.found), y = X.Entities.found)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_classic() + ylab('Entities found in pathway') + xlab(NULL)

generate_figs(barplot_significant_pathways, "./plots/barplot_significant_pathways", c(6,2))




### compare with STAT 3 GOF



reactome_STAT3_GOF_stim <- read.csv(paste0('./data/', 'reactome_result_STAT3_GOF_stim.csv'), header = TRUE)
reactome_STAT3_GOF_rest <- read.csv(paste0('./data/', 'reactome_result_STAT3_GOF_rest.csv'), header = TRUE)
reactome_STAT3_GOF_intersect <- read.csv(paste0('./data/', 'reactome_result_STAT3_GOF_intersect.csv'), header = TRUE)
reactome_STAT3_GOF_combined <- read.csv(paste0('./data/', 'reactome_result_STAT3_GOF_combined.csv'), header = TRUE)


reactome_STAT3_GOF_stim_sig_pathways <- select_significant_pathways(reactome_STAT3_GOF_stim, 0.05)
reactome_STAT3_GOF_rest_sig_pathways <- select_significant_pathways(reactome_STAT3_GOF_rest, 0.05)
reactome_STAT3_GOF_intersect_sig_pathways <- select_significant_pathways(reactome_STAT3_GOF_intersect, 0.05)
reactome_STAT3_GOF_combined_sig_pathways <- select_significant_pathways(reactome_STAT3_GOF_combined, 0.05)


reactome_venn_sig_STAT3 <- list(
  STAT3_GOF_stim = reactome_STAT3_GOF_stim_sig_pathways,
  STAT3_GOF_rest = reactome_STAT3_GOF_rest_sig_pathways,
  STAT3_GOF_combined = reactome_STAT3_GOF_combined_sig_pathways
)

ggvenn(reactome_venn_sig_STAT3, stroke_size = 0.5, set_name_size = 3, text_size = 3)

reactome_venn_sig_CARTEx_STAT3 <- list(
  CARTEx = reactome_CARTEx_200_sig_pathways,
  STAT3_GOF = reactome_STAT3_GOF_combined_sig_pathways
)

ggvenn(reactome_venn_sig_CARTEx_STAT3, stroke_size = 0.5, set_name_size = 3, text_size = 3)

### 


# analyze entities 
reactome_CARTEx_200_sig_paths$Submitted.entities.found




# create barchart of significant pathways in reactome for STAT3 GOF

reactome_STAT3_GOF_combined_sig_paths <- reactome_STAT3_GOF_combined %>% filter(Entities.pValue <= 0.05 & Entities.FDR <= 0.05)

# entities found - number of curated molecules in pathway
# bars_to_highlight <- intersect(reactome_CARTEx_200_sig_paths$Pathway.name, reactome_STAT3_GOF_combined_sig_paths$Pathway.name)
# bars_to_highlight <- intersect(reactome_CARTEx_200$Pathway.name, reactome_STAT3_GOF_combined$Pathway.name)
bars_to_highlight <- intersect((reactome_CARTEx_200 %>% filter(Entities.pValue <= 0.05))$Pathway.name, (reactome_STAT3_GOF_combined %>% filter(Entities.pValue <= 0.05))$Pathway.name)

barplot_significant_pathways_STAT3_GOF_combined <- ggplot(reactome_STAT3_GOF_combined_sig_paths, aes(x = reorder(Pathway.name, X.Entities.found), y = X.Entities.found)) + 
  geom_bar(stat = 'identity', aes(fill = ifelse(Pathway.name %in% bars_to_highlight, "CARTEx", "General"))) + coord_flip() + theme_classic() + ylab('Entities found in pathway') + xlab(NULL) +
  scale_fill_manual(values = c("CARTEx" = "indianred", "General" = "grey60"), name = "Pathway") + 
  theme(legend.position = c(0.85, 0.15))

generate_figs(barplot_significant_pathways_STAT3_GOF_combined, "./plots/barplot_significant_pathways_STAT3_GOF_combined", c(7,3.75))


## examine intersection instead of "combined"

reactome_STAT3_GOF_intersect_sig_paths <- reactome_STAT3_GOF_intersect %>% filter(Entities.pValue <= 0.05 & Entities.FDR <= 0.05)

# bars_to_highlight <- intersect(reactome_CARTEx_200_sig_paths$Pathway.name, reactome_STAT3_GOF_intersect_sig_paths$Pathway.name)
bars_to_highlight <- intersect((reactome_CARTEx_200 %>% filter(Entities.pValue <= 0.05))$Pathway.name, (reactome_STAT3_GOF_intersect %>% filter(Entities.pValue <= 0.05))$Pathway.name)
length(bars_to_highlight)

barplot_significant_pathways_STAT3_GOF_intersect <- ggplot(reactome_STAT3_GOF_intersect_sig_paths, aes(x = reorder(Pathway.name, X.Entities.found), y = X.Entities.found)) + 
  geom_bar(stat = 'identity', aes(fill = ifelse(Pathway.name %in% bars_to_highlight, "CARTEx", "General"))) + coord_flip() + theme_classic() + ylab('Entities found in pathway') + xlab(NULL) +
  scale_fill_manual(values = c("CARTEx" = "indianred", "General" = "grey60"), name = "Pathway") + 
  theme(legend.position = c(0.85, 0.15))

generate_figs(barplot_significant_pathways_STAT3_GOF_intersect, "./plots/barplot_significant_pathways_STAT3_GOF_intersect", c(7,5.5))


















