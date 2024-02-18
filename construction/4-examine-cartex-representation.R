#!/usr/bin/env Rscript

### Script name: 4-examine-cartex-representation.R
### Description: 
### Author: Vandon Duong

# track time
ptm <- proc.time()
print("The script is starting...")

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
setwd(PATH_CONSTRUCTION)



cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)

TSG <- rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = TRUE))


expt.list <- c("GSE125881", "GSE136874", "GSE136184", "GSE120575", "GSE235676", "GSE151511", "mdandersonTCM")

# create dictionary hash containing all genes for each experiment in expt.list
expt.allgenes <- c()
expt.CARTEx630 <- c()
expt.CARTEx200 <- c()
expt.CARTEx84 <- c()
expt.TSG <- c()
expt.allgenes.counts <- c()
expt.CARTEx630.counts <- c()
expt.CARTEx200.counts <- c()
expt.CARTEx84.counts <- c()
expt.TSG.counts <- c()

for (exptID in expt.list){
  print(exptID)
  expt.allgenes[[exptID]] <- read.csv(paste(PATH_EXPERIMENTS, exptID, "/data/", exptID, "_allgenes.csv", sep = ''), header = TRUE, row.names = 1)$x
  expt.allgenes.counts[[exptID]] <- length(expt.allgenes[[exptID]])
  print(paste0("all genes: ", expt.allgenes.counts[[exptID]]))
  
  expt.CARTEx630[[exptID]] <- intersect(expt.allgenes[[exptID]], rownames(cartex_630_weights))
  expt.CARTEx200[[exptID]] <- intersect(expt.allgenes[[exptID]], rownames(cartex_200_weights))
  expt.CARTEx84[[exptID]] <- intersect(expt.allgenes[[exptID]], rownames(cartex_84_weights))
  
  expt.CARTEx630.counts[[exptID]] <- length(expt.CARTEx630[[exptID]])
  expt.CARTEx200.counts[[exptID]] <- length(expt.CARTEx200[[exptID]])
  expt.CARTEx84.counts[[exptID]] <- length(expt.CARTEx84[[exptID]])
  
  print(paste0("CARTEx630 genes: ", expt.CARTEx630.counts[[exptID]]))
  print(paste0("CARTEx200 genes: ", expt.CARTEx200.counts[[exptID]]))
  print(paste0("CARTEx84 genes: ", expt.CARTEx84.counts[[exptID]]))
  
  expt.TSG[[exptID]] <- intersect(expt.allgenes[[exptID]], TSG)
  expt.TSG.counts[[exptID]] <- length(expt.TSG[[exptID]])
  print(paste0("all genes: ", expt.TSG.counts[[exptID]]))
  
  print("")
}

expt.matrix.counts <- do.call("cbind", list(expt.allgenes.counts, expt.CARTEx630.counts, expt.CARTEx200.counts, expt.CARTEx84.counts, expt.TSG.counts))
expt.matrix.counts

write.csv(expt.matrix.counts, paste0(PATH_CONSTRUCTION, "data/expt_matrix_counts.csv"), row.names=TRUE)


CARTEx_versions <- list(
  CARTEx_630 = rownames(cartex_630_weights),
  CARTEx_200 = rownames(cartex_200_weights),
  CARTEx_84 = rownames(cartex_84_weights)
)

# upset_CARTEx_versions <- as.grob(upset(fromList(CARTEx_versions), order.by = "freq"))


# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# https://github.com/krassowski/complex-upset



upset_expt_allgenes <- upset(fromList(expt.allgenes), intersect = names(expt.allgenes))
generate_figs(upset_expt_allgenes, "./plots/upset_expt_allgenes", c(16, 6))

upset_expt_CARTEx630 <- upset(fromList(expt.CARTEx630), intersect = names(expt.CARTEx630))
generate_figs(upset_expt_CARTEx630, "./plots/upset_expt_CARTEx630", c(6, 6))

upset_expt_CARTEx200 <- upset(fromList(expt.CARTEx200), intersect = names(expt.CARTEx200))
generate_figs(upset_expt_CARTEx200, "./plots/upset_expt_CARTEx200", c(6, 6))

upset_expt_CARTEx84 <- upset(fromList(expt.CARTEx84), intersect = names(expt.CARTEx84))
generate_figs(upset_expt_CARTEx84, "./plots/upset_expt_CARTEx84", c(6, 6))

upset_expt_TSG <- upset(fromList(expt.TSG), intersect = names(expt.TSG))
generate_figs(upset_expt_TSG, "./plots/upset_expt_TSG", c(6, 6))



