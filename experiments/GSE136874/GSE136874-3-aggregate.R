# aggregate seurat object
# Vandon Duong

set.seed(123)

experiment = 'GSE136874'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(data.table)
library(EnhancedVolcano)
library(SingleR)
library(scuttle)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/sc-mode.html#pseudo-bulk-aggregation
# https://rdrr.io/bioc/SingleR/man/aggregateReference.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html


query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

# integration.obj@meta.data[["split.ident"]] <- "Query"
# integration.obj@meta.data$split.ident[integration.obj@meta.data$identifier == "GSE164378"] <- "Reference"
# integration.list <- SplitObject(integration.obj, split.by = "split.ident")
# integration.list$Query@meta.data[["cluster_id"]] <- factor(integration.list$Query@meta.data$identifier2)

# sce.obj <- SingleCellExperiment(assays = list(counts = integration.list$Query@assays$RNA@counts ), colData = integration.list$Query@meta.data)
sce.obj <- SingleCellExperiment(assays = list(counts = query.obj@assays$RNA@counts ), colData = query.obj@meta.data)
sce.obj <- logNormCounts(sce.obj)

aggregated.n1.obj <- aggregateReference(sce.obj, sce.obj$identifier2, power = 0)
aggregated.n5.obj <- aggregateReference(sce.obj, sce.obj$identifier2, ncenters = 5)

aggregated.n1.obj@assays@data$logcounts
aggregated.n5.obj@assays@data$logcounts


# CARTEx with weights // 630 genes
cartex_630_weights <- read.csv("../../weights/cartex-630-weights.csv", header = TRUE, row.names = 1)
common <- intersect(rownames(cartex_630_weights), rownames(query.obj))
expr <- t(as.matrix(GetAssayData(query.obj))[match(common, rownames(as.matrix(GetAssayData(query.obj)))),])
weights <- cartex_630_weights[match(common, rownames(cartex_630_weights)),]
scores <- expr %*% as.matrix(weights)
query.obj@meta.data$CARTEx_630 <- Z(scores)
query.obj@meta.data$CARTEx_630i <- integerize(query.obj@meta.data$CARTEx_630)
