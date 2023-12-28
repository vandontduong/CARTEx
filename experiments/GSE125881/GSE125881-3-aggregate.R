### Script name: GSE125881-3-aggregate.R
### Description: aggregate seurat object for GSE125881
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE125881'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))

####################################################################################################
######################################## Load data and filter ######################################
####################################################################################################

# https://bioconductor.org/books/release/SingleRBook/sc-mode.html#pseudo-bulk-aggregation
# https://rdrr.io/bioc/SingleR/man/aggregateReference.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

query.obj <- readRDS(paste('./data/', experiment, '_query_scored.rds', sep = ''))

# https://github.com/satijalab/seurat/issues/6007



# pseudo_samples <- 3
# num_cells <- length(Cells(expt.obj))
# labels <- sample.int(pseudo_samples, num_cells, replace = TRUE)

expt.obj@meta.data$pblabels <- PseudoBulkLabels(expt.obj, 3)

expt.obj.agg <- AggregateExpression(expt.obj, group.by = c('Group', 'pblabels'), return.seurat = TRUE)

# filter by cell type?
# use identifier instead of Group when aggregating query data

