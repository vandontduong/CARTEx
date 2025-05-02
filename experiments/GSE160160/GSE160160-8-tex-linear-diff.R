# testing different annotation methods

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
experiment = 'GSE160160'
setwd(paste(PATH_EXPERIMENTS, experiment, sep = ''))


expt.obj <- readRDS(paste('./data/', experiment, '_annotated.rds', sep = ''))

FeaturePlot(expt.obj, features = c("PDCD1","TCF7","CXCR5","CX3CR1","HAVCR2"))

annotate_Tex_cells <- function(seurat_obj) {
  # Extract expression data from Seurat object
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Get expression status (TRUE/FALSE) for each gene
  tcf7_expr <- expr_data["TCF7", ] > 0
  cxcr5_expr <- expr_data["CXCR5", ] > 0
  pdcd1_expr <- expr_data["PDCD1", ] > 0
  cx3cr1_expr <- expr_data["CX3CR1", ] > 0
  havcr2_expr <- expr_data["HAVCR2", ] > 0
  
  # Initialize annotation vector
  cell_annotations <- rep("Unclassified", ncol(seurat_obj))
  
  # Define Tex_term (PDCD1+ HAVCR2+)
  tex_term_cells <- pdcd1_expr & havcr2_expr
  cell_annotations[tex_term_cells] <- "Tex_term"
  
  # Define Tex_int (CX3CR1+ PDCD1+)
  tex_int_cells <- cx3cr1_expr & pdcd1_expr & !tex_term_cells
  cell_annotations[tex_int_cells] <- "Tex_int"
  
  # Define Tex_prog (TCF7+ CXCR5+ PDCD1+)
  tex_prog_cells <- tcf7_expr & cxcr5_expr & pdcd1_expr & !tex_term_cells & !tex_int_cells
  cell_annotations[tex_prog_cells] <- "Tex_prog"
  
  # Add annotations to Seurat object metadata
  seurat_obj$Tex_linear_diff <- cell_annotations
  
  # Return modified Seurat object
  return(seurat_obj)
}

expt.obj <- annotate_Tex_cells(expt.obj)


DimPlot(expt.obj, group.by = "Tex_linear_diff", shuffle = TRUE, seed = 123, cols = c('red', 'blue', 'grey'))



annotate_Tex_linear_diff <- function(seurat_obj) {
  # Extract expression data from the "data" slot
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Define the required genes
  required_genes <- c("PDCD1", "SLAMF6", "CX3CR1")
  
  # Check if all required genes are present in the Seurat object
  if (!all(required_genes %in% rownames(expr_data))) {
    stop("Not all required genes (PDCD1, SLAMF6, CX3CR1) are present in the Seurat object.")
  }
  
  # Create logical vectors for gene expression (> 0 means expressed)
  pdcd1_expr <- expr_data["PDCD1", ] > 0
  slamf6_expr <- expr_data["SLAMF6", ] > 0
  cx3cr1_expr <- expr_data["CX3CR1", ] > 0
  
  # Initialize annotation vector with "Unclassified"
  annotations <- rep("Unclassified", ncol(seurat_obj))
  
  # Assign annotations based on the specified conditions
  annotations[pdcd1_expr & slamf6_expr & !cx3cr1_expr] <- "Tex_prog"
  annotations[pdcd1_expr & cx3cr1_expr & !slamf6_expr] <- "Tex_int"
  annotations[pdcd1_expr & !slamf6_expr & !cx3cr1_expr] <- "Tex_term"
  
  # Add the annotations to the Seurat object's metadata
  seurat_obj$Tex_linear_diff <- annotations
  
  # Return the modified Seurat object
  return(seurat_obj)
}


expt.obj <- annotate_Tex_linear_diff(expt.obj)
table(expt.obj$Tex_linear_diff)

DimPlot(expt.obj, group.by = "Tex_linear_diff", shuffle = TRUE, seed = 123, cols = c('green','red', 'blue', 'grey'))





FeaturePlot(expt.obj, features = c("PDCD1","TCF7","TOX","SLAMF6"))

annotate_Tex_linear_diff <- function(seurat_obj) {
  # Extract expression data from the "data" slot
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Define the required genes
  required_genes <- c("TCF7","TOX","SLAMF6")
  
  # Check if all required genes are present in the Seurat object
  if (!all(required_genes %in% rownames(expr_data))) {
    stop("Not all required genes (TCF7, TOX, SLAMF6) are present in the Seurat object.")
  }
  
  # Create logical vectors for gene expression (> 0 means expressed)
  tcf7_expr <- expr_data["TCF7", ] > 0
  slamf6_expr <- expr_data["SLAMF6", ] > 0
  tox_expr <- expr_data["TOX", ] > 0
  
  # Initialize annotation vector with "Unclassified"
  annotations <- rep("Unclassified", ncol(seurat_obj))
  
  # Assign annotations based on the specified conditions
  annotations[tcf7_expr & slamf6_expr & tox_expr] <- "Tex_prog"
  
  # Add the annotations to the Seurat object's metadata
  seurat_obj$Tex_linear_diff <- annotations
  
  # Return the modified Seurat object
  return(seurat_obj)
}



expt.obj <- annotate_Tex_linear_diff(expt.obj)

table(expt.obj$Tex_linear_diff)
DimPlot(expt.obj, group.by = "Tex_linear_diff", shuffle = TRUE, seed = 123, cols = c('red', 'grey'))





annotate_Tex_linear_diff <- function(seurat_obj) {
  # Extract expression data from the "data" slot
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Define the required genes
  required_genes <- c("PDCD1","IL7R")
  
  # Check if all required genes are present in the Seurat object
  if (!all(required_genes %in% rownames(expr_data))) {
    stop("Not all required genes (PDCD1, IL7R) are present in the Seurat object.")
  }
  
  # Create logical vectors for gene expression (> 0 means expressed)
  pdcd1_expr <- expr_data["PDCD1", ] > 0
  il7r_expr <- expr_data["IL7R", ] > 0
  
  # Initialize annotation vector with "Unclassified"
  annotations <- rep("Unclassified", ncol(seurat_obj))
  
  # Assign annotations based on the specified conditions
  annotations[pdcd1_expr & il7r_expr] <- "Tex_prog"
  
  # Add the annotations to the Seurat object's metadata
  seurat_obj$Tex_linear_diff <- annotations
  
  # Return the modified Seurat object
  return(seurat_obj)
}



expt.obj <- annotate_Tex_linear_diff(expt.obj)

table(expt.obj$Tex_linear_diff)
DimPlot(expt.obj, group.by = "Tex_linear_diff", shuffle = TRUE, seed = 123, cols = c('red', 'grey'))



FeaturePlot(expt.obj, features = c("PDCD1","IL7R"))






annotate_Tex_linear_diff <- function(seurat_obj) {
  # Extract expression data from the "data" slot
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Define the required genes
  required_genes <- c("TCF7","HAVCR2")
  
  # Check if all required genes are present in the Seurat object
  if (!all(required_genes %in% rownames(expr_data))) {
    stop("Not all required genes (TCF7, HAVCR2) are present in the Seurat object.")
  }
  
  # Create logical vectors for gene expression (> 0 means expressed)
  tcf7_expr <- expr_data["TCF7", ] > 0
  havcr2_expr <- expr_data["HAVCR2", ] > 0
  
  # Initialize annotation vector with "Unclassified"
  annotations <- rep("Unclassified", ncol(seurat_obj))
  
  # Assign annotations based on the specified conditions
  annotations[tcf7_expr & !havcr2_expr] <- "Tex_prog"
  annotations[!tcf7_expr & havcr2_expr] <- "Tex_term"
  
  # Add the annotations to the Seurat object's metadata
  seurat_obj$Tex_linear_diff <- annotations
  
  # Return the modified Seurat object
  return(seurat_obj)
}



expt.obj <- annotate_Tex_linear_diff(expt.obj)

table(expt.obj$Tex_linear_diff)
DimPlot(expt.obj, group.by = "Tex_linear_diff", shuffle = TRUE, seed = 123, cols = c('blue', 'red', 'grey'))

FeaturePlot(expt.obj, features = c("PDCD1","TCF7", "HAVCR2", "SLAMF6"))

