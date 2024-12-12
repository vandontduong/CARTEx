### Script name: cartex-utilities.R
### Description: load libraries, paths, and custom functions
### Author: Vandon Duong

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

# build package
# http://rcs.bu.edu/examples/r/tutorials/BuildingPackages/r-project-setup.html

# additional functions
# https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html

####################################################################################################
######################################### Library and paths ########################################
####################################################################################################

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)
library(SingleR)
library(scuttle)
library(destiny)
library(data.table)
library(patchwork)
library(EnhancedVolcano)
library(clustree)
library(ROGUE)
library(glmGamPoi)
library(ggbeeswarm)
# library(UpSetR) # need fromList() function
library(ComplexUpset)
library(ComplexHeatmap)
library(ca)
library(cowplot)

# absolute path to where the project directory resides
PATH_CARTEX <- '/oak/stanford/groups/cmackall/vandon/CARTEx/'

# relative paths
PATH_CELLANNOTATE <- paste(PATH_CARTEX, 'cellannotate/', sep = '')
PATH_CONSTRUCTION <- paste(PATH_CARTEX, 'construction/', sep = '')
PATH_EXPERIMENTS <- paste(PATH_CARTEX, 'experiments/', sep = '')
PATH_SIGNATURES <- paste(PATH_CARTEX, 'signatures/', sep = '')
PATH_WEIGHTS <- paste(PATH_CARTEX, 'weights/', sep = '')

# file calling
# paste(<PATH_NAME>, <FILE_NAME>, sep = '')

####################################################################################################
######################################### General functions ########################################
####################################################################################################

############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: Z score normalization
# @ s: score

Z <- function(s){
  s <- as.numeric(s)
  z <- (s - mean(s))/sd(s)
  return(z)
}




############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: generate_figs() saves figures in jpeg and pdf formats
# @ figure_object: an object that visualizes data
# @ file_name: desired output file name (excluding extension)
# @ dimensions: specify desired figure size using c(x,y)

generate_figs <- function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
    ggsave(filename = gsub(" ", "", paste(file_name,".png")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".png")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}






############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: round scores to nearest integer; cut-off at +/-4
# @ score: vector of floats

integerize <- function(score){
  score_mod <- round(score)
  score_mod[score_mod <= -4] <- -4
  score_mod[score_mod >= 4] <- 4
  return (score_mod)
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: check metadata levels() and return "NULL" if levels() does not exist
# need to update levels() for anything that results in "not NULL"
# @ atlas: Seurat object with metadata

check_levels <- function(atlas){
  for (i in colnames(atlas@meta.data)){
    j = levels(atlas@meta.data[[i]])
    if (is.null(j) == TRUE){
      print(paste(i, 'levels() is NULL'))
    } else {
      print(paste(i, 'levels() is not NULL'))
      print(j)
    }
    cat("\n")}
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: sort list and optionally shift the dataset

SortNumStrList <- function(num_str_list, shift){
  if(shift == TRUE){
    return(as.character(sort(as.numeric(num_str_list) - 1)))
  }
  return(as.character(sort(as.numeric(num_str_list))))
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: VlnPlot() customized to show cutoff thresholds for quality control filters
# only use on original dataset before QC filtering, otherwise it will not display the non-filtered datapoints
# @ atlas: Seurat object with features relevant for quality control
# @ metrics: list of features relevant for quality control
# @ low_cutoff: vector of lower-bound thresholds for filtering cells, in the order of metrics
# @ high_cutoff: vector of upper-bound thresholds for filtering cells, in the order of metrics
# @ identity: string describing the relevant metadata identity group
# @ ncols: number of columns to split the resulting figure into

ViolinPlotQC <- function(atlas, metrics, low_cutoff, high_cutoff, identity, ncols){
  plot.list <- list()
  counter = 0
  for (i in metrics){
    counter = counter + 1
    vlnplt <- VlnPlot(object = atlas, features = metrics[counter], group.by = identity) + NoLegend() + ggtitle(i) + scale_y_continuous(limits = c(0, NA))
    if (is.na(low_cutoff[counter]) == FALSE & is.na(high_cutoff[counter]) == FALSE){
      plot.list[[i]] <- vlnplt & geom_hline(yintercept = low_cutoff[counter], linetype='dashed', color=c('red')) & geom_hline(yintercept = high_cutoff[counter], linetype='dashed', color=c('red'))
    }
    if (is.na(low_cutoff[counter]) == FALSE & is.na(high_cutoff[counter]) == TRUE){
      plot.list[[i]] <- vlnplt & geom_hline(yintercept = low_cutoff[counter], linetype='dashed', color=c('red'))
    }
    if (is.na(low_cutoff[counter]) == TRUE & is.na(high_cutoff[counter]) == FALSE){
      plot.list[[i]] <- vlnplt & geom_hline(yintercept = high_cutoff[counter], linetype='dashed', color=c('red'))
    }
    if (is.na(low_cutoff[counter]) == TRUE & is.na(high_cutoff[counter]) == TRUE){
      plot.list[[i]] <- vlnplt
    }
  }
  # combined_plots <- Reduce(`+`, plot.list) + patchwork::plot_layout(ncol = ncols)
  combined_plots <- wrap_plots(plot.list, ncol = ncols)
  return(combined_plots)
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: ChoosePC() to select principal components to capture most of the variation
# @ atlas: Seurat object

ChoosePC <- function(atlas){
  # Determine percent of variation associated with each PC
  pct <- atlas[["pca"]]@stdev / sum(atlas[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # Minimum of the two calculation
  return(min(co1, co2))
}

# References
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html


############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: generate diffusion embeddings and attach to Seurat object
# @ atlas: Seurat object
# @ k_int: number of nearest neighbor to consider

RunDiffusion <- function(atlas, k_int){
  sce <- as.SingleCellExperiment(atlas)
  dm <- DiffusionMap(sce, k = k_int, verbose = TRUE)
  tmp <- data.matrix(data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], row.names = colnames(atlas)))
  atlas[["dm"]] <- CreateDimReducObject(embeddings = tmp, key = "DC_", assay = DefaultAssay(atlas))
  atlas@meta.data$DC1rank <- rank(atlas[['dm']]@cell.embeddings[,1])
  atlas@meta.data$DC2rank <- rank(atlas[['dm']]@cell.embeddings[,2])
  dpt_obj <- DPT(dm)
  atlas@meta.data$dpt <- dpt_obj@dm@eigenvec0
  return(list("atlas" = atlas, "dmap" = dm))
}

# http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/diffusionMaps
# https://rnabioco.github.io/cellar/previous/2019/docs/3_norm_viz_clustering.html
# https://github.com/satijalab/seurat/issues/1475
# https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette.html






############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: DimPlot() customized to highlight cells by identity group
# @ atlas: Seurat object with features relevant for quality control
# @ identity: string describing the relevant metadata identity group
# @ reduction_map: string describing the relevant dimensionality reduction method
# @ highlight_color: string describing the color of the highlighted datapoints
# @ pt_size: string describing the size of the datapoints
# @ ncols: number of columns to split the resulting figure into

DimPlotHighlightIdents <- function(atlas, identity, reduction_map, highlight_color, pt_size, ncols){
  plot.list <- list()
  for (i in sapply(unique(x = atlas[[deparse(substitute(identity))]]), levels)) {
    tryCatch({
      plot.list[[i]] <- DimPlot(
        object = atlas, reduction = reduction_map, raster = FALSE, cols.highlight = highlight_color, pt.size = pt_size, sizes.highlight = pt_size,
        cells.highlight = Cells(atlas[, atlas[[deparse(substitute(identity))]] == i])
      ) + NoLegend() + ggtitle(i)
    }, error = function(e){cat("ERROR :",conditionMessage(e), " for ", i, "\n")})
  }
  # combined_plots <- CombinePlots(plots = plot.list, ncol = ncols)
  combined_plots <- Reduce(`+`, plot.list) + patchwork::plot_layout(ncol = ncols)
  return(combined_plots)
}

# References
# https://github.com/satijalab/seurat/issues/1396
# https://stackoverflow.com/questions/27676404/list-all-factor-levels-of-a-data-frame
# https://bioinformatics.stackexchange.com/questions/18902/dimplot-how-to-highlight-cells-with-identity-colors

# Notes
# Here are two ways to capture cells which match; the latter works within defined function
# WhichCells(object = integration.obj, expression = identifier == experiment)
# Cells(integration.obj[, integration.obj[['identifier']] == experiment])







############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: correspondence analysis from contingency tables
# @ atlas: Seurat object
# @ identity: string describing the relevant metadata identity group
# @ reduction_map: string describing the relevant dimensionality reduction method
# @ highlight_color: string describing the color of the highlighted datapoints
# @ pt_size: string describing the size of the datapoints
# @ ncols: number of columns to split the resulting figure into

make.ca.plot.df <- function (ca.plot.obj,
                             row.lab = "Rows",
                             col.lab = "Columns") {
  df <- data.frame(Label = c(rownames(ca.plot.obj$rows),
                             rownames(ca.plot.obj$cols)),
                   Dim1 = c(ca.plot.obj$rows[,1], ca.plot.obj$cols[,1]),
                   Dim2 = c(ca.plot.obj$rows[,2], ca.plot.obj$cols[,2]),
                   Variable = c(rep(row.lab, nrow(ca.plot.obj$rows)),
                                rep(col.lab, nrow(ca.plot.obj$cols))))
  rownames(df) <- 1:nrow(df)
  df
}

CorrespondenceAnalysisPlot <- function(contingency_table, row_lab, col_lab){
  ca.fit <- ca(contingency_table)
  ca.plot <- plot(ca.fit)
  ca.plot.df <- make.ca.plot.df(ca.plot, row.lab = row_lab, col.lab = col_lab)
  ca.sum <- summary(ca.fit)
  dim.var.percs <- ca.sum$scree[,"values2"]
  
  p <- ggplot(ca.plot.df, aes(x = Dim1, y = Dim2, col = Variable, shape = Variable, label = Label)) +
    geom_vline(xintercept = 0, lty = "dashed", alpha = .5) + geom_hline(yintercept = 0, lty = "dashed", alpha = .5) + geom_point()
  
  p <- p +
    scale_x_continuous(limits = range(ca.plot.df$Dim1) + c(diff(range(ca.plot.df$Dim1)) * -0.2, diff(range(ca.plot.df$Dim1)) * 0.2)) +
    scale_y_continuous(limits = range(ca.plot.df$Dim2) + c(diff(range(ca.plot.df$Dim2)) * -0.2, diff(range(ca.plot.df$Dim2)) * 0.2)) +
    scale_size(range = c(4, 7), guide = F) +
    geom_label_repel(show.legend = F, segment.alpha = .5, point.padding = unit(5, "points"))
  
  p <- p +
    labs(x = paste0("Dimension 1 (", signif(dim.var.percs[1], 3), "%)"),
         y = paste0("Dimension 2 (", signif(dim.var.percs[2], 3), "%)"),
         col = "", shape = "") + theme_bw() + scale_colour_brewer(palette = "Set1") + theme(legend.position='bottom')
}

# References
# https://stats.stackexchange.com/questions/147721/which-is-the-best-visualization-for-contingency-tables
# https://www.r-bloggers.com/2019/08/correspondence-analysis-visualization-using-ggplot/







############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: clean up table and remove any rows or columns that are all zero values
# @ tab: table to clean up

CleanTable <- function(tab){
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
}




############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: heatmap
# @ atlas: Seurat object
# @ gene_list: gene list to select for

HeatmapWithMultiGroups <- function(atlas, gene_list, anno_colors, sort_by_group){
  meta_data <- atlas@meta.data
  # sort_by_group <- names(anno_colors)[anno_i]
  ordered_meta_data <- meta_data[order(meta_data[[sort_by_group]]), ]
  
  common <- intersect(gene_list, rownames(atlas))
  
  scale_data <- atlas@assays$RNA@scale.data
  scale_data <- scale_data[match(common, rownames(atlas)), rownames(ordered_meta_data)]
  
  ordered_meta_data <- select(ordered_meta_data, names(anno_colors))
  rownames(ordered_meta_data) <- NULL
  
  heatmap_anno = HeatmapAnnotation(df = ordered_meta_data,
                                   show_annotation_name = TRUE,
                                   col = anno_colors)
  
  col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  
  heatmap_plt <- Heatmap(scale_data, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                         column_order = NULL, show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE,
                         use_raster = TRUE, bottom_annotation = NULL, top_annotation = heatmap_anno)
  
  return(heatmap_plt)
}


# https://divingintogeneticsandgenomics.com/post/enhancement-of-scrnaseq-heatmap-using-complexheatmap/
# https://bioinformatics.stackexchange.com/questions/4851/how-i-can-reproduce-this-heat-map







############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: Bar plot stacked and split according to two different identity groups
# @ atlas: Seurat object
# @ x_identity: string describing the relevant metadata identity group to stack by
# @ y_identity: string describing the relevant metadata identity group to split by

BarPlotStackSplit <- function(atlas, x_identity, y_identity, color_set = NULL){
  md <- as.data.table(atlas@meta.data)[, .N, by = c(x_identity, y_identity)]
  md$percent <- round(100*md$N / sum(md$N), digits = 1)
  if (is.null(color_set)) {
    barplot <- ggplot(md, aes_string(x = y_identity, y = 'N', fill = x_identity)) + geom_col(position = "fill")
  }
  else {
    barplot <- ggplot(md, aes_string(x = y_identity, y = 'N', fill = x_identity)) + geom_col(position = "fill") +
      scale_fill_manual(values=color_set) + theme_classic() + theme(axis.title.x = element_blank()) + ylab('Fraction')
  }
  return(barplot)
}






############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: PercentageFeatureSet() customized to compute percentage of feature set detected, i.e. expression level is non-zero, rather than what percentage it is of all transcripts
# Examine CARTEx representation at single-cell resolution
# @ atlas: Seurat object
# @ feature_set: gene set to evaluate on

PercentageFeatureSetDetected <- function(atlas, feature_set){
  all.genes <- rownames(atlas)
  exp.mat <- atlas[["RNA"]]$counts
  exp.mat <- exp.mat[which(all.genes %in% feature_set),]
  PFSD <- diff(exp.mat@p) / length(feature_set) * 100
  return(PFSD)
}

# Notes
# equivalent
# intersect(all.genes, rownames(cartex_630_weights))
# rownames(cartex_630_weights)[rownames(cartex_630_weights) %in% all.genes]

# Regular PFS generates multiple duplicate vectors appended to each other, so need to subset
# expt.obj@meta.data$percent.CARTEx_630 <- PercentageFeatureSet(expt.obj, features = intersect(all.genes, rownames(cartex_630_weights)), assay = 'RNA')[1:length(Cells(expt.obj))]
# expt.obj@meta.data$PFSD_CARTEx_630 <- PercentageFeatureSetDetected(expt.obj, rownames(cartex_630_weights)) 





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: score cells using weighted signature
# @ atlas: Seurat object
# @ sig_weights: apply weights for scoring

SignatureScore <- function(atlas, sig_weights){
  common <- intersect(toupper(rownames(sig_weights)), toupper(rownames(atlas)))
  expr <- t(as.matrix(GetAssayData(atlas))[match(common, toupper(rownames(as.matrix(GetAssayData(atlas))))),])
  weights <- sig_weights[match(common, toupper(rownames(sig_weights))),]
  scores <- expr %*% as.matrix(weights)
  return(scores)
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: build entropy model based on ROGUE procedure
# @ atlas: Seurat object
# @ col_labels: 
# @ col_samples: 

EntropyScore <- function(atlas, col_labels, col_samples){
  expr <- atlas@assays$RNA@layers$data
  rownames(expr) <- Features(atlas@assays$RNA)
  meta <- atlas@meta.data
  results <- rogue(expr, labels = meta[[col_labels]], samples = meta[[col_samples]], platform = "UMI", span = 0.6)
  # results <- results[ , colSums(is.na(results))==0] # remove columns with any NA
  results <- results[colSums(!is.na(results)) > 0] # remove columns with all NA
  return(results)
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: construct a metadata vector to assign cells according to pseudo-bulk samples
# @ atlas: Seurat object
# @ pseudo_samples: 

PseudoBulkLabels <- function(atlas, pseudo_samples){
  num_cells <- length(Cells(atlas))
  labels <- sample.int(pseudo_samples, num_cells, replace = TRUE)
  return(labels)
}





############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: read tcsv files, e.g. for GSE120575 pre-processing

read.tcsv <- function(file, header=TRUE, sep=",", ...) {
  n <- max(count.fields(file, sep=sep), na.rm=TRUE)
  x <- readLines(file)
  
  .splitvar <- function(x, sep, n) {
    var <- unlist(strsplit(x, split=sep))
    length(var) <- n
    return(var)
  }
  
  x <- do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x <- apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out <- read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
}


############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: subroutine for scoring signatures, including CARTEx, cell states, other dysfunction
# @ atlas: Seurat object

ScoreSubroutine <- function(atlas) {
  
  # CARTEx with weights // 630 genes
  cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
  atlas@meta.data$CARTEx_630 <- Z(SignatureScore(atlas, cartex_630_weights))
  atlas@meta.data$CARTEx_630i <- integerize(atlas@meta.data$CARTEx_630)
  
  # CARTEx with weights // 200 genes
  cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
  atlas@meta.data$CARTEx_200 <- Z(SignatureScore(atlas, cartex_200_weights))
  atlas@meta.data$CARTEx_200i <- integerize(atlas@meta.data$CARTEx_200)
  
  # CARTEx with weights // 84 genes
  cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
  atlas@meta.data$CARTEx_84 <- Z(SignatureScore(atlas, cartex_84_weights))
  atlas@meta.data$CARTEx_84i <- integerize(atlas@meta.data$CARTEx_84)
  
  # load additional signatures to score
  activation.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = FALSE, row.names = 1)))
  anergy.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = FALSE, row.names = 1)))
  stemness.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = FALSE, row.names = 1)))
  senescence.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = FALSE, row.names = 1)))
  
  atlas <- AddModuleScore(atlas, features = list(activation.sig, anergy.sig, stemness.sig, senescence.sig), name="State", search = TRUE)
  
  # z score normalization
  atlas@meta.data$Activation <- scale(atlas@meta.data$State1)
  atlas@meta.data$Anergy <- scale(atlas@meta.data$State2)
  atlas@meta.data$Stemness <- scale(atlas@meta.data$State3)
  atlas@meta.data$Senescence <- scale(atlas@meta.data$State4)
  
  atlas@meta.data$Activationi <- integerize(atlas@meta.data$Activation)
  atlas@meta.data$Anergyi <- integerize(atlas@meta.data$Anergy)
  atlas@meta.data$Stemnessi <- integerize(atlas@meta.data$Stemness)
  atlas@meta.data$Senescencei <- integerize(atlas@meta.data$Senescence)
  
  atlas@meta.data$Activationi <- factor(atlas@meta.data$Activationi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
  atlas@meta.data$Anergyi <- factor(atlas@meta.data$Anergyi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
  atlas@meta.data$Stemnessi <- factor(atlas@meta.data$Stemnessi, levels = c(4,3,2,1,0,-1,-2,-3,-4))
  atlas@meta.data$Senescencei <- factor(atlas@meta.data$Senescencei, levels = c(4,3,2,1,0,-1,-2,-3,-4))
  
  atlas@meta.data$State1 <- NULL
  atlas@meta.data$State2 <- NULL
  atlas@meta.data$State3 <- NULL
  atlas@meta.data$State4 <- NULL
  
  # examine other signatures and score
  NK_like <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = FALSE, row.names = 1)))
  Wherry_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = FALSE, row.names = 1)))
  BBD_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = FALSE, row.names = 1)))
  PD1_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1)))
  TSR <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Chu_2023_Nat_Med_T_stress_response.csv", sep = ''), header = FALSE, row.names = 1)))
  Tex_Term <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Daniel_2022_Nat_Imm_TexTerm.csv", sep = ''), header = FALSE, row.names = 1)))
  Tex_KLR <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Daniel_2022_Nat_Imm_TexKLR.csv", sep = ''), header = FALSE, row.names = 1)))
  
  atlas <- AddModuleScore(atlas, features = list(NK_like, Wherry_Tex, BBD_Tex, PD1_Tex, TSR, Tex_Term, Tex_KLR), name="Signature", search = TRUE)
  
  atlas@meta.data$NKlike_Tex <- scale(atlas@meta.data$Signature1)
  atlas@meta.data$LCMV_Tex <- scale(atlas@meta.data$Signature2)
  atlas@meta.data$BBD_Tex <- scale(atlas@meta.data$Signature3)
  atlas@meta.data$PD1_Tex <- scale(atlas@meta.data$Signature4)
  atlas@meta.data$TSR <- scale(atlas@meta.data$Signature5)
  atlas@meta.data$Tex_Term <- scale(atlas@meta.data$Signature6)
  atlas@meta.data$Tex_KLR <- scale(atlas@meta.data$Signature7)
  
  atlas@meta.data$NKlike_Texi <- integerize(atlas@meta.data$NKlike_Tex)
  atlas@meta.data$LCMV_Texi <- integerize(atlas@meta.data$LCMV_Tex)
  atlas@meta.data$BBD_Texi <- integerize(atlas@meta.data$BBD_Tex)
  atlas@meta.data$PD1_Texi <- integerize(atlas@meta.data$PD1_Tex)
  atlas@meta.data$TSRi <- integerize(atlas@meta.data$TSR)
  atlas@meta.data$Tex_Termi <- integerize(atlas@meta.data$Tex_Term)
  atlas@meta.data$Tex_KLRi <- integerize(atlas@meta.data$Tex_KLR)
  
  atlas@meta.data$Signature1 <- NULL
  atlas@meta.data$Signature2 <- NULL
  atlas@meta.data$Signature3 <- NULL
  atlas@meta.data$Signature4 <- NULL
  atlas@meta.data$Signature5 <- NULL
  atlas@meta.data$Signature6 <- NULL
  atlas@meta.data$Signature7 <- NULL
  
  return(atlas)
}


############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: generate custom annotations
# @ de_genes: 
# @ select_genes: 
# @ select_genes_name: 

CustomKeyValPairsVolcanoPlot <- function(de_genes, select_genes, select_genes_name) {
  keyvals <- list()
  keyvals$shape <- ifelse(rownames(de_genes) %in% select_genes, 17, 1)
  names(keyvals$shape)[keyvals$shape == 17] <- select_genes_name
  names(keyvals$shape)[keyvals$shape == 1] <- 'Normal'
  
  keyvals$color <- ifelse(rownames(de_genes) %in% select_genes, 'darkgoldenrod', 'black')
  names(keyvals$color)[keyvals$color == 'darkgoldenrod'] <- select_genes_name
  names(keyvals$color)[keyvals$color == 'black'] <- 'Normal'
  
  keyvals$ptsize <- ifelse(rownames(de_genes) %in% select_genes, 2, 1)
  names(keyvals$ptsize)[keyvals$ptsize == 2] <- select_genes_name
  names(keyvals$ptsize)[keyvals$ptsize == 1] <- 'Normal'
  return(keyvals)
  }




############################
#^^^^^^^^^^^^^^^^^^^^^^^^^^#
# function: FeaturePlot with improved split.by functionality
# @ atlas: Seurat object
# @ features: string describing the relevant feature(s) to map
# @ split_identity: string describing the relevant metadata identity group to split by
# @ split_ids: vector of strings describing the unique identities within split_identity (custom names)
# @ color_scale: color_scale (optional)

FeaturePlotSplitBy <- function(atlas, features, split_identity, split_ids, color_scale = NULL, reduction_map = 'umap'){
  featplots <- FeaturePlot(expt.obj, reduction = reduction_map, features = features, order = FALSE, pt.size = 0.1, split.by = split_identity)
  if (is.null(color_scale)) {
    for(i in 1:length(split_ids)){featplots[[i]] <- featplots[[i]] + labs(title = split_ids[i])}
  }
  else {
    for(i in 1:length(split_ids)){featplots[[i]] <- featplots[[i]] + color_scale + labs(title = split_ids[i])}
  }
  
  return(featplots)
}

# https://github.com/satijalab/seurat/issues/5243
# https://divingintogeneticsandgenomics.com/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/





####################################################################################################
###################################### Construction functions ######################################
####################################################################################################

# These functions are used only for the /CARTEx/construction/ directory

# calculate CARTEx
cartex <- function(expression_df, metadata, weights){
  # Load files
  normgenes <- read.csv(expression_df, row.names=1, header=TRUE)
  meta <- read.csv(metadata, row.names=1, header = TRUE)
  cartex <- read.csv(weights, header = TRUE)
  
  cartex_200=cartex$PC1
  names(cartex_200) <- cartex$X
  common <- intersect(names(cartex_200), rownames(normgenes))
  expr=t(normgenes[match(common, rownames(normgenes)),])
  cartex_200=cartex_200[match(common, names(cartex_200))]
  weights = as.numeric(cartex_200)
  
  if(all(colnames(expr)==names(cartex_200)) & all(rownames(meta)==rownames(expr))){
    scores = apply(expr, 1, function(x) x %*% weights)
    scores = data.frame(cbind(scores, samples=rownames(expr), meta))
    rownames(scores) <- rownames(expr)
    scores$cartex_score <- Z(scores$scores)
    scores=scores[,-1]
    return(scores)
  }
  return(paste0("genes or samples don't match"))
}

prepare_weights <- function(cvfit, lambda){
  weights <- as.matrix(coef(cvfit, s = lambda))
  weights <- weights[apply(weights!=0, 1, all),]
  weights <- data.frame(weights[-1]) # remove intercept
  colnames(weights) <- "PC1"
  return(weights)
}

euler_fitting = function(nested_list){
  euler_fit <- euler(nested_list)
  euler_plot <- plot(euler_fit, shape = 'ellipse', quantities = list(type = 'counts'), fills = brewer.pal(length(nested_list), 'Set3'))
  euler_error <- error_plot(euler_fit)
  return(list(plot = euler_plot, error = euler_error))
}



