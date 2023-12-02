# CARTEx utilities - custom functions
# Vandon Duong

# build package
# http://rcs.bu.edu/examples/r/tutorials/BuildingPackages/r-project-setup.html

# additional functions
# https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html

####################################################################################################
######################################### General functions ########################################
####################################################################################################

#####
# Z score

Z <- function(s){
  s <- as.numeric(s)
  z <- (s - mean(s))/sd(s)
  return(z)
}

#####
# save figures in jpeg and pdf formats

generate_figs <- function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

#####
# round scores to nearest integer; cut-off at +/-5

integerize <- function(score){
  score_mod <- round(score)
  score_mod[score_mod < -4] <- -5
  score_mod[score_mod > 4] <- 5
  return (score_mod)
}

#####
# Check metadata levels() and return "NULL" if levels() does not exist
# need to update levels() for anything that results in "not NULL"

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


#####
# DimPlot() customized to highlight cells by identity group

DimPlotHighlightIdents <- function(atlas, identity, reduction_map, highlight_color, pt_size, ncols){
  plot.list <- list()
  for (i in sapply(unique(x = atlas[[deparse(substitute(identity))]]), levels)) {
    plot.list[[i]] <- DimPlot(
      object = atlas, reduction = reduction_map, raster = FALSE, cols.highlight = highlight_color, pt.size = pt_size, sizes.highlight = pt_size,
      cells.highlight = Cells(atlas[, atlas[[deparse(substitute(identity))]] == i])
    ) + NoLegend() + ggtitle(i)
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

#####
# Bar plot stacked and split according to two different identity groups

BarPlotSplit <- function(atlas, x_identity, y_identity){
  md <- as.data.table(atlas@meta.data)[, .N, by = c(x_identity, y_identity)]
  md$percent <- round(100*md$N / sum(md$N), digits = 1)
  barplot <- ggplot(md, aes_string(x = y_identity, y = 'N', fill = x_identity)) + geom_col(position = "fill")
  return(barplot)
}

#####
# read tcsv files, e.g. for GSE120575 pre-processing

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



