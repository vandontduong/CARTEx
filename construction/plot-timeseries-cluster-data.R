## Author: Hima Anbunathan
## Plot time series data

# ml R/4.2.0
# Rscript plot-timeseries-cluster-data.R -f ./data/Data_Figure1B.csv -c C1 -o ./plots/cluster1.png
# Rscript plot-timeseries-cluster-data.R -f ./data/Data_Figure1B.csv -c C2 -o ./plots/cluster2.png
# Rscript plot-timeseries-cluster-data.R -f ./data/Data_Figure1B.csv -c C3 -o ./plots/cluster3.png
# Rscript plot-timeseries-cluster-data.R -f ./data/Data_Figure1B.csv -c C4 -o ./plots/cluster4.png
# Rscript plot-timeseries-cluster-data.R -f ./data/Data_Figure1B.csv -c C5 -o ./plots/cluster5.png

library(optparse)
library(reshape2)
suppressMessages(library(dplyr))
library(ggplot2)

# Usage
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "Path to the gene cluster file"),
  make_option(c("-c", "--cluster_num"), type = "character", default = NULL,
              help = "Which cluster to plot. eg C1, C2 etc"), 
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output plot name"))


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

create_data_for_plot <- function(inputfile, cluster_num) {
  clusters <- read.csv(inputfile, header = TRUE, row.names = 1)
  cluster <- clusters[clusters$cluster == cluster_num,]
  print(paste0('Number of genes in cluster ',cluster_num," = ",dim(cluster)[1]))
  CAR = unlist(lapply(colnames(cluster)[-22], function(x) unlist(strsplit(x, "_"))[1]))
  Time = unlist(lapply(colnames(cluster)[-22], function(x) unlist(strsplit(x, "_"))[3]))
  pdata <- data.frame(CAR, Time)
  dat2 <- data.frame(cbind(t(cluster[,-22]), pdata))
  dat3 <- melt(dat2, id.vars = c("CAR", "Time"))
  dat3$Time <- as.factor(dat3$Time)
  dat3$CAR <- as.factor(dat3$CAR) 
  return(dat3)
  
}

if (!is.null(opt$file) && !is.null(opt$cluster_num) && !is.null(opt$output)) {
  
  dat3 <- create_data_for_plot(opt$file, opt$cluster_num)
  
  gd <- dat3 %>% 
    group_by(variable, Time, CAR) %>% 
    summarise(gene_mean_per_time = mean(value))
  
  gd2 <- gd %>% 
    group_by(CAR, Time) %>% 
    summarise(all_genes = mean(gene_mean_per_time))
  
  cluster_value = unlist(strsplit(opt$cluster_num, "C"))[2]
  
  finalplot = ggplot(gd, aes(x = Time, y = gene_mean_per_time)) + 
    geom_line(aes(group = variable), color="gray90") +
    geom_line(data = gd2, aes(x = Time, y = all_genes, group=CAR, color=CAR), linewidth=2) +
    ylab("z-score") +
    scale_color_manual(values = c("dodgerblue", "black", "indianred")) + 
    ggtitle(paste0('Cluster ',cluster_value)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=0.5, size=12), 
          legend.position="bottom", 
          legend.text = element_text(colour="black", size = 12), 
          legend.title = element_text(colour="black", size = 12, face="bold"), 
          text = element_text(size=12), 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  
  
  ggsave(opt$output, width = 3.5, height = 3)
  
} else {
  cat("\n")
  cat("Usage: Rscript plot-timeseries-cluster-data.R -f <gene-cluster_file> -c <cluster> -o <output filename>\n")
  cat("\n")
  cat("For example: Rscript plot-timeseries-cluster-data.R -f Data_Figure1B.csv -c C1 -o cluster1.pdf\n")
  cat("\n")
}