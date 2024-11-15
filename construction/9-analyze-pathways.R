# minimally modified from Hima's original script

setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")

## Generates plot for pathway enrichment analysis

library(ggplot2)
library(gprofiler2)

generate_cluster_table <- function(df, cluster_name) {
  
  cluster_subset <- df[df$cluster == cluster_name, ]
  cluster_genes <- rownames(cluster_subset)
  
  query_genes <- cluster_genes
  enrichment_result <- gost(query = query_genes,
                            organism = "hsapiens",
                            ordered_query = FALSE,
                            significant = TRUE,
                            sources = c("KEGG", "REAC"))
  
  result_df <- enrichment_result$result
  result_df <- result_df[order(result_df$p_value, decreasing = FALSE), ]
  result_df <- result_df[,colnames(result_df) %in% c("query", "p_value", "term_id","term_name")]
  result_df <- head(result_df, 5)
  result_df <- cbind(cluster=rep(cluster_name, length(result_df$query)), 
                     GO_BP=result_df$term_id, 
                     GO_pvalue=result_df$p_value, 
                     GO_term_name=result_df$term_name)
  
  return(result_df)
}


# Data
clusters <- read.csv("./data/Data_Figure_pathway_analysis.csv", header = TRUE, row.names = 1)
print(table(clusters$cluster))

data_frames <- list()
for (i in 1:5){
  cluster=paste0('C',i)
  print(cluster)
  df <- generate_cluster_table(clusters, cluster)
  data_frames[[i]] <- df
}

pathway_df <- do.call(rbind, lapply(data_frames, as.data.frame, stringsAsFactors = FALSE))
pathway_df$GO_pvalue2 <- -log10(as.numeric(as.character(pathway_df$GO_pvalue)))
pathway_df$GO_term_name <- factor(pathway_df$GO_term_name, levels = unique(as.character(pathway_df$GO_term_name[order(pathway_df$cluster)])))

# Plot
plot_pathways_analysis <- ggplot(pathway_df, aes(x=GO_term_name, y=GO_pvalue2)) + 
  geom_point(aes(color=cluster), size=2) + 
  coord_flip() + 
  xlab("") +
  ylab("Enrichment (-log10 P)") +  
  theme_classic() + 
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8), 
        axis.title.y = element_text(size=7),
        legend.position="bottom", 
        legend.text = element_text(colour="black", size = 9), 
        legend.title = element_text(colour="black", size = 8, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3), 
        panel.grid.major = element_line(linewidth=0.3)
  ) + 
  facet_grid(~cluster) + 
  scale_color_manual(values=c("lightsalmon", "lightpink2", "lightpink3", "plum3", "mediumorchid1"), name="")


generate_figs(plot_pathways_analysis, "./plots/plot_pathways_analysis", c(6,5))




