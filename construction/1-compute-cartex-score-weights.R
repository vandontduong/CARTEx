# minimally modified from Hima's original script

setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")

library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)

####################################################################################################
############################################# Functions ############################################
####################################################################################################


#------------
# Data
#------------
expression_data="./data/2021-04-26-cartex-experiment-log2-TPM-data.csv"
differential_analysis_results="./data/2021-04-26-cartex-gene-signature-3317-genes.csv"
samplefile="./data/CARTEx-experiment-metadata.csv"
precomputed_cluster5_genes = "./data/2021-04-30-cartex-gene-clusters.csv"

#------------------------------------------------------------------------------------------------
# Selecting 3317 CARTEx gene-list (from supervised analysis) to perform cluster analysis
#-------------------------------------------------------------------------------------------------
data <- read.csv(expression_data, header = TRUE, row.names = 1)
cart_exhaustion_genelist <- read.csv(differential_analysis_results, header = TRUE, row.names = 1)
cart_exhaustion_genelist <- cart_exhaustion_genelist$x
cart_exhaustion_genelist <- gsub("-", ".", cart_exhaustion_genelist)
missing_changed <- c("X43162", "X43346", "X43167", "X43351", "X43168")
cart_exhaustion_genelist <- c(cart_exhaustion_genelist, missing_changed)
select = data[rownames(data) %in% cart_exhaustion_genelist, ]
sampleTable <- read.csv(samplefile,header=T, row.names=1)
sampleTable=sampleTable[order(rownames(sampleTable)),]
all(rownames(sampleTable)==colnames(select)) # check if the samplenames and order match

write.csv(rownames(data), "./data/background_genes.csv", row.names=FALSE)

#------------
# Heatmap
#------------
mat <- t(scale(t(select))) # mean-center data
ann <- data.frame(Timepoint = as.factor(sampleTable$Timepoint), CAR = sampleTable$CAR)
colAnn <- HeatmapAnnotation(
  df = ann,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'), 
  annotation_legend_param = list(
    CAR = list(nrow = 2, title = 'CAR',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Timepoint = list(nrow = 3,title = 'Timepoint',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold'))
  ))
# Color heatmap
myCol <- colorRampPalette(c('blue', 'white', 'red'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

#-------------------------------------
# CLUSTER ANALYSIS
#-------------------------------------
pamClusters <- cluster::pam(mat[,-c(10:12)], k = 5)
pamClusters$clustering <- paste0('C', pamClusters$clustering)

# plot
set.seed(123)
split = c('CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','Controls','Controls','Controls','HA','HA','HA','HA','HA','HA','HA','HA','HA')
col_order = c('CD19_Donor76_Day11', 'CD19_Donor86_Day11', 'CD19_Donor90_Day11', 'CD19_Donor76_Day15', 'CD19_Donor86_Day15', 'CD19_Donor90_Day15', 'CD19_Donor76_Day21', 'CD19_Donor86_Day21', 'CD19_Donor90_Day21',
              'HA_Donor76_Day11', 'HA_Donor86_Day11', 'HA_Donor90_Day11', 'HA_Donor76_Day15', 'HA_Donor86_Day15', 'HA_Donor90_Day15', 'HA_Donor76_Day21', 'HA_Donor86_Day21', 'HA_Donor90_Day21',
              'Control_Donor76_Day0', 'Control_Donor86_Day0', 'Control_Donor90_Day0')

hmap <- Heatmap(mat, 
                name = 'Zscore',
                column_split = split, column_order = col_order,
                row_gap = unit(1.5, "mm"),
                row_split = pamClusters$clustering, 
                cluster_row_slices = FALSE,
                row_dend_reorder = TRUE, 
                cluster_rows = TRUE, 
                show_row_dend = TRUE,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'right',
                row_dend_width = unit(5,'mm'),
                column_dend_reorder = FALSE,
                column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                top_annotation = colAnn, 
                use_raster=FALSE)

# draw(hmap,  annotation_legend_side = 'left', heatmap_legend_side = 'left')
plot_heatmap_all_clusters <- as.grob(hmap)
# generate_figs(plot_heatmap_all_clusters, './plots/plot_heatmap_all_clusters', c(10,6))


# Save the clusters and gene info
dd <- data.frame(cbind(mat, cluster = pamClusters$cluster))
# write.csv(dd, "./data/2021-04-30-cartex-gene-clusters.csv")

#-----------------------------------------------------
# DERIVE CARTEx score based on cluster 5 genes
#-----------------------------------------------------
all(rownames(sampleTable)==colnames(data))
sampleTable2 = sampleTable[sampleTable$CAR %in% c("HA", "CD19"),]
sampleTable2 = droplevels(sampleTable2)
data_subset_car_of_interest = data[, colnames(data) %in% rownames(sampleTable2)]
# dd <- read.csv(precomputed_cluster5_genes, header=TRUE, row.names = 1)
# Subset cluster 5 genes from CART expression matrix
c_genes = dd[dd$cluster == "C5",]
c_genes = rownames(c_genes) # 630 genes
select_cluster_expression <- data_subset_car_of_interest[rownames(data_subset_car_of_interest) %in% c_genes,]
x=t(select_cluster_expression)
pca=prcomp(x, scale = TRUE)
scores = as.data.frame(pca$x)
# Save the weights
PC1_4=data.frame(pca$rotation[,1:4])
PC1_4$Var1=rownames(PC1_4)
weights_pc1_genes=PC1_4$PC1
save = weights_pc1_genes
names(save) <- rownames(PC1_4)
save = save[order(names(save))]
write.csv(save, file="./data/cartex-630-weights.csv")

#-----------------------------------------------------
# CARTEx200 geneset
#-----------------------------------------------------
var_genes <- apply(select_cluster_expression, 1, var)
topXmostvariable_genes=names(sort(var_genes, decreasing=TRUE))[1:200]
select_cluster_expression <- select_cluster_expression[match(topXmostvariable_genes, rownames(select_cluster_expression)),]
# Specify order of cartex based on most variable genes 
cartex_200=PC1_4[match(rownames(select_cluster_expression), rownames(PC1_4)), ]
write.csv(cartex_200, file="./data/cartex200_gene_ranks.csv")

rownames(cartex_200)
cartex_200$PC1
names(var_genes)


df_CARTEx_weights <- PC1_4[,c("PC1", "Var1")]
df_CARTEx_weights$gene_name <- df_CARTEx_weights$Var1
df_CARTEx_weights$Var1 <- NULL
df_CARTEx_weights$selected <- df_CARTEx_weights$gene_name %in% rownames(cartex_200)
# df_CARTEx_weights$gene_name <- factor(df_CARTEx_weights$gene_name, levels = df_CARTEx_weights$gene_name[order(df_CARTEx_weights$PC1)])
# reverse order
df_CARTEx_weights$gene_name <- factor(df_CARTEx_weights$gene_name, levels = df_CARTEx_weights$gene_name[order(-df_CARTEx_weights$PC1)])


plot_weights <- ggplot(df_CARTEx_weights, aes(x = PC1, y = gene_name, fill = selected)) +
  geom_bar(stat = "identity", position = "identity", width = 1) + labs(x = "Weights", y = "Genes from cluster 5") +
  theme(axis.text.y = element_blank()) + scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "orangered")) + guides(fill="none")

generate_figs(plot_weights, './plots/plot_weights', c(2,3))

ggplot(df_CARTEx_weights, aes(x = PC1, y = gene_name, fill = selected)) +
  geom_col() + labs(x = "Weights", y = "Genes from cluster 5") +
  theme(axis.text.y = element_blank()) + scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "orangered"))





#-----------------------------------------------------
# Initial PCA analysis
#-----------------------------------------------------

# Dependencies
library(DESeq2)
library(ggplot2)

# Define file paths
INPUT_FILE='./data/deseq2.RData'
OUTPUT_FILE='./plots/initial_pca.pdf'

# Load data
load(INPUT_FILE)
rld <- rlog(dds)

# Plot PCA (Figure 1A)
pcaData <- plotPCA(rld, intgroup = c("CAR", "Timepoint"), returnData = TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

Timepoint_cols <- setNames(colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4), c(0,11,15,21))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = CAR)) + 
  geom_point(size=5) + scale_color_manual(values = c('11' = '#B9D6DF', '15' = '#7390DD', '21' = '#0000CD')) +
  theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=1, size=14), 
                          legend.position="bottom", legend.text = element_text(colour="black", size = 9), 
                          legend.title = element_text(colour="black", size = 8, face="bold"), text = element_text(size=15), 
                          # panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) + 
  ylab(paste0("PC2 (", percentVar[2], "% variance)"))

# generate_figs(initial_pca, './plots/plot_initial_pca', c(5,5))


# Save figure
# ggsave(OUTPUT_FILE)
# generate_figs()



initial_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = CAR, alpha = Timepoint)) + 
  geom_point(size=5, shape = 21, color = 'black') + scale_fill_manual(values = c('CD19' = 'dodgerblue', 'HA' = 'indianred')) +
  scale_alpha_manual(values = c('11' = 0.5, '15' = 0.75, '21' = 1)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 0.2, hjust=1, size=14), 
                          legend.position="bottom", legend.text = element_text(colour="black", size = 9), 
                          legend.title = element_text(colour="black", size = 8, face="bold"), text = element_text(size=15), 
                          # panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) + 
  ylab(paste0("PC2 (", percentVar[2], "% variance)"))

generate_figs(initial_pca, './plots/plot_initial_pca', c(5,5))






