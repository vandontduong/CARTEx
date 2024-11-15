


setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")


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


### canonical exhaustion markers

intersect(rownames(data), c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1"))

expression_canonical_exhaustion_markers <- data[c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1"),]
expression_canonical_exhaustion_markers <- expression_canonical_exhaustion_markers[,c("Control_Donor76_Day0", "Control_Donor86_Day0", "Control_Donor90_Day0", "CD19_Donor76_Day11", "CD19_Donor86_Day11", "CD19_Donor90_Day11",  "CD19_Donor76_Day15", "CD19_Donor86_Day15", "CD19_Donor90_Day15",  "CD19_Donor76_Day21", "CD19_Donor86_Day21", "CD19_Donor90_Day21", "HA_Donor76_Day11", "HA_Donor86_Day11", "HA_Donor90_Day11", "HA_Donor76_Day15", "HA_Donor86_Day15", "HA_Donor90_Day15", "HA_Donor76_Day21", "HA_Donor86_Day21", "HA_Donor90_Day21")]

CAR <- c('Control','Control','Control','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','HA','HA','HA','HA','HA','HA','HA','HA','HA')
CAR <- factor(CAR, levels = c("CD19", "HA", "Control"))
Days <- c(0, 0, 0, 11, 11, 11, 15, 15, 15, 21, 21, 21, 11, 11, 11, 15, 15, 15, 21, 21, 21)
Days <- as.character(Days)
Donor <- c("Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90")

library(tidyr)

df_canonical_exhaustion_markers <- data.frame(t(expression_canonical_exhaustion_markers), CAR, Days, Donor)

melted_df <- reshape2::melt(df_canonical_exhaustion_markers, id.vars = c("Days", "CAR", "Donor"))
melted_df <- melted_df %>% group_by(CAR, Days, variable) %>% summarise_at(vars(value), list(name = mean, sd))

custom_colors <- c("CD19" = "dodgerblue", "HA" = "indianred", "Control" = "darkgoldenrod")

plot_canonical_exhaustion_markers <- ggplot(melted_df, aes(x = Days, y = name, fill = CAR)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~variable, scales = "free_y") +
  geom_errorbar(aes(ymin = name - fn1, ymax = name + fn1), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Days", y = "Gene Expression", fill = "CAR") + ylim(0,10) + 
  scale_fill_manual(values = custom_colors) + theme_classic()

# log2-transformed TPM values
generate_figs(plot_canonical_exhaustion_markers, './plots/plot_canonical_exhaustion_markers', c(4,2))




# check literature-identified markers
expression_literature_exhaustion_markers <- data[c("RGS16", "METRNL", "EBI3", "CXCL10", "SNX9", "TNFRSF9"),]
expression_literature_exhaustion_markers <- expression_literature_exhaustion_markers[,c("Control_Donor76_Day0", "Control_Donor86_Day0", "Control_Donor90_Day0", "CD19_Donor76_Day11", "CD19_Donor86_Day11", "CD19_Donor90_Day11",  "CD19_Donor76_Day15", "CD19_Donor86_Day15", "CD19_Donor90_Day15",  "CD19_Donor76_Day21", "CD19_Donor86_Day21", "CD19_Donor90_Day21", "HA_Donor76_Day11", "HA_Donor86_Day11", "HA_Donor90_Day11", "HA_Donor76_Day15", "HA_Donor86_Day15", "HA_Donor90_Day15", "HA_Donor76_Day21", "HA_Donor86_Day21", "HA_Donor90_Day21")]

df_literature_exhaustion_markers <- data.frame(t(expression_literature_exhaustion_markers), CAR, Days, Donor)

melted_df <- reshape2::melt(df_literature_exhaustion_markers, id.vars = c("Days", "CAR", "Donor"))
melted_df <- melted_df %>% group_by(CAR, Days, variable) %>% summarise_at(vars(value), list(name = mean, sd))

custom_colors <- c("CD19" = "dodgerblue", "HA" = "indianred", "Control" = "darkgoldenrod")

plot_literature_exhaustion_markers <- ggplot(melted_df, aes(x = Days, y = name, fill = CAR)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~variable, scales = "free_y") +
  geom_errorbar(aes(ymin = name - fn1, ymax = name + fn1), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Days", y = "Gene Expression", fill = "CAR") + ylim(0,10) + 
  scale_fill_manual(values = custom_colors) + theme_classic()






#------------
# Heatmap
#------------
sampleTable2 <- sampleTable[c("Control_Donor76_Day0", "Control_Donor86_Day0", "Control_Donor90_Day0", "CD19_Donor76_Day11", "CD19_Donor86_Day11", "CD19_Donor90_Day11",  "CD19_Donor76_Day15", "CD19_Donor86_Day15", "CD19_Donor90_Day15",  "CD19_Donor76_Day21", "CD19_Donor86_Day21", "CD19_Donor90_Day21", "HA_Donor76_Day11", "HA_Donor86_Day11", "HA_Donor90_Day11", "HA_Donor76_Day15", "HA_Donor86_Day15", "HA_Donor90_Day15", "HA_Donor76_Day21", "HA_Donor86_Day21", "HA_Donor90_Day21"),]

mat <- t(scale(t(select))) # mean-center data
ann <- data.frame(Timepoint = as.factor(sampleTable2$Timepoint), CAR = sampleTable2$CAR)

# Color heatmap
myCol <- colorRampPalette(c('blue', 'white', 'red'))(100)
myBreaks <- seq(-3, 3, length.out = 100)


#-------------------------------------
# CLUSTER ANALYSIS
#-------------------------------------
pamClusters <- cluster::pam(mat[,-c(10:12)], k = 5)
pamClusters$clustering <- paste0('C', pamClusters$clustering)
table(pamClusters$clustering)


cartex_630_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-630-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)
cartex_84_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-84-weights.csv", sep = ''), header = TRUE, row.names = 1)
NK_like <- rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = TRUE, row.names = 1))
Wherry_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = TRUE, row.names = 1))
BBD_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = TRUE, row.names = 1))
PD1_Tex <- rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = TRUE, row.names = 1))
activation.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = TRUE, row.names = 1))
anergy.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = TRUE, row.names = 1))
stemness.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = TRUE, row.names = 1))
senescence.sig <- rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = TRUE, row.names = 1))


STAT3_GOF_sig_stim <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_stim.csv"), header = TRUE, row.names = 1))
STAT3_GOF_sig_rest <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_rest.csv"), header = TRUE, row.names = 1))
STAT3_GOF_sig_combined <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_combined.csv"), header = TRUE, row.names = 1))




# reorder columns
mat <- mat[, c("Control_Donor76_Day0", "Control_Donor86_Day0", "Control_Donor90_Day0", "CD19_Donor76_Day11", "CD19_Donor86_Day11", "CD19_Donor90_Day11",  "CD19_Donor76_Day15", "CD19_Donor86_Day15", "CD19_Donor90_Day15",  "CD19_Donor76_Day21", "CD19_Donor86_Day21", "CD19_Donor90_Day21", "HA_Donor76_Day11", "HA_Donor86_Day11", "HA_Donor90_Day11", "HA_Donor76_Day15", "HA_Donor86_Day15", "HA_Donor90_Day15", "HA_Donor76_Day21", "HA_Donor86_Day21", "HA_Donor90_Day21")]

rownames(mat)
index_CARTEx_630 <- rownames(mat) %in% rownames(cartex_630_weights)
index_CARTEx_200 <- rownames(mat) %in% rownames(cartex_200_weights)
index_CARTEx_84 <- rownames(mat) %in% rownames(cartex_84_weights)
index_NKlike <- rownames(mat) %in% NK_like
index_Wherry <- rownames(mat) %in% Wherry_Tex
index_BBD <- rownames(mat) %in% BBD_Tex
index_PD1 <- rownames(mat) %in% PD1_Tex
index_activation <- rownames(mat) %in% activation.sig
index_anergy <- rownames(mat) %in% anergy.sig
index_stemness <- rownames(mat) %in% stemness.sig
index_senescence <- rownames(mat) %in% senescence.sig
index_STAT3_GOF_sig_stim <- rownames(mat) %in% STAT3_GOF_sig_stim
index_STAT3_GOF_sig_rest <- rownames(mat) %in% STAT3_GOF_sig_rest
index_STAT3_GOF_sig_combined <- rownames(mat) %in% STAT3_GOF_sig_combined

intersect_CARTEx_630 <- rownames(mat)[index_CARTEx_630]
intersect_CARTEx_200 <- rownames(mat)[index_CARTEx_200]
intersect_CARTEx_84 <- rownames(mat)[index_CARTEx_84]
intersect_NKlike <- rownames(mat)[index_NKlike]
intersect_Wherry <- rownames(mat)[index_Wherry]
intersect_BBD <- rownames(mat)[index_BBD]
intersect_PD1 <- rownames(mat)[index_PD1]
intersect_activation <- rownames(mat)[index_activation]
intersect_anergy <- rownames(mat)[index_anergy]
intersect_stemness <- rownames(mat)[index_stemness]
intersect_senescence <- rownames(mat)[index_senescence]
intersect_STAT3_GOF_sig_stim <- rownames(mat)[index_STAT3_GOF_sig_stim]
intersect_STAT3_GOF_sig_rest <- rownames(mat)[index_STAT3_GOF_sig_rest]
intersect_STAT3_GOF_sig_combined <- rownames(mat)[index_STAT3_GOF_sig_combined]


intersect(intersect_CARTEx_630, intersect_NKlike)
intersect(intersect_Wherry, intersect_NKlike)

# as.data.frame(rownames(mat))

mat_CARTEx_630 <- mat[index_CARTEx_630,]
mat_CARTEx_200 <- mat[index_CARTEx_200,]
mat_CARTEx_84 <- mat[index_CARTEx_84,]
mat_NKlike <- mat[index_NKlike,]
mat_Wherry <- mat[intersect_Wherry,]
mat_BBD <- mat[intersect_BBD,]
mat_PD1 <- mat[intersect_PD1,]
mat_activation <- mat[intersect_activation,]
mat_anergy <- mat[intersect_anergy,]
mat_senescence <- mat[intersect_senescence,]
mat_stemness <- mat[intersect_stemness,]
mat_STAT3_GOF_sig_stim <- mat[intersect_STAT3_GOF_sig_stim,]
mat_STAT3_GOF_sig_rest <- mat[intersect_STAT3_GOF_sig_rest,]
mat_STAT3_GOF_sig_combined <- mat[intersect_STAT3_GOF_sig_combined,]

# plot
set.seed(123)

split = c('Control','Control','Control','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','CD19','HA','HA','HA','HA','HA','HA','HA','HA','HA')
col_order = c('Control_Donor76_Day0', 'Control_Donor86_Day0', 'Control_Donor90_Day0',
              'CD19_Donor76_Day11', 'CD19_Donor86_Day11', 'CD19_Donor90_Day11', 'CD19_Donor76_Day15', 'CD19_Donor86_Day15', 'CD19_Donor90_Day15', 'CD19_Donor76_Day21', 'CD19_Donor86_Day21', 'CD19_Donor90_Day21',
              'HA_Donor76_Day11', 'HA_Donor86_Day11', 'HA_Donor90_Day11', 'HA_Donor76_Day15', 'HA_Donor86_Day15', 'HA_Donor90_Day15', 'HA_Donor76_Day21', 'HA_Donor86_Day21', 'HA_Donor90_Day21')


colAnn <- HeatmapAnnotation(
  df = ann,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'), 
  annotation_legend_param = list(
    CAR = list(nrow = 2, title = 'CAR',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Timepoint = list(nrow = 3,title = 'Timepoint',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold'))),
  col = list(CAR = c("Control" = "darkgoldenrod", "CD19" = "dodgerblue", "HA" = "indianred"),
             Timepoint = setNames(colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4), c(0,11,15,21)))
  )


# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# select_genes <- c('FOS', 'CSF1', 'TCF7', 'BTLA', 'ID3', 'NFATC1', 'KLRG1', 'CD160', 'NFKB1', 'TOX', 'GZMA', 'BATF', 'EOMES', 'PDCD1', 'ZEB2', 'CXCR5', 'JUN', 'IFNG', 'RUNX3', 'NR4A1', 'TNFRSF9', 'LAG3', 'GZMB', 'ENTPD1', 'IL21R', 'CTLA4', 'TIGIT')
# select_genes <- c('CSF1', 'CD44', 'IL2', 'IL3', 'STAT1', 'STAT3', 'STAT4', 'KLRG1', 'TOX', 'TCF7', "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1", 'FOS', 'BTAF', 'EOMES', 'JUN', 'CD160', 'BTLA', 'GZMA', 'GZMB', 'NFATC1', 'TBX21', 'RUNX3', 'IL21R', 'TNFRSF9', 'INFG', 'ZEB2', 'ICOSLG', 'NFKB1', 'NFKB2', 'REL', 'TNFRSF11A', 'CD40LG', 'SELL', 'SLAMF6', 'CD69')
select_genes <- c('CSF1', 'CD44', 'IL2', 'IL3', 'STAT3', 'KLRG1', 'TOX', 'TCF7', "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1", 'FOS', 'BTAF', 'EOMES', 'JUN', 'CD160', 'BTLA', 'GZMA', 'GZMB', 'NFATC1', 'TBX21', 'RUNX3', 'IL21R', 'TNFRSF9', 'INFG', 'ZEB2', 'ICOSLG', 'REL', 'TNFRSF11A', 'CD40LG', 'SELL', 'SLAMF6', 'CD69')
# select_genes <- c('MYC', 'BRD4', 'BRD3', 'BRD2', 'BRDT', 'CDK9', 'CCNT1', 'HIF1A', 'HIF2A', 'HIF1B', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'CD137', 'TNFRSF11A', 'TNFRSF11B', 'MAP3K14', 'ICOSLG', 'NFKBIA', 'NFKBIB', 'NFKBIE', 'EP300', 'CREBBP', 'TET2')

# examine NFKB pathways
# https://www.sciencedirect.com/science/article/pii/S1074761323002649
# Pichler, Andrea C., et al. "TCR-independent CD137 (4-1BB) signaling promotes CD8+-exhausted T cell proliferation and terminal differentiation." Immunity 56.7 (2023): 1631-1648.
# select_genes <-c('NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'CD137', 'TNFRSF11A', 'TNFRSF11B', 'MAP3K14', 'ICOSLG', 'NFKBIA', 'NFKBIB', 'NFKBIE', 'TOX', 'NFATC1', 'TRAF1', 'TRAF2', 'TRAF3', 'TNFRSF9')

# https://link.springer.com/article/10.1007/s13277-016-5331-4
# Slemc, Lucija, and Tanja Kunej. "Transcription factor HIF1A: downstream targets, associated pathways, polymorphic hypoxia response element (HRE) sites, and initiative for standardization of reporting in scientific literature." Tumor Biology 37 (2016): 14851-14861.
# HIF1A target genes
# select_genes <- c('HIF1A', 'VEGFA', 'ACAN', 'ADM', 'ALDOA', 'ALDOC', 'ANGPT2', 'ANKRD37', 'BHLHE40', 'BNIP3', 'CA9', 'COL2A1', 'ENO1', 'GAPDH', 'LDHA', 'LRP1', 'PFKFB4', 'PFKL', 'PGK1', 'PKM', 'SLC2A1', 'SLC2A3')

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4342852/
# https://en.wikipedia.org/wiki/Aryl_hydrocarbon_receptor_nuclear_translocator
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7182078/
# genes interacting with HIF1B or AhR
# select_genes <- c('HIF1B', 'ARNT', 'ARNT1', 'ARNT2', 'TRIP230', 'TACC3', 'AIP', 'EPAS1', 'NCOA2', 'SIM1', 'SIM2', 'AHR', 'IDO1', 'TDO2', 'HSP90AA1', 'HSP90AA2', 'HSP90AB1', 'HSP90B1', 'TRAP1')


index_select_genes <- which(rownames(mat) %in% select_genes)
select_genes_ordered <- rownames(mat)[index_select_genes]

rightAnn <- rowAnnotation(CAR = anno_mark(at = index_select_genes, labels = select_genes_ordered))


gene_cluster_pairs <- data.frame(rownames(pamClusters$data),pamClusters$clustering)
colnames(gene_cluster_pairs) <- c('gene', 'cluster')

select_gene_cluster_pairs <- gene_cluster_pairs %>% filter_all(any_vars(. %in% select_genes))
select_gene_cluster_pairs[order(select_gene_cluster_pairs$cluster),]

# check number of genes in each cluster
table(gene_cluster_pairs$cluster)

# search a specific gene for its corresponding cluster
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "STAT3"] # C5
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "SLAMF6"] # C1
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "CD69"] # C1

gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "SOCS3"] # C1 and inhibitor of STAT3 pathway
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "PIAS3"] # not present and inhibitor of STAT3 pathway

gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "STAT1"] # C1
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "IFNG"] # C5
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "TBX21"] # C3 T-Bet
gene_cluster_pairs$cluster[gene_cluster_pairs$gene == "SOCS3"] # C1


# list genes in a specific cluster
gene_cluster_pairs$gene[gene_cluster_pairs$cluster == "C2"]
write.csv(gene_cluster_pairs$gene[gene_cluster_pairs$cluster == "C2"], file = './data/cartex-cluster-2.csv', row.names = F)
# cartex_C2 <- rownames(read.csv('./data/cartex-cluster-2.csv', header = TRUE, row.names = 1))

write.csv(gene_cluster_pairs$gene[gene_cluster_pairs$cluster == "C3"], file = './data/cartex-cluster-3.csv', row.names = F)

hmap <- Heatmap(mat, 
                name = 'Zscore',
                column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
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
                # column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                show_column_names = FALSE,
                top_annotation = colAnn,
                right_annotation = rightAnn,
                use_raster=FALSE)

# cluster 2 is enriched with adhesion receptor CD44 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3038050/), enriched in activated T cells and T cell migration
# cluster 3 is enriched with stimulatory ligands IL2, IL3, ICOSLG; C3 also has T-bet (TBX21) and SOX2
# SOX2 and STAT3 are elevated

# dev.off()

# draw(hmap,  annotation_legend_side = 'left', heatmap_legend_side = 'left')
plot_heatmap_all_clusters <- as.grob(hmap)
generate_figs(plot_heatmap_all_clusters, './plots/plot_heatmap_all_clusters', c(6,7)) # due to anno_mark, which is not scalable for some reason, export manually

# pdf(plot_heatmap_all_clusters)
# png(filename = "./plots/plot_heatmap_all_clusters.png", width = 480, height = 480, units = "px", pointsize = 12)

# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html


# CREATE HEATMAP WITHOUT GENE ANNOTATIONS AND LEGENDS

colAnn_noAnno <- HeatmapAnnotation(
  df = ann,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'), 
  annotation_legend_param = list(
    CAR = list(nrow = 2, title = 'CAR',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Timepoint = list(nrow = 3,title = 'Timepoint',title_position = 'topcenter',legend_direction = 'vertical',title_gp = gpar(fontsize = 12, fontface = 'bold'),labels_gp = gpar(fontsize = 12, fontface = 'bold'))),
  col = list(CAR = c("Control" = "darkgoldenrod", "CD19" = "dodgerblue", "HA" = "indianred"),
             Timepoint = setNames(colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4), c(0,11,15,21))),
  show_legend = c(FALSE, FALSE, FALSE)
)

hmap_noAnno <- Heatmap(mat, 
                name = 'Zscore',
                column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
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
                # column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
                show_column_names = FALSE,
                top_annotation = colAnn_noAnno,
                use_raster=FALSE)

plot_heatmap_all_clusters_noAnno <- as.grob(hmap_noAnno)
generate_figs(plot_heatmap_all_clusters_noAnno, './plots/plot_heatmap_all_clusters_noAnno', c(8,7))



## 







hmap_CARTEx_630 <- Heatmap(mat_CARTEx_630, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_CARTEx_630], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_630), './plots/plot_heatmap_CARTEx_630', c(6,6))

hmap_CARTEx_200 <- Heatmap(mat_CARTEx_200, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_CARTEx_200], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_200), './plots/plot_heatmap_CARTEx_200', c(6,6))

hmap_CARTEx_84 <- Heatmap(mat_CARTEx_84, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_CARTEx_84], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_CARTEx_84), './plots/plot_heatmap_CARTEx_84', c(6,6))

hmap_NKlike <- Heatmap(mat_NKlike, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_NKlike], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_NKlike), './plots/plot_heatmap_NKlike', c(6,6))

hmap_Wherry <- Heatmap(mat_Wherry, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_Wherry], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_Wherry), './plots/plot_heatmap_Wherry', c(6,6))

hmap_BBD <- Heatmap(mat_BBD, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_BBD], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_BBD), './plots/plot_heatmap_BBD', c(6,6))

hmap_PD1 <- Heatmap(mat_PD1, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                    row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_PD1], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                    cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                    row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                    top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_PD1), './plots/plot_heatmap_PD1', c(6,6))

hmap_activation <- Heatmap(mat_activation, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                    row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_activation], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                    cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                    row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                    top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_activation), './plots/plot_heatmap_activation', c(6,6))

hmap_anergy <- Heatmap(mat_anergy, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_anergy], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_anergy), './plots/plot_heatmap_anergy', c(6,6))

hmap_senescence <- Heatmap(mat_senescence, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_senescence], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_senescence), './plots/plot_heatmap_senescence', c(6,6))

hmap_stemness <- Heatmap(mat_stemness, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_stemness], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_stemness), './plots/plot_heatmap_stemness', c(6,6))






STAT3_GOF_sig_stim <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_stim.csv"), header = TRUE, row.names = 1))
STAT3_GOF_sig_rest <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_rest.csv"), header = TRUE, row.names = 1))
STAT3_GOF_sig_combined <- rownames(read.csv(paste0(PATH_EXPERIMENTS, "GSE207935/data/GSE207935_STAT3_GOF_combined.csv"), header = TRUE, row.names = 1))


length(STAT3_GOF_sig_stim) - length(intersect(STAT3_GOF_sig_stim, rownames(mat)))
length(STAT3_GOF_sig_rest) - length(intersect(STAT3_GOF_sig_rest, rownames(mat)))
length(STAT3_GOF_sig_combined) - length(intersect(STAT3_GOF_sig_combined, rownames(mat)))



hmap_STAT3_GOF_sig_stim <- Heatmap(mat_STAT3_GOF_sig_stim, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                         row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_STAT3_GOF_sig_stim], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                         cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                         row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                         top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_STAT3_GOF_sig_stim), './plots/plot_heatmap_STAT3_GOF_sig_stim', c(6,6))

hmap_STAT3_GOF_sig_rest <- Heatmap(mat_STAT3_GOF_sig_rest, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                                   row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_STAT3_GOF_sig_rest], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                                   cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                                   row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                                   top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_STAT3_GOF_sig_rest), './plots/plot_heatmap_STAT3_GOF_sig_rest', c(6,6))

hmap_STAT3_GOF_sig_combined <- Heatmap(mat_STAT3_GOF_sig_combined, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                                   row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_STAT3_GOF_sig_combined], cluster_row_slices = FALSE, row_dend_reorder = TRUE, 
                                   cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                                   row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                                   top_annotation = colAnn, use_raster=FALSE)
generate_figs(as.grob(hmap_STAT3_GOF_sig_combined), './plots/plot_heatmap_STAT3_GOF_sig_combined', c(6,6))







CAR <- split
Days <- c(0, 0, 0, 11, 15, 21, 11, 15, 21, 11, 15, 21, 11, 15, 21, 11, 15, 21, 11, 15, 21)
Days <- as.character(Days)
Donor <- c("Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90", "Donor76", "Donor86", "Donor90")

mean_CARTEx_630 <- colMeans(mat_CARTEx_630)
mean_CARTEx_200 <- colMeans(mat_CARTEx_200)
mean_CARTEx_84 <- colMeans(mat_CARTEx_84)
mean_NKlike <- colMeans(mat_NKlike)
mean_Wherry <- colMeans(mat_Wherry)
mean_BBD <- colMeans(mat_BBD)
mean_PD1 <- colMeans(mat_PD1)
mean_activation <- colMeans(mat_activation)
mean_anergy <- colMeans(mat_anergy)
mean_senescence <- colMeans(mat_senescence)
mean_stemness <- colMeans(mat_stemness)

mean_dataset <- data.frame(mean_CARTEx_630, mean_CARTEx_200, mean_CARTEx_84, mean_NKlike, mean_Wherry, mean_BBD, mean_PD1,
                           mean_activation, mean_anergy, mean_senescence, mean_stemness, CAR, Days, Donor)



library(tidyr)
df_long <- gather(mean_dataset, key = "measurement", value = "value", -CAR, -Days, -Donor)

ggplot(df_long, aes(x = Days, y = value, group = paste(CAR, Donor), color = CAR)) +
  geom_line() +
  facet_wrap(~measurement, scales = "free_y") +
  labs(title = "Line Plot of Measurements by Category and Type",
       x = "Category",
       y = "Measurement")

mean_dataset_avgbydonor <- mean_dataset %>%
  group_by(CAR, Days, Donor) %>%
  summarize(across(starts_with("mean"), mean))


df_long_avgbydonor <- gather(mean_dataset_avgbydonor, key = "measurement", value = "value", -CAR, -Days, -Donor)
df_long_avgbydonor <- subset(df_long_avgbydonor, CAR != "Control")
df_long_avgbydonor <- subset(df_long_avgbydonor, measurement != "mean_CARTEx_84")
df_long_avgbydonor$measurement <- factor(df_long_avgbydonor$measurement, levels = c("mean_CARTEx_630", "mean_CARTEx_200", "mean_Wherry", "mean_NKlike", "mean_BBD", "mean_PD1", "mean_activation", "mean_anergy",  "mean_senescence", "mean_stemness"))

measurement_names <- c(mean_activation = "Activation", mean_anergy = "Anergy", mean_senescence = "Senescence", mean_stemness = "Stemness", mean_CARTEx_630 = "C5", mean_CARTEx_200 = "CARTEx",  mean_Wherry = "LCMV",  mean_NKlike = "NK-like", mean_BBD = "BBD",  mean_PD1 = "PD1")

plot_mean_zscore_by_signature <- ggplot(df_long_avgbydonor, aes(x = Days, y = value, group = CAR, color = CAR)) +
  geom_line() + geom_point() +
  facet_wrap(~measurement, scales = "free_y", labeller = labeller(measurement = measurement_names), ncol = 5) +
  labs(x = "Days", y = "Mean z-score") +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  scale_color_manual(values = c("CD19" = "dodgerblue", "HA" = "indianred")) +
  theme_classic() + theme(legend.position = "none")

generate_figs(plot_mean_zscore_by_signature, './plots/plot_mean_zscore_by_signature', c(5,2))


intersect(rownames(mat), c('HNF1A', 'HNF1B'))
mat[rownames(mat) %in% c('HNF1A', 'HNF1B'),]

HNF1B_score <- mat[rownames(mat) %in% c('HNF1B'),]
CAR <- split
Days <- c(11, 15, 21, 11, 15, 21, 11, 15, 21, 0, 0, 0, 11, 15, 21, 11, 15, 21, 11, 15, 21)
Days <- as.character(Days)

HNF1B_dataset <- data.frame(HNF1B_score, CAR, Days)

ggplot(data=HNF1B_dataset, aes(x=CAR, y=HNF1B_score)) + geom_boxplot() + geom_point(aes(color=Days)) + xlab("CAR") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','CD19')), label = "p.signif", label.y = 1.8) +
  theme_bw()

colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4)








#### 

qk_anno <- function(mat_select, sele_genes){
  idx_sele_genes <- which(rownames(mat_select) %in% sele_genes)
  ord_sele_genes <- rownames(mat_select)[idx_sele_genes]
  right_anno <- rowAnnotation(CAR = anno_mark(at = idx_sele_genes, labels = ord_sele_genes))
  return(right_anno)
}

qk_assemble <- function(gene_ID, cluster_ID){
  df <- data.frame(gene_ID, cluster_ID)
  colnames(df) <- c('Gene', 'Cluster')
  return(df)
}




rownames(mat_CARTEx_630)
setdiff(rownames(mat_CARTEx_630), rownames(mat_CARTEx_200))
# right_anno <- qk_anno(mat_CARTEx_630, c('ZEB2', 'JUN', 'REL', 'RUNX3', 'TNFRSF9', 'LAG3', 'GZMB', 'ENTPD1', 'STAT3', 'IL21R', 'CTLA4', 'TIGIT', 'SOX4', 'IFNG', 'CD70', 'DUSP4', 'LAIR2', 'GNLY', 'CFLAR', 'ATP9A', 'METRNL'))
right_anno <- qk_anno(mat_CARTEx_630, c('JUN', 'REL', 'RUNX3', 'LAG3', 'ENTPD1', 'STAT3', 'IL21R', 'CTLA4', 'CD70', 'CFLAR', 'IL26', 'RGS2', 'IL1A', 'GZMK', 'HOXB6', 'CDKN1A', 'IL36G', 'RAB30', 'CD276', 'CDK6', 'FASLG', 'NINJ1', 'ICAM1', 'RUNX2'))
hmap_CARTEx_630 <- Heatmap(mat_CARTEx_630, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_CARTEx_630], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_CARTEx_630

# novel genes: XKRX, ATP9A, METRNL

qk_assemble(rownames(mat_CARTEx_200), pamClusters$clustering[index_CARTEx_200])
right_anno <- qk_anno(mat_CARTEx_200, c('KLRC1', 'KLRD1', 'ZEB2', 'TNFRSF9', 'GZMB', 'TIGIT', 'SOX4', 'IFNG', 'KIR2DL4', 'XKRX', 'LAIR2', 'GNLY', 'CFLAR', 'ATP9A', 'METRNL', 'CXCL10', 'CCL4', 'CCL3', 'KLRC2', 'FOSB', 'KIR3DX1', 'MAGI1', 'TWIST1', 'RGS16'))
hmap_CARTEx_200 <- Heatmap(mat_CARTEx_200, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_CARTEx_200], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_CARTEx_200
sum(rownames(mat_CARTEx_200) == 'CFLAR')


rownames(mat_Wherry)
pamClusters$clustering[index_Wherry]
right_anno <- qk_anno(mat_Wherry, c('CD160', 'RGS16', 'EOMES', 'NFATC1', 'PDCD1', 'CCL3', 'CCL4', 'CXCL10', 'SPP1', 'CTLA4', 'CCRL2', 'LAG3', 'TNFRSF9', 'TNFRSF1A', 'ENTPD1', 'CASP4', 'LILRB4', 'PAWR', 'TCEAL9', 'TCEA2', 'SPOCK2'))
hmap_Wherry <- Heatmap(mat_Wherry, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_Wherry], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_Wherry


rownames(mat_NKlike)
right_anno <- qk_anno(mat_NKlike, c('PLS3', 'CDK6', 'ID3', 'KLRD1', 'KLRC1', 'KLRC2', 'RGS16', 'CSF1', 'TNFRSF9', 'TNFRSF18', 'NDFIP2', 'SRGAP3', 'KLRB1', 'CCL3', 'CCL4', 'SQLE', 'LYST', 'GZMB', 'IL2RA', 'SOX4', 'GNLY'))
hmap_NKlike <- Heatmap(mat_NKlike, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_NKlike], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_NKlike


rownames(mat_BBD)
right_anno <- qk_anno(mat_BBD, c('KLRC3', 'KLRD1', 'KIT', 'IL9R', 'KLF6', 'JUN', 'CCR8', 'IGFLR1', 'ATP8B4', 'SERPINE1', 'CD69', 'GZMK', 'GNLY', 'CXCR6', 'ITGA1', 'TXNIP', 'CD27', 'CD8B', 'KLRK1', 'NCR3', 'TSPAN32', 'FYB1', 'YPEL3', 'GADD45B'))
hmap_BBD <- Heatmap(mat_BBD, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_BBD], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_BBD


rownames(mat_PD1)
right_anno <- qk_anno(mat_PD1, c('RGS1', 'LAG3', 'TNFRSF9', 'SYNGR3', 'CD200', 'ATP8B4', 'TOX', 'BTLA', 'TIGIT', 'KLRB1', 'TNFSF4', 'CD80', 'CCR2', 'CD70', 'RGS13', 'CTLA4', 'VCAM1', 'CXCL13', 'BATF', 'IRF5', 'ZBED2', 'GLDC', 'L1CAM', 'CARD16', 'KCNN4'))
hmap_PD1 <- Heatmap(mat_PD1, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                    row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_PD1], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                    cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                    row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                    top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_PD1

rownames(mat_activation)
right_anno <- qk_anno(mat_activation, c('CD80', 'JUN', 'FOS', 'PIK3CG', 'PIK3CB', 'PIK3R1', 'ITPR1', 'MAPK8', 'NFKBIA', 'SOS1'))
hmap_activation <- Heatmap(mat_activation, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                    row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_activation], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                    cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                    row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                    top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_activation

rownames(mat_anergy)
right_anno <- qk_anno(mat_anergy, c('LAG3', 'CCL1', 'CSF1', 'NOTCH1', 'CCL3', 'NFATC1', 'TNFSF9', 'TNFSF11', 'NOCT', 'CASP4', 'DDR1', 'FURIN', 'DUSP6', 'CD40LG', 'ZNF629', 'EGR2', 'CDC14A', 'GADD45B', 'ING4', 'HLF', 'HSPA1A', 'JUP'))
hmap_anergy <- Heatmap(mat_anergy, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_anergy], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_anergy

qk_assemble(rownames(mat_senescence), pamClusters$clustering[index_senescence])
right_anno <- qk_anno(mat_senescence, c('CDKN1A', 'CYP1B1', 'STAT1', 'IFNG', 'CD44', 'RHOB', 'CDKN2B', 'IGFBP7', 'FILIP1L', 'ALDH1A3', 'IRF7', 'HSPA2', 'TGFB1I1', 'IGFBP4', 'SERPINE1', 'MAP2K3', 'RAB13', 'TFAP2A'))
hmap_senescence <- Heatmap(mat_senescence, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                       row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_senescence], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                       cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                       row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                       top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_senescence


qk_assemble(rownames(mat_stemness), pamClusters$clustering[index_stemness])
right_anno <- qk_anno(mat_stemness, c('ODF1', 'PROM1', 'CD40', 'TWIST1', 'KLF6', 'SNAP25', 'CASP4', 'S1PR2', 'BRD3', 'LRRN3', 'TTC22', 'MARCKS', 'MMP10', 'ADO', 'CHST2', 'ZNF202'))
hmap_stemness <- Heatmap(mat_stemness, name = 'Zscore', column_split = factor(split, levels=c("Control", "CD19", "HA")), column_order = col_order,
                           row_gap = unit(1.5, "mm"), row_split = pamClusters$clustering[index_stemness], cluster_row_slices = FALSE, row_dend_reorder = TRUE,
                           cluster_rows = TRUE, show_row_dend = TRUE,show_row_names = FALSE, row_names_gp = gpar(fontsize = 10, fontface = 'bold'), row_names_side = 'right',
                           row_dend_width = unit(5,'mm'), column_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 8, fontface = 'bold'), show_column_names = FALSE,
                           top_annotation = colAnn_noAnno, right_annotation = right_anno, use_raster=FALSE)
hmap_stemness





#### count up cluster representation

library(tidyverse)
library(forcats)

cluster_representation <- list(
  C5 = pamClusters$clustering[index_CARTEx_630],
  CARTEx = pamClusters$clustering[index_CARTEx_200],
  LCMV = pamClusters$clustering[index_Wherry],
  NK_like = pamClusters$clustering[index_NKlike],
  BBD = pamClusters$clustering[index_BBD],
  PD1 = pamClusters$clustering[index_PD1],
  Activation = pamClusters$clustering[index_activation],
  Anergy = pamClusters$clustering[index_anergy],
  Senescence = pamClusters$clustering[index_senescence],
  Stemness = pamClusters$clustering[index_stemness]
)


# Convert the list into a dataframe
df_cluster_representation <- cluster_representation %>% enframe(name = "Signature", value = "Cluster") %>% unnest(Cluster)

# Reverse the order of clusters 
df_cluster_representation$Cluster <- fct_rev(factor(df_cluster_representation$Cluster))

# Set the custom order for the vectors 
df_cluster_representation$Signature <- factor(df_cluster_representation$Signature, levels = rev(c('C5', 'CARTEx', 'LCMV', 'NK_like', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness')))

# Plot the relative distribution of the elements across clusters for each vector 
plot_signature_cluster_representation <- ggplot(df_cluster_representation, aes(x = factor(Signature), fill = factor(Cluster))) +
  geom_bar(position = "fill", colour="black") + labs(x = NULL, y = "Proportion", fill = "Cluster") +
  scale_y_continuous(labels = scales::percent_format()) + scale_fill_manual(values = rev(c("lightsalmon", "lightpink2", "lightpink3", "plum3", "mediumorchid1"))) +
  theme_classic() + coord_flip()

generate_figs(plot_signature_cluster_representation, './plots/plot_signature_cluster_representation', c(2.75,3))












