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

TSG <- toupper(rownames(read.csv(paste0(PATH_SIGNATURES, "tumor-suppressor-genes.csv"), row.names = 1, header = FALSE)))
activation.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "panther-activation.csv", sep = ''), header = FALSE, row.names = 1)))
anergy.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "SAFFORD_T_LYMPHOCYTE_ANERGY.csv", sep = ''), header = FALSE, row.names = 1)))
stemness.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", sep = ''), header = FALSE, row.names = 1)))
senescence.sig <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "M9143_FRIDMAN_SENESCENCE_UP.csv", sep = ''), header = FALSE, row.names = 1)))
NK_like <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "NK-like-dysfunction.csv", sep = ''), header = FALSE, row.names = 1)))
Wherry_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", sep = ''), header = FALSE, row.names = 1)))
BBD_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Selli_2023_Blood_TBBDex.csv", sep = ''), header = FALSE, row.names = 1)))
PD1_Tex <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Cai_2020_Pathology_PD1_Tex.csv", sep = ''), header = FALSE, row.names = 1)))
TSR <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Chu_2023_Nat_Med_T_stress_response.csv", sep = ''), header = FALSE, row.names = 1)))

Daniel_Tex_Term <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Daniel_2022_Nat_Imm_TexTerm.csv", sep = ''), header = FALSE, row.names = 1)))
Daniel_Tex_KLR <- toupper(rownames(read.csv(paste(PATH_SIGNATURES, "Daniel_2022_Nat_Imm_TexKLR.csv", sep = ''), header = FALSE, row.names = 1)))


# check overlap with Tex subsets from Daniel et al. 2022
intersect(Daniel_Tex_Term, Daniel_Tex_KLR)
intersect(toupper(rownames(cartex_200_weights)),Daniel_Tex_Term)
intersect(toupper(rownames(cartex_200_weights)), Daniel_Tex_KLR)





expt.list <- c("GSE125881", "GSE136874", "GSE136184", "GSE126030", "GSE120575", "GSE146264", "GSE151511", "GSE153931", "GSE160160", "GSE196606", "GSE207935", "Zenodo3993994", "mdandersonTCM")

# , "GSE235676"

# create dictionary hash containing all genes for each experiment in expt.list
expt.allgenes <- c()
expt.CARTEx630 <- c()
expt.CARTEx200 <- c()
expt.CARTEx84 <- c()
expt.TSG <- c()
expt.activation <- c()
expt.anergy <- c()
expt.stemness <- c()
expt.senescence <- c()
expt.NK_like <- c()
expt.Wherry_Tex <- c()
expt.BBD_Tex <- c()
expt.PD1_Tex <- c()
expt.TSR <- c()
expt.Tex_Term <- c()
expt.Tex_KLR <- c()

expt.allgenes.counts <- c()
expt.CARTEx630.counts <- c()
expt.CARTEx200.counts <- c()
expt.CARTEx84.counts <- c()
expt.TSG.counts <- c()
expt.activation.counts <- c()
expt.anergy.counts <- c()
expt.stemness.counts <- c()
expt.senescence.counts <- c()
expt.NK_like.counts <- c()
expt.Wherry_Tex.counts <- c()
expt.BBD_Tex.counts <- c()
expt.PD1_Tex.counts <- c()
expt.TSR.counts <- c()
expt.Tex_Term.counts <- c()
expt.Tex_KLR.counts <- c()


for (exptID in expt.list){
  print(exptID)
  expt.allgenes[[exptID]] <- toupper(read.csv(paste(PATH_EXPERIMENTS, exptID, "/data/", exptID, "_allgenes.csv", sep = ''), header = TRUE, row.names = 1)$x)
  expt.allgenes.counts[[exptID]] <- length(expt.allgenes[[exptID]])
  print(paste0("all genes: ", expt.allgenes.counts[[exptID]]))
  
  expt.CARTEx630[[exptID]] <- intersect(expt.allgenes[[exptID]], toupper(rownames(cartex_630_weights)))
  expt.CARTEx200[[exptID]] <- intersect(expt.allgenes[[exptID]], toupper(rownames(cartex_200_weights)))
  expt.CARTEx84[[exptID]] <- intersect(expt.allgenes[[exptID]], toupper(rownames(cartex_84_weights)))
  
  expt.CARTEx630.counts[[exptID]] <- length(expt.CARTEx630[[exptID]])
  expt.CARTEx200.counts[[exptID]] <- length(expt.CARTEx200[[exptID]])
  expt.CARTEx84.counts[[exptID]] <- length(expt.CARTEx84[[exptID]])
  
  print(paste0("CARTEx630 genes: ", expt.CARTEx630.counts[[exptID]]))
  print(paste0("CARTEx200 genes: ", expt.CARTEx200.counts[[exptID]]))
  print(paste0("CARTEx84 genes: ", expt.CARTEx84.counts[[exptID]]))
  
  expt.TSG[[exptID]] <- intersect(expt.allgenes[[exptID]], TSG)
  expt.TSG.counts[[exptID]] <- length(expt.TSG[[exptID]])
  print(paste0("TSG genes: ", expt.TSG.counts[[exptID]]))
  
  expt.activation[[exptID]] <- intersect(expt.allgenes[[exptID]], activation.sig)
  expt.activation.counts[[exptID]] <- length(expt.activation[[exptID]])
  print(paste0("Activation genes: ", expt.activation.counts[[exptID]]))
  
  expt.anergy[[exptID]] <- intersect(expt.allgenes[[exptID]], anergy.sig)
  expt.anergy.counts[[exptID]] <- length(expt.anergy[[exptID]])
  print(paste0("Anergy genes: ", expt.anergy.counts[[exptID]]))
  
  expt.stemness[[exptID]] <- intersect(expt.allgenes[[exptID]], stemness.sig)
  expt.stemness.counts[[exptID]] <- length(expt.stemness[[exptID]])
  print(paste0("Stemness genes: ", expt.stemness.counts[[exptID]]))
  
  expt.senescence[[exptID]] <- intersect(expt.allgenes[[exptID]], senescence.sig)
  expt.senescence.counts[[exptID]] <- length(expt.senescence[[exptID]])
  print(paste0("Senescence genes: ", expt.senescence.counts[[exptID]]))
  
  expt.NK_like[[exptID]] <- intersect(expt.allgenes[[exptID]], NK_like)
  expt.NK_like.counts[[exptID]] <- length(expt.NK_like[[exptID]])
  print(paste0("NK_like genes: ", expt.NK_like.counts[[exptID]]))
  
  expt.Wherry_Tex[[exptID]] <- intersect(expt.allgenes[[exptID]], Wherry_Tex)
  expt.Wherry_Tex.counts[[exptID]] <- length(expt.Wherry_Tex[[exptID]])
  print(paste0("Wherry_Tex genes: ", expt.Wherry_Tex.counts[[exptID]]))
  
  expt.BBD_Tex[[exptID]] <- intersect(expt.allgenes[[exptID]], BBD_Tex)
  expt.BBD_Tex.counts[[exptID]] <- length(expt.BBD_Tex[[exptID]])
  print(paste0("BBD_Tex genes: ", expt.BBD_Tex.counts[[exptID]]))
  
  expt.PD1_Tex[[exptID]] <- intersect(expt.allgenes[[exptID]], PD1_Tex)
  expt.PD1_Tex.counts[[exptID]] <- length(expt.PD1_Tex[[exptID]])
  print(paste0("PD1_Tex genes: ", expt.PD1_Tex.counts[[exptID]]))
  
  expt.TSR[[exptID]] <- intersect(expt.allgenes[[exptID]], TSR)
  expt.TSR.counts[[exptID]] <- length(expt.TSR[[exptID]])
  print(paste0("TSR genes: ", expt.TSR.counts[[exptID]]))
  
  expt.Tex_Term[[exptID]] <- intersect(expt.allgenes[[exptID]], Daniel_Tex_Term)
  expt.Tex_Term.counts[[exptID]] <- length(expt.Tex_Term[[exptID]])
  print(paste0("Tex_Term genes: ", expt.Tex_Term.counts[[exptID]]))
  
  expt.Tex_KLR[[exptID]] <- intersect(expt.allgenes[[exptID]], Daniel_Tex_KLR)
  expt.Tex_KLR.counts[[exptID]] <- length(expt.Tex_KLR[[exptID]])
  print(paste0("Tex_KLR genes: ", expt.Tex_KLR.counts[[exptID]]))
  
  print("")
}

expt.matrix.counts <- do.call("cbind", list(expt.allgenes.counts, expt.CARTEx630.counts, expt.CARTEx200.counts, expt.CARTEx84.counts, expt.TSG.counts, 
                                            expt.activation.counts, expt.anergy.counts, expt.stemness.counts, expt.senescence.counts, 
                                            expt.NK_like.counts, expt.Wherry_Tex.counts, expt.BBD_Tex.counts, expt.PD1_Tex.counts, 
                                            expt.TSR.counts, expt.Tex_Term.counts, expt.Tex_KLR.counts))
expt.matrix.counts

write.csv(expt.matrix.counts, paste0(PATH_CONSTRUCTION, "data/expt_matrix_counts.csv"), row.names=TRUE)




####################################################################################################
################################## Generate gene representation maps ###############################
####################################################################################################

library(tidyverse)
# https://stackoverflow.com/questions/47733031/how-to-plot-dataframe-in-r-as-a-heatmap-grid

expt.df.counts <- as.data.frame(expt.matrix.counts)
colnames(expt.df.counts) <- c('All', 'C5', 'CARTEx', 'CARTEx_84', 'TSG', 'Activation', 'Anergy', 'Stemness', 'Senescence', 'NK-like', 'LCMV', 'BBD', 'PD1', 'Stress Response', 'Tex-Term', 'Tex-KLR')

# store all counts
expt.df.counts.all <- expt.df.counts['All']

# rearrange and remove CARTEx_84 and TSG
expt.df.counts <- expt.df.counts[c('C5', 'CARTEx', 'NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness', 'Stress Response', 'Tex-Term', 'Tex-KLR')]

# define a vector with the denominator values and divide each column by the corresponding element in the vector
divisors <- c(630, 200, length(NK_like), length(Wherry_Tex), length(BBD_Tex), length(PD1_Tex), length(activation.sig), length(anergy.sig), length(senescence.sig), length(stemness.sig), length(TSR), length(Daniel_Tex_Term), length(Daniel_Tex_KLR))
expt.df.counts[] <- lapply(expt.df.counts, as.numeric)
expt.df.scaled <- sweep(expt.df.counts, 2, divisors, "/")


expt.heatmap <- expt.df.scaled %>% rownames_to_column() %>% gather(colname, value, -rowname)
expt.heatmap$colname <- factor(expt.heatmap$colname, ordered=TRUE, levels = c('C5', 'CARTEx', 'NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness', 'Stress Response', 'Tex-Term', 'Tex-KLR'))
expt.heatmap$rowname <- factor(expt.heatmap$rowname, ordered=TRUE, levels = c('GSE120575', 'GSE125881', 'GSE126030', 'GSE136184', 'GSE136874', 'GSE146264', 'GSE151511', 'GSE153931', 'GSE160160', 'GSE196606', 'GSE207935', 'Zenodo3993994', 'mdandersonTCM'))

x_labels <- paste0(levels(expt.heatmap$rowname), " (", as.character(unlist(expt.df.counts.all, use.names = FALSE)), ")")
y_labels <- paste0(levels(expt.heatmap$colname), " (", as.character(divisors), ")")



all_genes_repmap <- ggplot(expt.heatmap, aes(x = rowname, y = colname, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colours = c("white", "lightblue"), values = c(0,1)) +
  geom_text(aes(fill = expt.heatmap$value, label = round(expt.heatmap$value,2))) +
  theme_classic() + theme(axis.text.x = element_text(angle=25, vjust=1, hjust=1, size = 10, color = 'black'), axis.text.y = element_text(size = 10, color = 'black'), legend.position="none") + 
  xlab(NULL) + ylab(NULL) + scale_x_discrete(labels = x_labels) + scale_y_discrete(labels = y_labels)

generate_figs(all_genes_repmap, "./plots/all_genes_repmap", c(7, 3.5))




###
###








CARTEx_versions <- list(
  CARTEx_630 = rownames(cartex_630_weights),
  CARTEx_200 = rownames(cartex_200_weights),
  CARTEx_84 = rownames(cartex_84_weights)
)

# upset_CARTEx_versions <- as.grob(upset(fromList(CARTEx_versions), order.by = "freq"))


# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# https://github.com/krassowski/complex-upset



upset_expt_allgenes <- ComplexUpset::upset(UpSetR::fromList(expt.allgenes), intersect = names(expt.allgenes),
                                           base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_expt_allgenes, "./plots/upset_expt_allgenes", c(16, 6))

upset_expt_CARTEx630 <- ComplexUpset::upset(UpSetR::fromList(expt.CARTEx630), intersect = names(expt.CARTEx630),
                                            base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_expt_CARTEx630, "./plots/upset_expt_CARTEx630", c(6, 6))

upset_expt_CARTEx200 <- ComplexUpset::upset(UpSetR::fromList(expt.CARTEx200), intersect = names(expt.CARTEx200),
                                            base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_expt_CARTEx200, "./plots/upset_expt_CARTEx200", c(6, 6))

upset_expt_CARTEx84 <- ComplexUpset::upset(UpSetR::fromList(expt.CARTEx84), intersect = names(expt.CARTEx84),
                                           base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_expt_CARTEx84, "./plots/upset_expt_CARTEx84", c(6, 6))

upset_expt_TSG <- ComplexUpset::upset(UpSetR::fromList(expt.TSG), intersect = names(expt.TSG),
                                      base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_expt_TSG, "./plots/upset_expt_TSG", c(6, 6))


# check HNF1 detection
for (exptID in expt.list){
  print(exptID)
  print(intersect(expt.allgenes[[exptID]], "HNF1A"))
  print(intersect(expt.allgenes[[exptID]], "HNF1B"))
}


# detection_CARTEx_genes

expt.matrix.counts
plot(expt.matrix.counts[,1], expt.matrix.counts[,2], col = 'red', pch = 19, ylim = c(0, 650), xlab = 'All genes', ylab = 'CARTEx genes') +
  points(expt.matrix.counts[,1], expt.matrix.counts[,3], col = 'purple', pch = 19) +
  points(expt.matrix.counts[,1], expt.matrix.counts[,4], col = 'blue', pch = 19) +
  abline(h = 630, col = "black", lty = 2) + abline(h = 200, col = "black", lty = 2) + abline(h = 84, col = "black", lty = 2)

# detection_state_genes

plot(expt.matrix.counts[,1], expt.matrix.counts[,6], col = 'red', pch = 19, ylim = c(0, 200), xlab = 'All genes', ylab = 'State genes') +
  points(expt.matrix.counts[,1], expt.matrix.counts[,7], col = 'purple', pch = 19) +
  points(expt.matrix.counts[,1], expt.matrix.counts[,8], col = 'blue', pch = 19) +
  points(expt.matrix.counts[,1], expt.matrix.counts[,9], col = 'orange', pch = 19) +
  abline(h = 61, col = "black", lty = 2) + abline(h = 90, col = "black", lty = 2) + abline(h = 199, col = "black", lty = 2) + abline(h = 77, col = "black", lty = 2)

length(activation.sig)
length(anergy.sig)
length(stemness.sig)
length(senescence.sig)

# detection_exhaustion_genes

plot(expt.matrix.counts[,1], expt.matrix.counts[,10], col = 'red', pch = 19, ylim = c(0, 150), xlab = 'All genes', ylab = 'Exhaustion genes') +
  points(expt.matrix.counts[,1], expt.matrix.counts[,11], col = 'purple', pch = 19) +
  points(expt.matrix.counts[,1], expt.matrix.counts[,12], col = 'blue', pch = 19) +
  points(expt.matrix.counts[,1], expt.matrix.counts[,13], col = 'orange', pch = 19) +
  abline(h = 29, col = "black", lty = 2) + abline(h = 107, col = "black", lty = 2) + abline(h = 145, col = "black", lty = 2) + abline(h = 78, col = "black", lty = 2)

length(NK_like)
length(Wherry_Tex)
length(BBD_Tex)
length(PD1_Tex)


exhaustion_sigs <- list(
  C5 = toupper(rownames(cartex_630_weights)),
  CARTEx = toupper(rownames(cartex_200_weights)),
  NK_like = NK_like,
  LCMV = Wherry_Tex,
  BBD = BBD_Tex,
  PD1 = PD1_Tex
)

upset_exhaustion_sigs <- ComplexUpset::upset(UpSetR::fromList(exhaustion_sigs), intersect = names(exhaustion_sigs), width_ratio=0.15,
                                             base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_exhaustion_sigs, "./plots/upset_exhaustion_sigs", c(6, 3.5))


state_sigs <- list(
  C5 = toupper(rownames(cartex_630_weights)),
  CARTEx = toupper(rownames(cartex_200_weights)),
  Activation = activation.sig,
  Anergy = anergy.sig,
  Senescence = senescence.sig,
  Stemness = stemness.sig
)

upset_state_sigs <- ComplexUpset::upset(UpSetR::fromList(state_sigs), intersect = names(state_sigs), width_ratio=0.25,
                                        base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_state_sigs, "./plots/upset_state_sigs", c(4.5, 3.5))



exhaustion_state_sigs <- list(
  C5 = toupper(rownames(cartex_630_weights)),
  CARTEx = toupper(rownames(cartex_200_weights)),
  NK_like = NK_like,
  LCMV = Wherry_Tex,
  BBD = BBD_Tex,
  PD1 = PD1_Tex,
  Activation = activation.sig,
  Anergy = anergy.sig,
  Senescence = senescence.sig,
  Stemness = stemness.sig
)

upset_exhaustion_state_sigs <- ComplexUpset::upset(UpSetR::fromList(exhaustion_state_sigs), intersect = names(exhaustion_state_sigs), width_ratio=0.08, height_ratio = 1.25,
                                        base_annotations=list('Intersection size'=(intersection_size(counts=FALSE))))
generate_figs(upset_exhaustion_state_sigs, "./plots/upset_exhaustion_state_sigs", c(11, 3.5))




####################################################################################################
#################################### Generate representation charts ################################
####################################################################################################

# progress bar charts?

exhaustion_state_sigs <- list(
  CARTEx = rownames(cartex_200_weights),
  NKlike_Tex = NK_like,
  LCMV_Tex = Wherry_Tex,
  BBD_Tex = BBD_Tex,
  PD1_Tex = PD1_Tex,
  Activation = activation.sig,
  Anergy = anergy.sig,
  Senescence = senescence.sig,
  Stemness = stemness.sig
)

gene_counts <- sapply(exhaustion_state_sigs, length)

count_intersections <- function(vector, given_vector) {
  length(intersect(vector, given_vector))
}

reference_genes <- rownames(cartex_630_weights)

intersection_counts <- sapply(exhaustion_state_sigs, count_intersections, reference_genes)


# fraction of signature genes represented in C5
fraction_signature <- intersection_counts/gene_counts
fraction_signature_df <- as.data.frame(fraction_signature[-1])
colnames(fraction_signature_df) <- 'value'

fraction_sig_repchart <- ggplot(fraction_signature_df, aes(x=factor(rownames(fraction_signature_df), levels = rownames(fraction_signature_df)), y=value)) + 
  geom_bar(stat = "identity") + ylim(c(0,1)) + theme_classic() + ylab("Fraction") + xlab(NULL) + 
  scale_x_discrete(labels = c('NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness')) +
  geom_text(aes(label=round(value,2)), vjust=-1)
generate_figs(fraction_sig_repchart, "./plots/fraction_sig_repchart", c(6,2))


# fraction of signature genes represented in CARTEx

fraction_signature <- sapply(exhaustion_state_sigs, count_intersections, rownames(cartex_200_weights))/gene_counts
fraction_signature_df <- as.data.frame(fraction_signature[-1])
colnames(fraction_signature_df) <- 'value'

fraction_sig_repchart_CARTEx <- ggplot(fraction_signature_df, aes(x=factor(rownames(fraction_signature_df), levels = rownames(fraction_signature_df)), y=value)) + 
  geom_bar(stat = "identity") + ylim(c(0,1)) + theme_classic() + ylab("Fraction") + xlab(NULL) + 
  scale_x_discrete(labels = c('NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness')) +
  geom_text(aes(label=round(value,2)), vjust=-1)
generate_figs(fraction_sig_repchart_CARTEx, "./plots/fraction_sig_repchart_CARTEx", c(6,2))


# fraction of C5 genes represented by signature

fraction_CARTEx_630 <- intersection_counts/length(reference_genes)
fraction_CARTEx_630_df <- as.data.frame(fraction_CARTEx_630)
colnames(fraction_CARTEx_630_df) <- 'value'

fraction_C5_repchart <- ggplot(fraction_CARTEx_630_df, aes(x=factor(rownames(fraction_CARTEx_630_df), levels = rownames(fraction_CARTEx_630_df)), y=value)) + 
  geom_bar(stat = "identity") + ylim(c(0,1)) + theme_classic() + ylab("Fraction") + xlab(NULL) + 
  scale_x_discrete(labels = c('CARTEx', 'NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness')) +
  geom_text(aes(label=round(value,2)), vjust=-1)
generate_figs(fraction_C5_repchart, "./plots/fraction_C5_repchart", c(6,2))


# fraction of CARTEx genes represented in signature

fraction_CARTEx_200 <- intersection_counts/length(rownames(cartex_200_weights))
fraction_CARTEx_200_df <- as.data.frame(fraction_CARTEx_200[-1])
colnames(fraction_CARTEx_200_df) <- 'value'

fraction_CARTEx_repchart <- ggplot(fraction_CARTEx_200_df, aes(x=factor(rownames(fraction_CARTEx_200_df), levels = rownames(fraction_CARTEx_200_df)), y=value)) + 
  geom_bar(stat = "identity") + ylim(c(0,1)) + theme_classic() + ylab("Fraction") + xlab(NULL) + 
  scale_x_discrete(labels = c('NK-like', 'LCMV', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness')) +
  geom_text(aes(label=round(value,2)), vjust=-1)
generate_figs(fraction_CARTEx_repchart, "./plots/fraction_CARTEx_repchart", c(6,2))


# Increase margin size
# https://r-graph-gallery.com/210-custom-barplot-layout.html#las
# par(mar=c(7,5,1,1))






