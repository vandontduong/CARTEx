
setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")


cartex_200_weights <- read.csv(paste(PATH_WEIGHTS, "cartex-200-weights.csv", sep = ''), header = TRUE, row.names = 1)

library(readxl)
library(stringr)

filename <- "GOrilla_function_CARTEx_200_results.xlsx"
GOrilla_function_results <- read_excel(paste0('./data/', filename), sheet = 'results', col_names = TRUE)

# inspect modules
GOrilla_function_results
GOrilla_function_results$`GO term`
GOrilla_function_results[c('GO term', 'Description')]




# check all genes

list_of_module_genes_all <- list()
list_of_module_genesymbols_all <- list()

for (sheet in c(GOrilla_function_results$`GO term`)){
  sheet_name <- gsub(":", "_", sheet)
  list_of_module_genes_all[[sheet]] <- c(read_excel(paste0('./data/', filename), sheet = sheet_name, col_names = FALSE))
  list_of_module_genesymbols_all[[sheet]] <- c(sapply(list_of_module_genes_all[[sheet]], function(x) str_extract(x, "^[^ -]+")))
}

# there are 185 genes in the whole map
all_function_module_genes <- unique(unlist(list_of_module_genesymbols_all))

missing_function_module_genes <- setdiff(rownames(cartex_200_weights), all_function_module_genes)



# select modules (TERMINAL)
select_function_modules <- c("GO:0004666", "GO:0023024", "GO:0008009", "GO:0051019", "GO:0004879", "GO:0001228")
# GO:0048020 is a subset of GO:0008009
# GO:0017017 and GO:0008330 are a subset of GO:0051019
# GO:0046965 and GO:0035259 a subset of GO:0004879
# note, there is some overlap between GO:0004879 and GO:0001228
select_function_modules_mod <- gsub(":", "_", select_function_modules)

# reinspect select modules
GOrilla_function_results[GOrilla_function_results$`GO term` %in% select_function_modules, ]

list_of_module_genes <- list()
list_of_module_genesymbols <- list()

for (sheet in select_function_modules){
  sheet_name <- gsub(":", "_", sheet)
  list_of_module_genes[[sheet]] <- c(read_excel(paste0('./data/', filename), sheet = sheet_name, col_names = FALSE))
  list_of_module_genesymbols[[sheet]] <- c(sapply(list_of_module_genes[[sheet]], function(x) str_extract(x, "^[^ -]+")))
}

# 32 unique genes from the select modules
select_function_module_genes <- unique(unlist(list_of_module_genesymbols))




# examine what is missing

# 153 genes are not in the terminal modules
missing_function_module_genes <- setdiff(all_function_module_genes, select_function_module_genes)

select_toplvl_function_modules <- c('GO:0005488', 'GO:0098772', 'GO:0005102', 'GO:0048018')

genes_binding <- unname(unlist(c(list_of_module_genesymbols_all['GO:0005488']))) # 183
genes_molecular_function_regulator <- unname(unlist(c(list_of_module_genesymbols_all['GO:0098772']))) # 40
genes_signaling_receptor_binding <- unname(unlist(c(list_of_module_genesymbols_all['GO:0005102']))) # 39
genes_receptor_ligand_activity <- unname(unlist(c(list_of_module_genesymbols_all['GO:0048018']))) # 39


intersect(missing_function_module_genes, genes_binding) # 151
intersect(missing_function_module_genes, genes_molecular_function_regulator) # 33

intersect(genes_binding, genes_molecular_function_regulator) # 38
intersect(genes_binding, genes_signaling_receptor_binding) # 39
intersect(genes_molecular_function_regulator, genes_signaling_receptor_binding) # 21

# check overlap
test <- list(
  binding = genes_binding,
  MFR = genes_molecular_function_regulator,
  SRB = genes_signaling_receptor_binding, 
  RLA = genes_receptor_ligand_activity
)

library(ggvenn)
ggvenn(test)

# create gene lists from venn diagram based on visual inspection and grouping

# RLA has 19 genes overlapping with SRB
intersect(genes_receptor_ligand_activity, genes_signaling_receptor_binding)

# SRB has 20 genes that are not part of RLA
setdiff(genes_signaling_receptor_binding, genes_receptor_ligand_activity)


# group MFR and SRB as cell signal regulators - 58 total genes
genes_cell_signal_regulators <- union(genes_molecular_function_regulator, genes_signaling_receptor_binding)

# genes not in terminal modules
genes_other_CSR <- intersect(missing_function_module_genes, genes_cell_signal_regulators) # 44
genes_other_binding <- intersect(missing_function_module_genes, setdiff(genes_binding, genes_cell_signal_regulators)) # 109

# check that there are 185 genes accounted for from the GOrilla function map
length(select_function_module_genes) + length(genes_other_CSR) + length(genes_other_binding)


# generate representation map

GOrilla_function_results[GOrilla_function_results$`GO term` %in% select_function_modules, ]


CARTEx_200_composite_genes <- list(
  Chemokines = unname(unlist(list_of_module_genesymbols['GO:0008009'])),
  MHC_class_I = unname(unlist(list_of_module_genesymbols['GO:0023024'])),
  Transcription = union(unname(unlist(list_of_module_genesymbols['GO:0001228'])), unname(unlist(list_of_module_genesymbols['GO:0004879']))),
  Oxidoreductase = unname(unlist(list_of_module_genesymbols['GO:0004666'])),
  MAPK = unname(unlist(list_of_module_genesymbols['GO:0051019'])),
  other_CSR = genes_other_CSR,
  other_binding = genes_other_binding
)

# check again for 185 genes from GOrilla map
length(unique(unlist(CARTEx_200_composite_genes)))

gene_counts <- sapply(CARTEx_200_composite_genes, length)


# fraction of signature genes represented in C5
fraction_GOrilla_function <- gene_counts / 185
fraction_GOrilla_function_df <- as.data.frame(fraction_GOrilla_function)
colnames(fraction_GOrilla_function_df) <- 'value'

fraction_GOrilla_function_repchart <- ggplot(fraction_GOrilla_function_df, aes(x=factor(rownames(fraction_GOrilla_function_df), levels = rownames(fraction_GOrilla_function_df)), y=value)) + 
  geom_bar(stat = "identity") + ylim(c(0,1)) + theme_classic() + ylab("Fraction") + xlab(NULL) + 
  scale_x_discrete(labels = c('Chemokine\nactivity', 'Natural killer\nreceptors', 'Transcription\nactivity', 'Oxidoreductase\nactivity', 'MAPK\nactivity', 'Other signal\nregulation', 'Other\nbinding')) +
  geom_text(aes(label=round(value,2)), vjust=-1)
# generate_figs(fraction_sig_repchart, "./plots/fraction_sig_repchart", c(6,2))










