
setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")

# load signatures

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



# Function to compute Jaccard Index
jaccard_index <- function(set1, set2) {
  intersection_length <- length(intersect(set1, set2))
  union_length <- length(union(set1, set2))
  if (union_length == 0) return(0)  # Avoid division by zero
  return(intersection_length / union_length)
}

# Permutation test function for gene sets
perm_test_gene_sets <- function(set1, set2, num_permutations = 10000) {
  # Observed Jaccard Index
  observed_jaccard <- jaccard_index(set1, set2)
  
  # Combine all genes for shuffling
  all_genes <- union(set1, set2)
  
  # Perform permutations
  permutations <- numeric(num_permutations)
  for (i in 1:num_permutations) {
    # Shuffle one set while maintaining its size
    permuted_set2 <- sample(all_genes, size = length(set2), replace = FALSE)
    permutations[i] <- jaccard_index(set1, permuted_set2)
  }
  
  # Compute p-value (one-tailed test for simplicity; adjust if two-tailed is needed)
  p_value <- sum(permutations >= observed_jaccard) / num_permutations
  
  # Return results
  return(list(
    observed_jaccard = observed_jaccard,
    p_value = p_value,
    permutations = permutations
  ))
}


# test scenario
test_result <- perm_test_gene_sets(rownames(cartex_200_weights), NK_like)
test_result <- perm_test_gene_sets(rownames(cartex_630_weights), rownames(cartex_200_weights))


print(paste("Observed Jaccard Index:", test_result$observed_jaccard))
print(paste("P-value:", test_result$p_value))




sig_list <- data.frame(
  Signature = c('C5', 'CARTEx', 'LCMV', 'NK_like', 'BBD', 'PD1', 'Activation', 'Anergy', 'Senescence', 'Stemness', 'Stress_Response', 'Tex_Term', 'Tex_KLR'),
  Genes = I(list(rownames(cartex_630_weights), rownames(cartex_200_weights), Wherry_Tex, NK_like, BBD_Tex, PD1_Tex, activation.sig, anergy.sig, senescence.sig, stemness.sig, TSR, Daniel_Tex_Term, Daniel_Tex_KLR))
)


# Function to compute Jaccard and p-value for all pairs
compute_metrics <- function(df) {
  n <- nrow(df)
  jaccard_matrix <- matrix(0, nrow = n, ncol = n)
  p_value_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {  # Don't compare a set to itself
        result <- perm_test_gene_sets(df$Genes[[i]], df$Genes[[j]])
        jaccard_matrix[i, j] <- result$observed_jaccard
        p_value_matrix[i, j] <- result$p_value
      } else {
        jaccard_matrix[i, j] <- 1  # Jaccard of a set with itself is 1
        p_value_matrix[i, j] <- 1  # P-value for self-comparison can be set to 1 or NA
      }
    }
  }
  
  # Convert matrices to data frames for joining with original df
  jaccard_df <- as.data.frame(jaccard_matrix) %>% 
    setNames(paste0("Jaccard_", df$Signature))
  p_value_df <- as.data.frame(p_value_matrix) %>% 
    setNames(paste0("P_value_", df$Signature))
  
  # Combine with original data
  df <- cbind(df, jaccard_df, p_value_df)
  
  return(df)
}


sig_list <- compute_metrics(sig_list)

print(sig_list)

library(pheatmap)
library(dplyr)

prep_for_heatmaps <- function(df) {
  # Assuming df has columns for Jaccard Index and p-values named like "Jaccard_SignatureX", "P_value_SignatureX"
  
  # Extract Jaccard Index data
  jaccard_cols <- grep("^Jaccard_", names(df), value = TRUE)
  jaccard_data <- df[, jaccard_cols]
  
  # Extract p-value data
  p_value_cols <- grep("^P_value_", names(df), value = TRUE)
  p_value_data <- df[, p_value_cols]
  
  # Convert to matrices for heatmap functions
  jaccard_matrix <- as.matrix(jaccard_data)
  p_value_matrix <- as.matrix(p_value_data)
  
  # Naming rows and columns
  rownames(jaccard_matrix) <- df$Signature
  colnames(jaccard_matrix) <- gsub("Jaccard_", "", colnames(jaccard_matrix))
  rownames(p_value_matrix) <- df$Signature
  colnames(p_value_matrix) <- gsub("P_value_", "", colnames(p_value_matrix))

  return(list(jaccard_matrix = jaccard_matrix, p_value_matrix = p_value_matrix))
}

result <- prep_for_heatmaps(sig_list)



# The Jaccard Index measures the similarity between sets based on their overlapping genes, while the p-value is derived from a permutation test to assess whether the observed overlap is statistically significant. 




# Generate Jaccard Index heatmap
JIheatmap <- pheatmap(result$jaccard_matrix, 
         # main = "Jaccard Index Heatmap",
         color = colorRampPalette(c("white", "blue"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 10, legend = FALSE)

generate_figs(JIheatmap, "./plots/JIheatmap", c(8,5))

###
###


# Generate p-value heatmap
pheatmap(result$p_value_matrix, 
         main = "P-value Heatmap",
         color = colorRampPalette(c("yellow", "red"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.3f",
         fontsize_number = 6)






