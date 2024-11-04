

####################################################################################################
####################################### Initialize environment #####################################
####################################################################################################

set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")
setwd(PATH_EXPERIMENTS)


Identifier <- c("GSE125881", "GSE136874", "GSE136184", "GSE126030", "GSE120575", "GSE146264", "GSE151511", "GSE153931", "GSE160160", "GSE196606", "GSE207935", "Zenodo3993994", "mdandersonTCM")

Description <- c("CD19 CAR T cells from cancer patients across 16 weeks", "Functional CD19 CAR T cells vs exhaustion-prone GD2 CAR T cells", "Aging CD8+ T cells from healthy donors across 9 decades", "T cell activation", "TILs from cancer patients after receiving immune checkpoint blockade", "CD8+ T cell from autoimmune patients with psoriasis", "CD19 CAR T cell infusion product from cancer patients", "SARS-CoV-2 reactive T cells from COVID-19 patients", "Continuous antigen exposure in CAR T cells", "T cells from NBEAL2-deficient patients with thrombocytopenia", "T cells from patients with STAT3 gain-of-function syndrome", "T cells from patients with Parkinson's disease", "Integrated atlas of T cells from numerous cancer types")

Study <- c("Sheih et al. 2020 Nature Communications", "Lynn et al. 2019 Nature", "Lu et al. 2022 Nature Communications", "Szabo et al. 2019 Nature Communications", "Sade-Feldman et al. 2018 Cell", "Liu et al. 2021 Journal of Allergy and Clinical Immunology", "Deng et al. 2020 Nature Medicine", "Kusnadi et al. 2021 Science Immunology", "Good et al. 2021 Cell", "Delage et al. 2023 Nature Communications", "Schmitt et al. 2022 JCI Insight", "Wang et al. 2021 Cell Discovery", "Chu et al. 2023 Nature Medicine")

URL <- c("https://pubmed.ncbi.nlm.nih.gov/31924795/", "https://pubmed.ncbi.nlm.nih.gov/31802004/", "https://pubmed.ncbi.nlm.nih.gov/36050300/", "https://pubmed.ncbi.nlm.nih.gov/31624246/", "https://pubmed.ncbi.nlm.nih.gov/30388456/", "https://pubmed.ncbi.nlm.nih.gov/33309739/", "https://pubmed.ncbi.nlm.nih.gov/33020644/", "https://pubmed.ncbi.nlm.nih.gov/33478949/", "https://pubmed.ncbi.nlm.nih.gov/34861191/", "https://pubmed.ncbi.nlm.nih.gov/37349339/", "https://pubmed.ncbi.nlm.nih.gov/36136607/", "https://pubmed.ncbi.nlm.nih.gov/34282123/", "https://pubmed.ncbi.nlm.nih.gov/37248301/")

df_experiments <- data.frame(Identifier, Description, Study, URL)


fetch_qc_data <- function(expt_name){
  return(read.csv(paste('./' ,expt_name, '/data/', expt_name, '_prepare_qc_review.csv', sep = '')))
}

# "","genes","cells"
# "All",5000,16291
# "preQC",5000,4619
# "preQC",5000,4619


# qc_review


reformat <- data.frame()

for (ID in Identifier){
  test <- fetch_qc_data(ID)
  test <- cbind(ID, unname(test[1,2:3]), test[2,3], test[3,3])
  reformat <- rbind(reformat, test)
}

colnames(reformat) <- c('Identifier', 'Genes', 'Total cells', 'Pre-QC CD8+ T cells', 'Post-QC CD8+ T cells')


summarize_expts <- merge(df_experiments, reformat, by = 'Identifier')

write.csv(summarize_expts, "./summarize_expts.csv", row.names=FALSE)










