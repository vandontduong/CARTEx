# examine-cartex-biology
# Vandon Duong

setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)

library(ggplotify)
library(eulerr)
library(UpSetR)
library(ggvenn)
library(RColorBrewer)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

euler_fitting = function(nested_list){
  euler_fit <- euler(nested_list)
  euler_plot <- plot(euler_fit, shape = 'ellipse', quantities = list(type = 'counts'), fills = brewer.pal(length(nested_list), 'Set3'))
  euler_error <- error_plot(euler_fit)
  return(list(plot = euler_plot, error = euler_error))
}

####################################################################################################
########################################## Subset overlaps #########################################
####################################################################################################

CARTEx_subsets <- list(
  Ridge = rownames(read.csv("./data/weights_ridge_reg.csv", row.names=1, header = TRUE)),
  Quasi_Ridge = rownames(read.csv("./data/weights_quasi_ridge_reg.csv", row.names=1, header = TRUE)),
  Elastic = rownames(read.csv("./data/weights_elastic_reg.csv", row.names=1, header = TRUE)),
  Lasso = rownames(read.csv("./data/weights_lasso_reg.csv", row.names=1, header = TRUE))
)

euler_CARTEx_subsets <- euler_fitting(CARTEx_subsets)
generate_figs(euler_CARTEx_subsets$plot, "./plots/euler_CARTEx_subsets_plot")
generate_figs(euler_CARTEx_subsets$error, "./plots/euler_CARTEx_subsets_error")

####################################################################################################
########################################## Biology overlap #########################################
####################################################################################################

CARTEx_630 = rownames(read.csv("../weights/cartex-630-weights.csv", row.names=1, header = TRUE))
CARTEx_200 = rownames(read.csv("../weights/cartex-200-weights.csv", row.names=1, header = TRUE))
CARTEx_84 = rownames(read.csv("../weights/cartex-84-weights.csv", row.names=1, header = TRUE))

GO_hypoxia = rownames(read.csv("../signatures/GO_0071456_cellular_response_to_hypoxia.csv", row.names=1, header = TRUE))
GO_apoptosis = rownames(read.csv("../signatures/GO_0042981_regulation_of_apoptotic_process.csv", row.names=1, header = TRUE))
GO_necrosis = rownames(read.csv("../signatures/GO_0001906_cell_killing.csv", row.names=1, header = TRUE))
GO_muscle = rownames(read.csv("../signatures/GO_0060537_muscle_tissue_development.csv", row.names = 1, header = TRUE))
M_matrisome = rownames(read.csv("../signatures/M5885_NABA_MATRISOME_ASSOCIATED.csv", row.names = 1, header = TRUE))
D_isoproterenol = rownames(read.csv("../signatures/D007545_isoproterenol.csv", row.names = 1, header = TRUE))

NK_like = rownames(read.csv("../signatures/NK-like-dysfunction.csv", row.names=1, header = TRUE))
Wherry_Tex = rownames(read.csv("../signatures/Wherry_2007_Immunity_LCMV_Tex_humanized_version.csv", row.names = 1, header = TRUE))
BBD_Tex = rownames(read.csv("../signatures/Selli_2023_Blood_TBBDex.csv", row.names = 1, header = TRUE))
PD1_Tex = rownames(read.csv("../signatures/Cai_2020_Pathology_PD1_Tex.csv", row.names = 1, header = TRUE))

Activation = rownames(read.csv("../signatures/panther-activation.csv", row.names=1, header = TRUE))
Anergy = rownames(read.csv("../signatures/SAFFORD_T_LYMPHOCYTE_ANERGY.csv", row.names = 1, header = TRUE))
Stemness = rownames(read.csv("../signatures/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.csv", row.names = 1, header = TRUE))
Senescence = rownames(read.csv("../signatures/M9143_FRIDMAN_SENESCENCE_UP.csv", row.names = 1, header = TRUE))

Treg_UP = rownames(read.csv("../signatures/M5725_GSE7852_TREG_VS_TCONV_LN_UP.csv", row.names = 1, header = TRUE))
Treg_DN = rownames(read.csv("../signatures/M5674_GSE7460_TCONV_VS_TREG_LN_DN.csv", row.names = 1, header = TRUE))
Treg = intersect(Treg_UP, Treg_DN)

IFNg <- c("FNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "HLA-DRA1")
TNFa <- rownames(read.csv("../signatures/M5890_HALLMARK_TNFA_SIGNALING_VIA_NFKB.csv", row.names = 1, header = TRUE))

Zheng_2017_Cell <- rownames(read.csv("../signatures/Zheng_2017_Cell_Tex_signature.csv", row.names = 1, header = TRUE))
Li_2019_Cell <- rownames(read.csv("../signatures/Li_2019_Cell_Tex_signature.csv", row.names = 1, header = TRUE))
Guo_2018_NatMed <- rownames(read.csv("../signatures/Guo_2018_NatMed_Tex_signature.csv", row.names = 1, header = TRUE))
Zhang_2018_Nature <- rownames(read.csv("../signatures/Zhang_2018_Nature_Tex_signature.csv", row.names = 1, header = TRUE))

TSG = rownames(read.csv("../signatures/tumor-suppressor-genes.csv", row.names = 1, header = TRUE))


CARTEx_composite <- list(
  CARTEx = CARTEx_84,
  Hypoxia = intersect(GO_hypoxia, CARTEx_84),
  Apoptosis = intersect(GO_apoptosis, CARTEx_84),
  Necrosis = intersect(GO_necrosis, CARTEx_84),
  Adrenergic_stress = intersect(D_isoproterenol, CARTEx_84)
)

euler_CARTEx_composite <- euler_fitting(CARTEx_composite)
generate_figs(euler_CARTEx_composite$plot, "./plots/euler_CARTEx_composite_plot")
generate_figs(euler_CARTEx_composite$error, "./plots/euler_CARTEx_composite_error")

upset_CARTEx_composite <- as.grob(upset(fromList(CARTEx_composite), order.by = "freq"))
generate_figs(upset_CARTEx_composite, "./plots/upset_CARTEx_composite", c(10, 6))


T_state <- list(
  CARTEx = CARTEx_84,
  Activation = Activation,
  Anergy = Anergy,
  Stemness = Stemness,
  Senescence = Senescence
)

euler_T_state <- euler_fitting(T_state)
generate_figs(euler_T_state$plot, "./plots/euler_T_state_plot")
generate_figs(euler_T_state$error, "./plots/euler_T_state_error")

upset_T_state <- as.grob(upset(fromList(T_state), order.by = "freq"))
generate_figs(upset_T_state, "./plots/upset_T_state", c(10, 6))


Tex_signatures <- list(
  CARTEx = CARTEx_84,
  NK_like = NK_like,
  LCMV_Tex = Wherry_Tex,
  BBD_Tex = BBD_Tex
)

euler_Tex_signatures <- euler_fitting(Tex_signatures)
generate_figs(euler_Tex_signatures$plot, "./plots/euler_Tex_signatures_plot")
generate_figs(euler_Tex_signatures$error, "./plots/euler_Tex_signatures_error")

upset_Tex_signatures <- as.grob(upset(fromList(Tex_signatures), order.by = "freq"))
generate_figs(upset_Tex_signatures, "./plots/upset_Tex_signatures", c(10, 6))


Tex_clinicalTILs <- list(
  CARTEx = CARTEx_84,
  HCC = Zheng_2017_Cell,
  Melanoma = Li_2019_Cell,
  NSCLC = Guo_2018_NatMed,
  Colorectal = Zhang_2018_Nature
)

euler_Tex_clinicalTILs <- euler_fitting(Tex_clinicalTILs)
generate_figs(euler_Tex_clinicalTILs$plot, "./plots/euler_Tex_clinicalTILs_plot")
generate_figs(euler_Tex_clinicalTILs$error, "./plots/euler_Tex_clinicalTILs_error")

upset_Tex_clinicalTILs <- as.grob(upset(fromList(Tex_clinicalTILs), order.by = "freq"))
generate_figs(upset_Tex_clinicalTILs, "./plots/upset_Tex_clinicalTILs", c(10, 6))


CARTEx_TSG <- list(
  CARTEx_630 = CARTEx_630,
  CARTEx_200 = CARTEx_200,
  CARTEx_84 = CARTEx_84,
  Tumor_suppressor_genes = TSG
)

euler_CARTEx_TSG <- euler_fitting(CARTEx_TSG)
generate_figs(euler_CARTEx_TSG$plot, "./plots/euler_CARTEx_TSG_plot")
generate_figs(euler_CARTEx_TSG$error, "./plots/euler_CARTEx_TSG_error")



