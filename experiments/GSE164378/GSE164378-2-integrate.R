# integrate seurat object
# Vandon Duong

experiment = 'GSE136874'

setwd(paste("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/", experiment, sep = ''))

library(Seurat)
library(ggpubr)
library(ggplotify)
library(stringr)
library(dplyr)






# https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/articles/integration_large_datasets.html
# https://satijalab.org/seurat/articles/covid_sctmapping

