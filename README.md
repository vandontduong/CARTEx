# CARTEx

## About

This is the github repository for the CARTEx signature project. Pre-processing and analyses were performed using [Seurat](https://satijalab.org/seurat/) and custom pipelines written in [R](https://www.r-project.org/).



## Datasets

We derived the CARTEx signature from in-house experiments comparing highly functional CD19 CAR T cells, exhaustion-prone HA CAR T cells, and control T cells. We then applied the CARTEx signature to both in-house and publicly available experiments.

| Identifier | Experiment | 
|:----------:|------------|
| GSE136874  | Functional CD19 CAR T cells vs exhaustion-prone GD2 CAR T cells |
| GSE120575  | TILs from cancer patients after receiving immune checkpoint blockade |
| GSE235676  | TILs from cancer patients before receiving immune checkpoint blockade |
| GSE125881  | CD19 CAR T cells from cancer patients across 16 weeks |
| GSE151511  | CD19 CAR T cell infusion product from cancer patients |
| GSE153931  | SARS-CoV-2 reactive T cells from COVID-19 patients |
| GSE136184  | Aging CD8 T cells from healthy donors across 9 decades |


Note, we do not include the raw or processed data for these experiments, as the file sizes are far too large. However, our scripts for preparing and exploring the datasets are shared in this repository. Publicly-available experiments can be downloaded from the Gene Expression Omnibus (GEO).




[More notes](REF.md)

