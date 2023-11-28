# CARTEx

- [About](#about)
- [Experiments](#experiments)
  - [Datasets](#datasets)
  - [Utilities](#utilities)
- [Literature](#literature)
  - [Reviews](#reviews)
  - [Signatures](#signatures)
- 

## About <a name = "about"></a>

This is the github repository for the CARTEx signature project. Pre-processing and analyses were performed using [Seurat](https://satijalab.org/seurat/) and custom pipelines written in [R](https://www.r-project.org/). The purpose of this research was to develop a robust signature and procedure for quantifying T cell exhaustion states. 

Historically, T cell exhaustion was described based on a few canonical markers. Recent advances in single-cell sequencing has enabled high-resolution of transcriptomics and deeper insights underpinning cellular behavior, but many groups still identified "exhausted" T cells based on simple correlations with canonical markers (e.g. "PDCD1 is highly expressed in this cluster of cells, so these may be exhausted"). Using a robust model of CAR T cell exhaustion, we created a transcriptional signature to quantify T cell exhaustion. This signature was derived from differentially expressed genes between highly functional CD19 CAR T cells, exhaustion-prone HA CAR T cells, and control T cells. The HA CAR spontaneously aggregates and causes tonic signaling, which has been validated as a robust model of exhaustion in numerous studies ([Long et al. 2015 Nature Medicine](https://pubmed.ncbi.nlm.nih.gov/25939063/), [Lynn et al. 2019 Nature](https://pubmed.ncbi.nlm.nih.gov/31802004/), [Gennert et al. 2021 PNAS](https://pubmed.ncbi.nlm.nih.gov/34285077/), [Weber et al. 2021 Science](https://pubmed.ncbi.nlm.nih.gov/33795428/)).

Preprint: coming soon...

Publication: coming soon...

## Experiments <a name = "experiments"></a>

We derived the CARTEx signature from in-house experiments comparing highly functional CD19 CAR T cells, exhaustion-prone HA CAR T cells, and control T cells. We then applied the CARTEx signature to both in-house and publicly available experiments.

General procedure:

1. Prepare datasets and filter cells based on
2. Calculate signatures
3. Visualize dataset and organize based on metadata

### Datasets <a name = "datasets"></a>

| Identifier | Experiment | Study |
|:----------:|------------|-------|
| [GSE136874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136874) | Functional CD19 CAR T cells vs exhaustion-prone GD2 CAR T cells | [Lynn et al. 2019 Nature](https://pubmed.ncbi.nlm.nih.gov/31802004/) |
| [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) | TILs from cancer patients after receiving immune checkpoint blockade | [Sade-Feldman et al. 2018 Cell](https://pubmed.ncbi.nlm.nih.gov/30388456/) |
| [GSE235676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235676) | TILs from cancer patients before receiving immune checkpoint blockade | [Mei et al. 2023 Nature Cancer](https://pubmed.ncbi.nlm.nih.gov/37460871/) |
| [GSE125881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125881) | CD19 CAR T cells from cancer patients across 16 weeks | [Sheih et al. 2020 Nature Communications](https://pubmed.ncbi.nlm.nih.gov/31924795/) |
| [GSE151511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151511) | CD19 CAR T cell infusion product from cancer patients | [Deng et al. 2020 Nature Medicine](https://pubmed.ncbi.nlm.nih.gov/33020644/) |
| [GSE153931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153931) | SARS-CoV-2 reactive T cells from COVID-19 patients | [Kusnadi et al. 2021 Science Immunology](https://pubmed.ncbi.nlm.nih.gov/33478949/) |
| [GSE136184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136184) | Aging CD8 T cells from healthy donors across 9 decades | [Lu et al. 2022 Nature Communications](https://pubmed.ncbi.nlm.nih.gov/36050300/) |
| [GSE207935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207935) | T cells from patients with STAT3 gain-of-function syndrome | [Schmitt et al. 2022 JCI Insight](https://pubmed.ncbi.nlm.nih.gov/36136607/) |
| [GSE212217](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212217) | T cells from cancer patients with mismatch-repair deficiency | [Chow et al. 2023 Cancer Discovery](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9905265/) |
| [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378) | Azimuth reference atlas comprised of human peripheral blood mononuclear cells (PBMCs) | [Hao et al. 2021 Cell](https://pubmed.ncbi.nlm.nih.gov/34062119/) |
| [GSE107011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011) | Monaco refererence atlas | [Monaco et al. 2019 Cell Reports](https://pubmed.ncbi.nlm.nih.gov/30726743/) |
| [Zenodo 3993994](https://zenodo.org/record/3993994) | T cells from patients with Parkinson's disease | [Wang et al. 2021 Cell Discovery](https://www.nature.com/articles/s41421-021-00280-3) |

Note, we do not include the raw or processed data for these experiments, as the file sizes are far too large. However, our scripts for preparing and exploring the datasets are shared in this repository. Publicly-available experiments can be downloaded from the Gene Expression Omnibus (GEO).




[More notes](REF.md)


### Utilities <a name = "utilities"></a>

We built custom functions and archived these within [`cartex-utilities.R`](cartex-utilities.R). Some of the key functions included:
- `generate_figs()` to save figures in jpeg and pdf formats
- `integerize()` to round scores to nearest integer
- `DimPlotHighlightIdents()` to customize Seurat's [`DimPlot()`](https://satijalab.org/seurat/reference/dimplot) for highlighting cells by identity group

## Literature <a name = "literature"></a>

### Reviews <a name = "reviews"></a>

[Pauken and Wherry 2015 Trends in Immunology](https://pubmed.ncbi.nlm.nih.gov/25797516/) discusses the involvement of T cell exhaustion in infection and cancer.

[Collier et al. 2021 Nature Immunology](https://pubmed.ncbi.nlm.nih.gov/34140679/) discusses the involvement of T cell exhaustion in autoimmune diseases.

### Signatures <a name = "signatures"></a>

[Wherry et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17950003/) describes CD8 T cell exhaustion based on chronic infection by lymphocytic choriomeningitis virus (LCMV).

[Good et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34861191/) describes NK-like CAR T cell exhaustion based on continuous antigen exposure of mesothelin CAR T cells.







