
## Color schemes

### Construction of CARTEx using bulk RNAseq on sorted CD8 T cells

| CAR | Color |
|:---:|:-----:|
| Control | `#darkgoldenrod` |
| CD19 | `#dodgerblue` |
| HA | `#indianred` |

### Cell cycle phase

| Phase | Color |
|:-----:|:-----:|
| G1 | `#089392` |
| S | `#EAE29C` |
| G2M | `#CF597E` |

Use this code: `hcl.colors(n = 3, palette = "Temps")`

### Monaco cell type differentiation

| Type | Color |
|:----:|:-----:|
| Naive CD8 T cells | `deepskyblue` |
| Central memory CD8 T cells | `seagreen` |
| Effector memory CD8 T cells | `darkgoldenrod` |
| Terminal effector CD8 T cells | `plum3` |

### GSE136184 aging T cells

| Age Group | Color |
|:---------:|:-----:|
| Newborn | `#D3D3D3` |
| Under 30 | `#C0D5DC` |
| Under 50 | `#ADD8E6` |
| Under 70 | `#566CD9` |
| Elderly | `#0000CD` |

Use this code: `colorRampPalette(c("lightgrey","lightblue","mediumblue"))(5)`

| Control | Color |
|:-------:|:-----:|
| YoungNaive | `#royalblue` |
| OldTerminal | `#orchid` |

### GSE136874 CD19 vs GD2 CAR T cells at day 10 of cell culture

| CAR | Color |
|:---:|:-----:|
| CD19 | `dodgerblue` |
| GD2 | `indianred` |


### GSE120575 Cancer TILs receiving checkpoint inhibitor

| Response | Color |
|:---:|:-----:|
| NR | `firebrick` |
| R | `seagreen` |







## CD8+ subtype notes

[Lückel et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32697883/)

[St Paul and Ohashi 2020](https://pubmed.ncbi.nlm.nih.gov/32624246/)

## Navigating Sherlock

[Sherlock](https://www.sherlock.stanford.edu/docs/#welcome-to-sherlock) is a high-performance computing (HPC) cluster at Stanford University. It offers interactive compute environment and storage for researchers.

RStudio Server Settings

R version: R/4.2.0

Additional modules: biology physics geos/3.6.2 hdf5/1.12.2

`sh_quota -f HOME` to examine account storage

`du -hs $(ls -A)` to examine directory storage

`~/ondemand/data/sys/dashboard/batch_connect/sys/sh_rstudio`

`~/.rstudio/ctx/ctx-8787/environmentTmp` 

`squeue -u vandon` to examine job lists

Install `devtools` on Sherlock requires:

`ml system harfbuzz fribidi`

`ml devel cmake libgit2`

`ml system openssl`

Install `python-pptx` on Sherlock requires:

`pip install --user python-pptx`

## Housekeeping

### Organizing R code
https://www.r-bloggers.com/2018/09/r-code-best-practices/

### Saving /construction/data/ files
https://stackoverflow.com/questions/1967370/git-replacing-lf-with-crlf

### Installing reference data
library(SeuratData)

InstallData("pbmcref")

### Bypassing university proxy to install SeuratData datasets and use Azimuth

https://support.posit.co/hc/en-us/articles/200488488-Configuring-R-to-Use-an-HTTP-or-HTTPS-Proxy

https://stackoverflow.com/questions/6467277/proxy-setting-for-r

https://stackoverflow.com/questions/38276457/error-in-curlcurl-fetch-memoryurl-handle-handle-failure-when-receiving






