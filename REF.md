


## Navigating Sherlock

[Sherlock](https://www.sherlock.stanford.edu/docs/#welcome-to-sherlock) is a high-performance computing (HPC) cluster at Stanford University. It offers interactive compute environment and storage for researchers.

RStudio Server Settings

R version: R/4.2.0

Additional modules: biology physics geos/3.6.2 hdf5/1.12.2

`sh_quota -f HOME` to examine account storage

`du -hs $(ls -A)` to examine directory storage

`~/.rstudio/ctx/ctx-8787/environmentTmp` 

`squeue -u vandon` to examine job lists

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

