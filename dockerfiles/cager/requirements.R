#################################
## List of R packages required ##
#################################

## CRAN packages:
required_packages_cran = c(
  "optparse",         ## Read in data
  "rlang",            ## Error handling
  "tidyr",            ## Data formatting and handling
  "tidyverse",        ## Data formatting and handling
  "stringr",          ## Data formatting and handling
  "dplyr",            ## Data formatting and handling
  "purrr",            ## Parallel processing
  "magrittr",         ## Code formatting
  "viridis",          ## Plotting
  "ggplot2",          ## Plotting
  "gplots",           ## Plotting
  "ggseqlogo",        ## Plotting
  "memoise")          ## Caching

message(
  "; Installing these R packages from CRAN repository: ",
  required_packages_cran)
install.packages(
  required_packages_cran,
  repos="https://cran.uib.no/")

install.packages(
  'BiocManager',
  repos='https://cloud.r-project.org/')

BiocManager::install("remotes")

## Bioconductor packages:
required_packages_bioconductor <- c(
  "Bioconductor/BiocArchive",
  "BSgenome",
  "ChIPseeker",
  "rtracklayer",
  "txdbmaker")

message(
  "; Installing these R Bioconductor packages: ",
  required_packages_bioconductor)
BiocManager::install(
  required_packages_bioconductor)

BiocManager::install(
  "CAGEr",
  version="3.20")
