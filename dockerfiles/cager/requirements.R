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
  "gplots",           ## Plotting
  "ggseqlogo",        ## Plotting
  "devtools")         ## Install other packages from github

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

remotes::install_version(
  "ggplot2",
  version = "3.4.4",
  repos = "https://cloud.r-project.org/")

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES="never")

BiocManager::install("rtracklayer", ask = FALSE)

## Bioconductor packages:
required_packages_bioconductor <- c(
  "BSgenome",
  "ChIPseeker",
  "txdbmaker")

message(
  "; Installing these R Bioconductor packages: ",
  required_packages_bioconductor)
BiocManager::install(
  required_packages_bioconductor)

devtools::install_github("charles-plessy/CAGEr", ref="devel")
