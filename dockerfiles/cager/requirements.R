#################################
## List of R packages required ##
#################################

## CRAN packages:
required_packages_cran = c(
<<<<<<< HEAD
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
    "ggseqlogo")        ## Plotting
=======
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
>>>>>>> dev

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
<<<<<<< HEAD
    "Bioconductor/BiocArchive",
    "BSgenome",
    "ChIPseeker",
    "rtracklayer",
    "txdbmaker")
=======
  "BSgenome",
  "ChIPseeker",
  "txdbmaker")
>>>>>>> dev

message(
    "; Installing these R Bioconductor packages: ",
    required_packages_bioconductor)
BiocManager::install(
    required_packages_bioconductor)

<<<<<<< HEAD
BiocManager::install(
    "CAGEr",
    version="3.20")
=======
devtools::install_github("charles-plessy/CAGEr", ref="devel")
>>>>>>> dev
