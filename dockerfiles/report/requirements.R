#################################
## List of R packages required ##
#################################

## CRAN packages:
required_packages_cran = c(
    "knitr",
    "rmarkdown",
    "ggplot2",
    "gplots",
    "ggrepel")

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

BiocManager::install(
    "Bioconductor/BiocArchive")

BiocManager::install(
    "CAGEr",
    version="3.21")
