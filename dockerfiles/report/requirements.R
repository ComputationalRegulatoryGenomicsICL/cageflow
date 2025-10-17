#################################
## List of R packages required ##
#################################

## CRAN packages:
required_packages_cran = c(
  "knitr",
  "rmarkdown",
  "ggplot2",
  "gplots",
  "ggrepel",
  "devtools")

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

devtools::install_github("charles-plessy/CAGEr", ref="devel")
