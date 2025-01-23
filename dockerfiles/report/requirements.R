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
  "memoise")

message(
  "; Installing these R packages from CRAN repository: ",
  required_packages_cran)
install.packages(
  required_packages_cran,
  repos="https://cran.uib.no/")

