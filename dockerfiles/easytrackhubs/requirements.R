#################################
## List of R packages required ##
#################################

## CRAN packages:
required_packages_cran = c(
  "devtools")        ## Installing from github

message(
  "; Installing these R packages from CRAN repository: ",
  required_packages_cran)
install.packages(
  required_packages_cran,
  repos="https://cran.uib.no/")

devtools::install_github("Novartis/easyTrackHubs")


