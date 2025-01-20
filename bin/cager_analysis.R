#!/usr/bin/env Rscript

# 
# analyse CTSS with CAGEr: expression profiling, differential expression analysis, shifting promoters, and enhancer calling
# 

# Load libraries
required.libraries <- c(
    "optparse",
    "rlang",
    "CAGEr",
    "dplyr",
    "purrr",
    "magrittr",
    "stringr",
    "tidyr",
    "tibble",
    "data.table"
    )

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}
