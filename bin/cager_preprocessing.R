#!/usr/bin/env Rscript

# 
# Process data with CAGEr: merging of replicates and normalization, clustering of tags
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

# parse options
option_list = list(
    make_option(
        c("-i", "--cageexp_object"),
        type = "character",
        default = NULL,
        help = "Path to the CAGEexp object with tags (Mandatory)"),
    make_option(
        c("-n", "--range_min"),
        type = "integer",
        default = 10,
        help = "CAGE tag range minimum for normalization calculation (Optional)"),
    make_option(
        c("-m", "--range_max"),
        type = "integer",
        default = 10000,
        help = "CAGE tag range maximum for normalization calculation (Optional)"),
    make_option(
        c("-e", "--method"),
        type = "character",
        default = "powerLaw",
        help = "Method for normalization (Optional)"),
    make_option(
        c("-t", "--total_tag_num"),
        type = "integer",
        default = 1*10^6,
        help = "Total number of tags. Setting it to 1 million (default) results in normalized tags per million (tpm) values (Optional)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run."),
    make_option(
        c("-c", "--num_core"),
        type = "integer",
        default = 0,
        help = "Number of cores to use (Optional), defaults to 0 (no parallelization)")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names

ce_path         <- opt$cageexp_object
range_min       <- opt$range_min
range_max       <- opt$range_max
method          <- opt$method
total_tag_num   <-opt$total_tag_num
project_dir     <- opt$project_dir
num_core        <- opt$num_core

# import functions

# for analysis
source(file.path(project_dir, "bin/cager_normalization.R"))
source(file.path(project_dir, "bin/plot_number_of_ctss.R"))
source(file.path(project_dir, "bin/cager_modified_plots.R"))
source(file.path(project_dir, "bin/cager_clustering.R"))

# Read in CAGEexp object
ce <- readRDS(ce_path)

# Merging of replicates and Normalization
# uses functions from cager_modified_plots.R
ce <- cager_normalization(
    ce=ce,
    rangeMin=range_min,
    rangeMax=range_max,
    method=method,
    total_tag_num=total_tag_num,
    cager_folder="cager_results/") # = cager_folder

# CTSS clustering
# uses functions from cager_modified_plots.R
ce <- cager_clustering(
    ce=ce,
    iqw_plot_lim=c(0, 150),
    num_core=num_core)

# save output
# RDS
# CTSS count matrix
# consensus cluster count matrix

