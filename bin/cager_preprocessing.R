#!/usr/bin/env Rscript

# 
# Process data with CAGEr: normalization, clustering of tags, consensus cluster calling, and track export
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
        c("-s", "--sample_num_thr"),
        type = "integer",
        default = 1,
        help = "Number of samples in which the CTSS Tpm threshold (ctss_thr) should be passed (Default = 1)"),
    make_option(
        c("-r", "--ctss_thr"),
        type = "integer",
        default = 1,
        help = "CTSS Tpm threshold which should be passed in sample_num_thr number of samples (Default = 1)"),
    make_option(
        c("-u", "--consensus_ctss_thr"),
        type = "integer",
        default = 5,
        help = "CTSS Tpm threshold for consensus clustering (Default = 5)"),
    make_option(
        c("-d", "--consensus_ctss_dist"),
        type = "integer",
        default = 100,
        help = "Distance threshold for consensus clustering (Default = 100)"),
    make_option(
        c("-a", "--annotation"),
        type = "character",
        default = NULL,
        help = "SQLite file with a TxDb genome annotation package (Mandatory)"),
    make_option(
        c("-k", "--pca_rank"),
        type = "integer",
        default = 50,
        help = "Rank of PCAs for analysis. (Default = 50) "),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run."),
    make_option(
        c("-b", "--bsgenome"),
        type = "character",
        default = NULL,
        help = "Name of the BSgenome version to be used (Mandatory)"),
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

ce_path             <- opt$cageexp_object
range_min           <- opt$range_min
range_max           <- opt$range_max
method              <- opt$method
total_tag_num       <- opt$total_tag_num
sample_num_thr      <- opt$sample_num_thr
ctss_thr            <- opt$ctss_thr
consensus_ctss_thr  <- opt$consensus_ctss_thr
consensus_ctss_dist <- opt$consensus_ctss_dist
tx_annotation       <- opt$annotation
pca_rank             <- opt$pca_rank
project_dir         <- opt$project_dir
bsgenome            <- opt$bsgenome
num_core            <- opt$num_core

# import functions
# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))

# for analysis
# source(file.path(project_dir, "bin/cager_merge_replicates.R"))
source(file.path(project_dir, "bin/cager_normalization.R"))
source(file.path(project_dir, "bin/plot_saving.R"))
source(file.path(project_dir, "bin/plot_number_and_pca_of_ctss.R"))
source(file.path(project_dir, "bin/cager_modified_plots.R"))
source(file.path(project_dir, "bin/cager_clustering.R"))
source(file.path(project_dir, "bin/cager_consensus_clustering.R"))
source(file.path(project_dir, "bin/cager_consensus_qc.R"))

reference_name <- install_bsgenome(bsgenome)

# Read in CAGEexp object
ce <- readRDS(ce_path)

# Normalization
# uses functions from cager_modified_plots.R
ce <- cager_normalization(
    ce=ce,
    rangeMin=range_min,
    rangeMax=range_max,
    method=method,
    total_tag_num=total_tag_num)

# CTSS clustering
# uses functions from cager_modified_plots.R and plot_number_of_ctss.R
ce <- cager_clustering(
    ce=ce,
    iqw_plot_lim=c(0, 150),
    sample_num_thr=sample_num_thr,
    ctss_thr=ctss_thr,
    num_core=num_core)

# save output
# RDS
saveRDS(ce, file = "normalized_clustered_cagexp.rds")

# Consensus clustering of clustered CTSS

ce <- consensus_clustering(
    ce=ce,
    tpmThreshold=consensus_ctss_thr,
    maxDist=consensus_ctss_dist,
    tx_annotation=tx_annotation,
    num_core=num_core)

# Consensus clustered CTSS quality plots
# uses functions from plot_number_and_pca_of_ctss.R
consensus_qc(ce=ce, pcarank=pca_rank)


