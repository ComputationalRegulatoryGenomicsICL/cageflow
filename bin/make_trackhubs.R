#!/usr/bin/env Rscript

# 
# Script to convert bigwigs to trackhubs
# 

# Load libraries
required.libraries <- c(
  "optparse",
  "easyTrackHubs"
)

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

# parse options
option_list = list(
    make_option(
        c("-b", "--bigwigs"),
        type = "character",
        default = NULL,
        help = "Path to the bigwig files, separated by commas (Mandatory)"),
    make_option(
        c("-r", "--ref_genome"),
        type = "character",
        default = NULL,
        help = "Name of the reference genome (Mandatory)"),
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

normalized_samples <- opt$bigwigs[grep("normalized", opt$bigwigs)]
sample_names <- sapply(strsplit(normalized_samples, "_", fixed=TRUE), "[", 1)

cage_bigwig_df <- data.frame(
    sample_name = sample_names,
    file_path = normalized_samples,
    file_format = "bigWig",
    reference_genome = opt$ref_genome,
    data_type = "CAGE",
    stringsAsFactors = FALSE
)

easyTrackHub(
    cage_bigwig_df,
    trackhub_name="cageflow",
    trackhub_path="trackhubs",
    maintainer_email="none")
