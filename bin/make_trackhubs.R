#!/usr/bin/env Rscript

# 
# Script to convert bigwigs to trackhubs
# 

# Load libraries
required.libraries <- c(
  "optparse",
  "stringr",
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
        help = "Path to the bigwig files, separated by whitespace (Mandatory)"),
    make_option(
        c("-r", "--ref_genome"),
        type = "character",
        default = NULL,
        help = "Name of the reference genome (Mandatory)")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

bigwigs <- unlist(strsplit(opt$bigwigs, " "))
normalized_samples <- bigwigs[grep("normalized", bigwigs)]
sample_names <- sapply(strsplit(normalized_samples, "_", fixed=TRUE), "[", 1)
stranded_names <- sapply(strsplit(normalized_samples, ".", fixed=TRUE), "[", 1)
color_names <- ifelse(stringr::str_detect(stranded_names,"plus"),"red","blue")

cage_bigwig_df <- data.frame(
    sample_name = sample_names,
    file_path = normalized_samples,
    file_format = "bigWig",
    reference_genome = opt$ref_genome,
    shortLabel=stranded_names,
    longLabel=stranded_names,
    color=color_names,
    data_type="CAGE",
    stringsAsFactors = FALSE
)

easyTrackHubs::easyTrackHub(
    cage_bigwig_df,
    trackhub_name="cageflow",
    trackhub_path="trackhubs",
    maintainer_email="none")
