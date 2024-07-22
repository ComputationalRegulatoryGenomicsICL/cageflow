#!/usr/bin/env Rscript

# 
# File to process data with CAGEr
# 

# Load libraries
required.libraries <- c(
    "optparse",
    "rlang",
    "CAGEr",
    "BSgenome",
    "dplyr",
    "purrr",
    "magrittr",
    "stringr",
    "tidyr"
    )

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

# parse options
option_list = list(
    make_option(
        c("-b", "--bsgenome"),
        type = "character",
        default = NULL,
        help = "Name of the BSgenome version to be used (Mandatory)"),
    make_option(
        c("-w", "--bigwig_list"),
        type = "complex",
        default = NULL,
        help = " List of bigwig files (Optional)"),
    make_option(
        c("-s", "--sample_list"),
        type = "complex",
        default = NULL,
        help = "Tsv file with information from the input channel with [id, pairedness, path] (Optional)"),
    make_option(
        c("-c", "--num_core"),
        type = "integer",
        default = 0,
        help = "Number of cores to use (Optional), defaults to 0 (no parallelization)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run.")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names
bsgenome    <- opt$bsgenome
bigwig_list <- opt$bigwig_list
sample_list <- opt$sample_list
num_core    <- opt$num_core
project_dir <- opt$project_dir

# import functions
source(file.path(project_dir, "bin/install_bsgenome.R"))
source(file.path(project_dir, "bin/parse_input.R"))
source(file.path(project_dir, "bin/cager_bam.R"))
source(file.path(project_dir, "bin/cager_bigwig.R"))

reference_name <- install_bsgenome(bsgenome)
reference_id <- unlist(strsplit(reference_name, "\\."))[4]

if (length(sample_list) > 0) {
    sample_table <- parse_input(sample_list)
    single_end_uniq <- unique(sample_table$single_end)
    if (length(single_end_uniq) == 1) {
        bam_type <- ifelse(single_end_uniq == "true",
                        "bam", "bamPairedEnd")
    } else {
        stop("Sample table contains both single-end and paired-end reads.")
    }
    ce <- read_in_bam(
        bsgenome_name=reference_name,
        input_files=sample_table$path,
        bam_pairedness=bam_type,
        sample_names=sample_table$id,
        cpus=num_core
    )
}else if(length(bigwig_list) > 0) {
    ce <- read_in_bigwig(
        bsgenome_name=reference_name,
        bigwig_str=bigwig_list,
        cpus=num_core
    )
} else {
    stop("Either bigwig or bam files should be provided")
}

saveRDS(ce, paste0(reference_id, "_CAGEexp_CTSS.rds"))



