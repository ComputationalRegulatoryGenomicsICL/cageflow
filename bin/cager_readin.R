#!/usr/bin/env Rscript


# 
# Script to read in data to CAGEr from BAM or BigWig format
# 

# Load libraries
required.libraries <- c(
    "optparse",
    "BSgenome",
    "CAGEr",
    "stringr",
    "purrr",
    "dplyr",
    "tidyr",
    "magrittr"
)

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

# parse options
option_list = list(
    make_option(
        c("-t", "--data_type"),
        type = "character",
        default = NULL,
        help = "Whether BAM (bam) or BigWig (bigwig) is provided (Mandatory)"),
    make_option(
        c("-s", "--sample_table_list"),
        type = "character",
        default = NULL,
        help = "Csv with information from the input channel with [id, pairedness, bigwig or bam path] (Mandatory)"),
    make_option(
        c("-b", "--bsgenome"),
        type = "character",
        default = NULL,
        help = "Name of the BSgenome version to be used (Mandatory)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = NULL,
        help = "Project directory, from which the analysis is run. (Mandatory)"),
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
data_type           <- opt$data_type
sample_table_list   <- opt$sample_table_list
bsgenome            <- opt$bsgenome
project_dir         <- opt$project_dir
num_core            <- opt$num_core

# import functions

# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# reading in bam / bigwig data
source(file.path(project_dir, "bin/parse_input.R"))
source(file.path(project_dir, "bin/cager_bam.R"))
source(file.path(project_dir, "bin/cager_bigwig.R"))

# Create folders for organized analysis
dir.create(file.path("plots"))
dir.create(file.path("tracks"))
dir.create(file.path("tables"))
dir.create(file.path("intermediate_cagerobj"))

reference_name <- install_bsgenome(bsgenome)

sample_table <- parse_input(sample_table_list, data_type)
single_end_uniq <- unique(sample_table$single_end)
if (length(single_end_uniq) < 1) {
    print(sample_table)
    stop("Sample table is empty or the header is missing.")
} else if (length(single_end_uniq) > 1) {
    print(sample_table)
    stop("Sample table contains both single-end and paired-end reads.")
} else {
    bam_type <- ifelse(
        single_end_uniq == "true",
        "bam",
        "bamPairedEnd")
}

if (tolower(data_type) == "bam"){
    ce <- read_in_bam(
        bsgenome_name=reference_name,
        bam_paths=sample_table$path,
        bam_pairedness=bam_type,
        sample_names=sample_table$id,
        cpus=num_core
    )
}else if(tolower(data_type) == "bigwig") {
    sample_names_files_dict <- list()
    for(idx in 1:nrow(sample_table)) {
        row <- sample_table[idx,]
        path1 <- trimws(unlist(strsplit(row$path, ","))[1])
        path2 <- trimws(unlist(strsplit(row$path, ","))[2])
        if (grepl("str1", path1)){
            sample_names_files_dict[[path1]] <- paste0(trimws(row$id), "_str1")
            sample_names_files_dict[[path2]] <- paste0(trimws(row$id), "_str2")
        } else if (grepl("str2", path1)){
            sample_names_files_dict[[path1]] <- paste0(trimws(row$id), "_str2")
            sample_names_files_dict[[path2]] <- paste0(trimws(row$id), "_str1")
        } else {
            print(path1)
            print(path2)
            stop("Unexpected path names")
        }
        
    }
    ce <- read_in_bigwig(
        bsgenome_name=reference_name,
        bigwig_paths=sample_table$path,
        sample_names=sample_names_files_dict
    )
} else {
    stop("Either bigwig or bam files should be provided")
}

# save the initial CAGEexp object
saveRDS(ce, "intermediate_cagerobj/initial_cagexp.rds")
