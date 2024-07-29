#!/usr/bin/env Rscript

# 
# Script to process data with CAGEr
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
        c("-a", "--annotation"),
        type = "character",
        default = NULL,
        help = "Genome annotation package (Mandatory)"),
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
bsgenome        <- opt$bsgenome
bigwig_list     <- opt$bigwig_list
sample_list     <- opt$sample_list
num_core        <- opt$num_core
project_dir     <- opt$project_dir
tx_annotation   <- opt$annotation
range_min       <- opt$range_min
range_max       <- opt$range_max
method          <- opt$method
total_tag_num   <-opt$total_tag_num

# import functions
# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# reading in bam / bigwig data
source(file.path(project_dir, "bin/parse_input.R"))
source(file.path(project_dir, "bin/cager_bam.R"))
source(file.path(project_dir, "bin/cager_bigwig.R"))
# for quality control
source(file.path(project_dir, "bin/annotation_helper.R"))
source(file.path(project_dir, "bin/updated_plots.R"))
source(file.path(project_dir, "bin/cager_qc.R"))
# for analysis
source(file.path(project_dir, "bin/cager_normalization.R"))
source(file.path(project_dir, "bin/plot_settings.R"))
source(file.path(project_dir, "bin/cager_clustering.R"))
source(file.path(project_dir, "bin/cager_cluster_annotation.R"))

cager_folder <- "cager"

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

# Quality controls and preliminary analyses
cager_qc(
    ce=ce,
    tx_annotation=tx_annotation,
    cager_folder=cager_folder)

# Merging of replicates and Normalization
ce <- cager_normalization(
    ce=ce,
    rangeMin=range_min,
    rangeMax=range_max,
    method=method,
    total_tag_num=total_tag_num,
    cager_folder=cager_folder)

# CTSS clustering
# ce <- cager_clustering(
#     ce=ce,
#     threshold=threshold,
#     thresholdIsTpm=thresholdIsTpm,
#     nrPassThreshold=nrPassThreshold,
#     method=method,
#     maxDist=maxDist,
#     removeSingletons=removeSingletons,
#     keepSingletonsAbove=keepSingletonsAbove,
#     qLow=qLow, qUp=qUp,
#     tpmThreshold=tpmThreshold,
#     plot_lim=plot_lim,
#     num_core=num_core,
#     cagerFolder=cager_folder
# )

# QC 2: Nucleotide and dinucleotide composition and Promoter width

# Creating consensus promoters across samples

# Track export for genome browsers

# Expression profiling and Differential expression analysis

# Shifting promoters

# Enhancers

#  Unibind enrichment?