#!/usr/bin/env Rscript

# 
# QC 2: Tag cluster annotation, nucleotide and dinucleotide composition
# 


# Load libraries
required.libraries <- c(
    "optparse",
    "rlang",
    "CAGEr",
    "GenomicFeatures",
    "ChIPseeker",
    "Biostrings",
    "tidyr",
    "viridis",
    "tidyverse",
    "ggplot2"
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
        help = "Path to the CAGEexp object with tag clusters (Mandatory)"),
    make_option(
        c("-a", "--annotation"),
        type = "character",
        default = NULL,
        help = "Genome annotation package, eg TxDb.Hsapiens.UCSC.hg38.knownGene (Mandatory)"),
    make_option(
        c("-b", "--bsgenome"),
        type = "character",
        default = NULL,
        help = "Name of the BSgenome version to be used (Mandatory)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run."),
    make_option(
        c("-t", "--tpm_threshold"),
        type = "double",
        default = 1,
        help = "Threshold to filter CTSS that has few tags per million (Optional), defaults to 1 "),
    make_option(
        c("-k", "--pca_rank"),
        type = "integer",
        default = 50,
        help = "Rank of PCAs for analysis. (Default = 50) ")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names

ce_path         <- opt$cageexp_object
tx_annotation   <- opt$annotation
bsgenome        <- opt$bsgenome
project_dir     <- opt$project_dir
tpmThreshold    <- opt$tpm_threshold
pca_rank        <- opt$pca_rank

# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# Read in TxDb object
tx_annotation_obj <- loadDb(tx_annotation)
# import functions for second quality control and plotting
source(file.path(project_dir, "bin/plot_saving.R"))
source(file.path(project_dir, "bin/cager_nucleotide_composition_functions.R"))
source(file.path(project_dir, "bin/cager_consensus_qc.R"))
source(file.path(project_dir, "bin/plot_number_and_pca_of_ctss.R"))


reference_name <- install_bsgenome(bsgenome)

# Read in CAGEexp object
ce <- readRDS(ce_path)

# extract tag clusters to GRanger object
sampleNames <- unname(sampleLabels(ce))
tag_clusters <- lapply(
    sampleNames,
    function (x) CAGEr::tagClustersGR(
        ce,
        sample = x,
        qLow = 0.1,
        qUp = 0.9))

names(tag_clusters) <- sampleNames

# annotate peaks with chipseeker
peakAnno_list <- lapply(
    tag_clusters,
    function(x) ChIPseeker::annotatePeak(
        x,
        TxDb = tx_annotation_obj,
        tssRegion = c(-3000, 3000),
        sameStrand = TRUE)
)
chipannot_plot <- ChIPseeker::plotAnnoBar(peakAnno_list)
save_plot(
    "chipseeker_tagCluster_annotation_plot.pdf",
    chipannot_plot
)

# # nucleotide composition
normalized_ctss_list <- extract_ctss_normalized_tmp_per_sample(ce, tpmThreshold)
ctss_sequences <- extract_ctss_sequences(normalized_ctss_list, reference_name)
outlist <- calculate_nucleotide_frequency(ctss_sequences)
ctss_nucl_freq_df_tidy <- outlist[[1]]
sample_names <- outlist[[2]]
nuclfreq_plot <- plot_nucleotide_frequency(
    ctss_nucl_freq_df_tidy,
    sample_names
)
save_plot(
    "nucleotide_frequencies_plot.pdf",
    nuclfreq_plot
)

# dinculeotide composition
expanded_ctss_list <- expand_ctss_regions(normalized_ctss_list, reference_name)
ctss_sequences <- extract_ctss_sequences(expanded_ctss_list, reference_name)
ctss_dinuc_freq_df_tidy <- count_dinucleotide_frequency(ctss_sequences)
dinuclfreq_plot <- plot_dinucleotide_frequency(
    ctss_dinuc_freq_df_tidy)
save_plot(
    "dinucleotide_frequencies_plot.pdf",
    dinuclfreq_plot
)

# Consensus clustered CTSS quality plots
# uses functions from plot_number_and_pca_of_ctss.R
# consensus_qc(ce=ce, pcarank=pca_rank)
