#!/usr/bin/env Rscript

# 
# QC 2: Tag cluster annotation, nucleotide and dinucleotide composition
# 


# Load libraries
required.libraries <- c(
    "optparse",
    "rlang",
    "CAGEr",
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
        c("-d", "--annot_db"),
        type = "character",
        default = NULL,
        help = "Genome annotation database package, eg org.Hs.eg.db (Mandatory)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run."),
    make_option(
        c("-t", "--tpm_threshold"),
        type = "double",
        default = 1,
        help = "Threshold to filter CTSS that has few tags per million (Optional), defaults to 1 ",
        metavar = "double"),
    make_option(
        c("-e", "--pdf_height"),
        type = "double",
        default = 5,
        help = "Height of the plot, depends on the number of samples (Optional), defaults to 5 ",
        metavar = "double"),
    make_option(
        c("-w", "--pdf_width"),
        type = "double",
        default = 4,
        help = "Width of the plot (Optional), defaults to 4 ",
        metavar = "double"),
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names

ce_path         <- opt$cageexp_object
tx_annotation   <- opt$annotation
annot_db        <- opt$annot_db
project_dir     <- opt$project_dir
tpmThreshold  <- opt$tpm_threshold
pdfWidth      <- opt$pdf_width
pdfHeight     <- opt$pdf_height

# import functions for second quality control
source(file.path(project_dir, "bin/cager_nucleotide_helper.R"))

# Read in CAGEexp object
ce <- readRDS(ce_path)

# extract tag clusters to GRanger object
sampleNames <- unname(sampleLabels(ce))
tag_clusters <- lapply(
    sampleNames,
    function (x) CAGEr::tagClustersGR(
        ce,
        sample = x,
        returnInterquantileWidth = TRUE,
        qLow = 0.1,
        qUp = 0.9))

names(tag_clusters) <- sampleNames

# annotate peaks with chipseeker
peakAnno_list <- lapply(
    tag_clusters,
    function(x) ChIPseeker::annotatePeak(
        x,
        TxDb = tx_annotation,
        tssRegion = c(-3000, 3000),
        annoDb = annot_db,
        sameStrand = TRUE)
)
pdf("chipseeker_tagCluster_annotation.pdf")
print(ChIPseeker::plotAnnoBar(peakAnno_list))
dev.off()

# nucleotide composition
normalized_ctss_list <- extract_ctss_normalized_tmp_per_sample(ce, tpmThreshold)
ctss_sequences <- extract_ctss_sequences(normalized_ctss_list, seq_data)
outlist <- calculate_nucleotide_frequency(ctss_sequences)
ctss_nucl_freq_df_tidy <- outlist[[1]]
sample_names <- outlist[[2]]
print(plot_nucleotide_frequency(
    ctss_nucl_freq_df_tidy,
    sample_names,
    "nucleotide_freq.pdf",
    pdfheight = pdfHeight,
    pdfwidth = pdfWidth
))

# dinculeotide composition
expanded_ctss_list <- expand_ctss_regions(normalized_ctss_list)
ctss_sequences <- extract_ctss_sequences(expanded_ctss_list, seq_data)
outlist <- count_dinucleotide_frequency(ctss_sequences)
ctss_dinuc_freq_df_tidy <- outlist[[1]]
dinuc_levels <- outlist[[2]]
column_names <- outlist[[3]]
print(plot_dinucleotide_frequency(
    ctss_dinuc_freq_df_tidy,
    dinuc_levels,
    column_names,
    "dinucleotide_freq.pdf",
    pdfheight = pdfHeight,
    pdfwidth = pdfWidth))