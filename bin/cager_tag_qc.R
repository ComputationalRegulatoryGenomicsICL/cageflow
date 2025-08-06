#!/usr/bin/env Rscript

#
# Initial quality control of CAGE reads
#


# Load libraries
required.libraries <- c(
    "optparse",
    "CAGEr",
    "GenomicFeatures",
    "gplots",
    "ggplot2",
    "ggseqlogo")

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
        c("-a", "--annotation"),
        type = "character",
        default = NULL,
        help = "SQLite file with a TxDb genome annotation package (Mandatory)"),
    make_option(
        c("-b", "--bsgenome"),
        type = "character",
        default = NULL,
        help = "Name of the BSgenome version to be used (Mandatory)"),
    make_option(
        c("-t", "--corrplot_tagCountThreshold"),
        type = "integer",
        default = 1,
        help = "Threshold for considering tags when calculating correlations (Default = 1)"),
    make_option(
        c("-p", "--project_dir"),
        type = "character",
        default = NULL,
        help = "Project directory, from which the analysis is run.")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names
ce_path         <- opt$cageexp_object
tx_annotation   <- opt$annotation
bsgenome        <- opt$bsgenome
corrplot_tagCountThreshold <- opt$corrplot_tagCountThreshold
project_dir     <- opt$project_dir

# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# import functions for quality control and plotting
source(file.path(project_dir, "bin/plot_saving.R"))
source(file.path(project_dir, "bin/cager_modified_plots.R"))
source(file.path(project_dir, "bin/qc_plots.R"))

reference_name <- install_bsgenome(bsgenome)

# Create folders for organized analysis
dir.create(file.path("plots"))
dir.create(file.path("tracks"))
dir.create(file.path("tables"))
dir.create(file.path("intermediate_cagerobj"))

# Read in CAGEexp object
ce <- readRDS(ce_path)

# Read in TxDb object
tx_annotation_obj <- loadDb(tx_annotation)
ce <- CAGEr::annotateCTSS(ce, tx_annotation_obj)

# Save intermediate annotated object
saveRDS(ce, "intermediate_cagerobj/annotated_cagexp.rds")

annotations <- CAGEr::plotAnnot(ce, "counts")
save_plot(
    "tag_region_annotation_plot.pdf",
    annotations
)

# to compare raw counts CTSStagCountDF is used
# bypassing the automatic selection of this assay
# uses function from cager_modified_plots.R and qc_plots.R

plot_correlation(
    datatype="CTSS",
    dataframe=CTSStagCountDF(ce),
    corrplot_tagCountThreshold=corrplot_tagCountThreshold)
print("CTSS correlation plotted")
