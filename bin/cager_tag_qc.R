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
    "ggseqlogo"
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
        c("-p", "--project_dir"),
        type = "character",
        default = 0,
        help = "Project directory, from which the analysis is run.")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names
ce_path         <- opt$cageexp_object
tx_annotation   <- opt$annotation
bsgenome        <- opt$bsgenome
project_dir     <- opt$project_dir

# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# import functions for quality control
source(file.path(project_dir, "bin/annotation_from_txdb_functions.R"))
source(file.path(project_dir, "bin/cager_modified_plots.R"))

reference_name <- install_bsgenome(bsgenome)

# Read in CAGEexp object
ce <- readRDS(ce_path)

# Read in TxDb object
tx_annotation_obj <- loadDb(tx_annotation)
ce <- annotate_gene_regions(ce, tx_annotation_obj, debugMode = F) # tx_annotation

pdf("tag_region_annotation.pdf")
annotations <- CAGEr::plotAnnot(ce, "counts")
print(annotations)
dev.off()

# uses function from cager_modified_plots.R
corr_m <- plotCorrelation2_local(
    CTSStagCountDF(ce),
    samples = "all",
    tagCountThreshold = 1,
    applyThresholdBoth = FALSE,
    method = "pearson",
    digits = 3,
    toPlot = FALSE)

# save intermediate file
saveRDS(corr_m, "corr_m.rds")

# plot correlations in heatmap format
pdf("correlations_heatmap.pdf")
gplots::heatmap.2(corr_m, trace="none", margins=c(12, 12),cexRow=0.2)
dev.off()

# Plot sequence distribution at the TSS
pdf("TSSLogos.pdf")
CAGEr::TSSlogo(
    CAGEr::CTSScoordinatesGR(ce) |> subset(annotation == "promoter"),
    upstream = 35)
dev.off()
