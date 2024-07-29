#!/usr/bin/env Rscript

# 
# Initial quality control of CAGE reads
# 


# Load libraries
required.libraries <- c(
    "optparse",
    "CAGEr",
    "gplots"
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
        help = "Genome annotation package, eg TxDb.Hsapiens.UCSC.hg38.knownGene (Mandatory)"),
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
project_dir     <- opt$project_dir

# import functions for quality control
source(file.path(project_dir, "bin/annotation_helper.R"))
source(file.path(project_dir, "bin/updated_plots.R"))

# Read in CAGEexp object
ce <- readRDS(ce_path)

ce <- annotate_gene_regions(ce, tx_annotation)

pdf("tag_region_annotation.pdf")
annotations <- CAGEr::plotAnnot(ce, "counts")
print(annotations)
dev.off()

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
