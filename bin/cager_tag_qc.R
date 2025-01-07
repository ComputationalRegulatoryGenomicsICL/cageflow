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
    "ggseqlogo",
    "memoise"
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
# import functions for quality control and plotting
source(file.path(project_dir, "bin/plot_saving.R"))
source(file.path(project_dir, "bin/cager_modified_plots.R"))

reference_name <- install_bsgenome(bsgenome)

# Read in CAGEexp object
ce <- readRDS(ce_path)

# Read in TxDb object
tx_annotation_obj <- loadDb(tx_annotation)
ce <- CAGEr::annotateCTSS(ce, tx_annotation_obj)

# Save intermediate annotated object
saveRDS(ce, "annotated_cagexp.rds")

annotations <- CAGEr::plotAnnot(ce, "counts")
save_plot(
    "tag_region_annotation_plot.pdf",
    annotations
)

if (length(sampleLabels(ce)) > 10){
    # uses function from cager_modified_plots.R
    corr_m <- plotCorrelation2_local(
        CTSStagCountDF(ce),
        samples = "all",
        tagCountThreshold = 1,
        applyThresholdBoth = FALSE,
        method = "pearson",
        digits = 3,
        plot_pairs = FALSE)

    # save intermediate file
    saveRDS(corr_m, "corr_m.rds")

    # plot correlations in heatmap format
    hm <- gplots::heatmap.2(corr_m, trace="none", margins=c(12, 12),cexRow=0.2)
    pdf("correlations_heatmap_plot.pdf")
    eval(hm$call)
    dev.off()
    saveRDS(hm, "correlations_heatmap_plot.rds")
} else {
    corr_m <- plotCorrelation2_local(
        CTSStagCountDF(ce),
        samples = "all",
        tagCountThreshold = 1,
        applyThresholdBoth = FALSE,
        method = "pearson",
        digits = 3,
        plot_pairs=TRUE)
    # save intermediate file
    saveRDS(corr_m, "corr_m.rds")
}


# Plot sequence distribution at the TSS
# tsslogo_plot <- CAGEr::TSSlogo(
#     CAGEr::CTSScoordinatesGR(ce) |> subset(annotation == "promoter"),
#     upstream = 35)
tsslogo_plot <- TSSlogo_local(
    CAGEr::CTSScoordinatesGR(ce) |> subset(annotation == "promoter"),
    genome_name=reference_name,
    upstream = 35)
save_plot(
    "TSSlogos_plot.pdf",
    tsslogo_plot
)
