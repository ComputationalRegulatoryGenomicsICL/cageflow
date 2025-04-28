#!/usr/bin/env Rscript

# #
# # Call enhancers with CAGEfightR
# #

# # Load libraries
required.libraries <- c(
    "optparse",
    "CAGEfightR",
    "GenomicRanges"
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
        help = "Path to the CAGEexp object with normalized tags and consensus clusters (Mandatory)"),
    make_option(
        c("-x", "--annotation"),
        type = "character",
        default = NULL,
        help = "SQLite file with a TxDb genome annotation package (Mandatory)"),
    make_option(
        c("-b", "--cfBalanceThreshold"),
        type = "double",
        default = 0.95,
        help = "threshold for the cagefightr balance score (Default=0.95)")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names
ce_path             <- opt$cageexp_object
tx_annotation       <- opt$annotation
cfBalanceThreshold   <- opt$cfBalanceThreshold

# import functions
# installing BSgenome
source(file.path(project_dir, "bin/install_bsgenome.R"))
# for analysis
source(file.path(project_dir, "bin/enhancer_functions.R"))
source(file.path(project_dir, "bin/qc_plots.R"))
source(file.path(project_dir, "bin/plot_saving.R"))

# Create folders for organized analysis
dir.create(file.path("plots"))
dir.create(file.path("tracks"))
dir.create(file.path("tables"))
dir.create(file.path("intermediate_cagerobj"))

# Read in CAGEexp object
ce <- readRDS(ce_path)

# call enhancers with CAGEfightR
supported_enhancers <- cagefightr_enhancers(
    ce=ce,
    cfBalanceThreshold=cfBalanceThreshold)

# exclude enhancers overlapping promoters defined by consensus clusters
true_enhancers <- exclude_enhancers_overlapping_promoters(
    BCs=supported_enhancers,
    ce=ce)

# annotate enhancers with transcript database information
annotated_enhancers <- annotate_enhancers(
    enhancers=true_enhancers,
    txdb=tx_annotation)

saveRDS(annotated_enhancers, file = "intermediate_cagerobj/enhancers.rds")
return("Annotated enhancers rds file saved")

# assign enhancers to samples
enhancer_expr_per_sample <- identify_sample_specific_enhancers(
    true_enhancers=true_enhancers,
    ce=ce)
# save enhancer expression per sample
outFileNameSamples <- file.path(
    "tables", "enhancer_expression_per_sample.tsv")
write.table(
    enhancer_expr_per_sample,
    file=outFileNameSamples,
    quote=FALSE,
    sep='\t')
print("Enhancer expressions per sample saved to file")

# plot PCA of enhancer expression per sample
pca_plot <- plot_pcs(
    count_matrix=enhancer_expr_per_sample)
save_plot(
    "enhancer_expression_pca_plot.pdf",
    pca_plot)
print("PCA plot of enhancer expression per sample saved")

# count and plot number of enhancers per sample
sample_enhancer_count <- count_number_of_enhancers(
    enhancer_expr_per_sample=enhancer_expr_per_sample)
# function from qc_plots.R
enhancer_count_plot <- plot_number_of_tag_clusters(
    sample_tag_count=sample_enhancer_count,
    yaxistitle="Number of enhancers per sample",
    mytitle="Number of enhancers per sample",
    myfilename="enhancer_count_per_sample")
save_plot(
    "enhancer_count_per_sample.pdf",
    enhancer_count_plot)
print("Enhancer counts plotted")


