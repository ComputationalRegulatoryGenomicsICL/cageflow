#!/usr/bin/env Rscript

# 
# Script to convert GTF into TxDb object
# 

# Load libraries
required.libraries <- c(
    "txdbmaker"
)

for (lib in required.libraries) {
  suppressPackageStartupMessages(library(lib, character.only=TRUE, quietly = T))
}

# parse options
option_list = list(
    make_option(
        c("-g", "--gtf_file"),
        type = "character",
        default = NULL,
        help = "Path to the GTF file")
)

message("; Reading arguments from command line.")
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)

# set variable names
gtf.file = opt$gtf_file

# create and save the TxDB object
txdb.filename = "annotation_from_gtf.sqlite"

txdb = makeTxDbFromGFF(gtf.file)

saveDb(txdb, txdb.filename)
