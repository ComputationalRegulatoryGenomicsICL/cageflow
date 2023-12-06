#!/usr/bin/env Rscript

library(CAGEr)
library(BSgenome)

# assume we have 2 variables as inputs: input_folder where sorted bam files are and BSgenome
# This script also assumes single end reads, so it should be adapted to determine if samples are 
# single-end or paired-end

args = commandArgs()
# input_folder = args[6]
# arg2 = args[9]
# saveRDS(input_folder, "test.RDS")
# write(input_folder, "test.RDS")
bsgenome = args[6]

# install.packages(bsgenome, repos = NULL, type="source")
# BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer1")
# library("BSgenome.Scerevisiae.UCSC.sacCer1")
# write(seqlengths(Scerevisiae), "test.txt")
# write(packages(), "test.txt")

sample.list = args[7]
sample.table = read.delim(sample.list, header = FALSE, sep = "\t")
names(sample.table) = c("id", "single_end", "path")

sample.names = sample.table$id
read.mode = sample.table$single_end
single_end_uniq = unique(sample.table$single_end)
if (length(single_end_uniq) == 1) {
    bam.type = ifelse(single_end_uniq == "true",
                      "bam", "bamPairedEnd")
} else {
    stop("Sample table contains both single-end and paired-end reads.")
}

input.files = sample.table$path

write(sample.names, "sample.names.txt")
write(read.mode, "read.mode.txt")
write(input.files, "input.files.txt")

# inputFiles <- list.files( input_folder,
#                           ".sorted.bam$",
#                           full.names = TRUE )

# print(inputFiles)

# # I split them here by the lane label "_L00" but other delimiters could be used depending on the name of the files

# sampleNames <- unlist(lapply(strsplit(basename(inputFiles), "_L00"), function(x) x[1]))

# # I am adding some things for the BSgenome (by all mean, these are not fully tested)

# if (dir.exists(system.file(package=BSgenome))){ # if BSgenome is already installed which probably isn't
#     message("BSgenome")
# } else if (length(BiocManager::available(BSgenome)) > 0){
#   message(paste0(BSgenome, " available in Bioconductor. Proceeding to install..."))
#   BiocManager::install(BSgenome)
# } else { # if the genome is not available in Bioconductor, it can either be installed from the file or forged
#  # check if install.packages(BSgenome) successful. If not, forge the BSgenome packages using BSgenome function forgeBSgenomeDataPkg
#  # it requires the BSgenome-seed file 
# }

# ce <- CAGEexp(genomeName     = BSgenome,
#               inputFiles     = input.files,
#               inputFilesType = "bam",
#               sampleLabels   = sample.names)

# ce <- getCTSS(ce, removeFirstG = T, useMulticore = T, nrCores = 8)

# saveRDS(ce, "/whichever/output/folder/CAGEexp_v1_readCTSS.RDS")