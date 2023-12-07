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

# start.packages = installed.packages()[,1]
# install.packages(bsgenome, repos = NULL, type="source")
# a = setdiff(installed.packages()[,1], start.packages)
# write(a, "test.txt")

install.packages(bsgenome, repos = NULL, type="source")
# BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer1")
library("BSgenome.Scerevisiae.UCSC.sacCer1")
write(seqlengths(Scerevisiae), "test.txt")

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
write(bam.type, "single.end.uniq.txt")

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

# write.csv(BSgenome::installed.genomes(), "bsgenome.txt")

ce = CAGEexp(genomeName     = "BSgenome.Scerevisiae.UCSC.sacCer1",
              inputFiles     = input.files,
              inputFilesType = bam.type,
              sampleLabels   = sample.names)

ce = getCTSS(ce, removeFirstG = T)
# write.table(colData(ce), "ce.txt", sep='\t')
saveRDS(ce, "CAGEexp_v1_readCTSS.RDS")

# write.table(colData(ce), "ce.txt", quote=FALSE, sep='\t')

# ce <- getCTSS(ce, removeFirstG = T, useMulticore = T, nrCores = 8)


# ce2 = getCTSS(ce, removeFirstG = T)
# write.table(colData(ce2), "ce2.txt", quote=FALSE, sep='\t')

# ce3 = CTSStagCountDF(ce2)
# write.table(ce3, file='ce3.tsv', sep='\t')

# ce4 = CTSScoordinatesGR(ce2)
# write.table(ce4, file='ce4.tsv', sep='\t')

