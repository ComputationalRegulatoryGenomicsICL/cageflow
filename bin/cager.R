#!/usr/bin/env Rscript

library(CAGEr)
library(BSgenome)

args = commandArgs()
bsgenome    = args[6]
sample.list = args[7]
cpus        = args[8]

install.packages(bsgenome, repos = NULL, type="source")

path = file.path(bsgenome)
ref.name = unlist(strsplit(basename(path), "_"))[1]
ref.id = unlist(strsplit(ref.name, "\\."))[4]

library(ref.name, character.only = TRUE)

sample.table = read.delim(sample.list, header = FALSE, sep = "\t")
names(sample.table) = c("id", "single_end", "path")

sample.names = sample.table$id
single_end_uniq = unique(sample.table$single_end)
if (length(single_end_uniq) == 1) {
    bam.type = ifelse(single_end_uniq == "true",
                      "bam", "bamPairedEnd")
} else {
    stop("Sample table contains both single-end and paired-end reads.")
}

input.files = sample.table$path

ce = CAGEexp(genomeName     = ref.name,
             inputFiles     = input.files,
             inputFilesType = bam.type,
             sampleLabels   = sample.names)

ce = getCTSS(ce, removeFirstG = T, useMulticore = T, nrCores = cpus)

saveRDS(ce, paste0(ref.id, "_CAGEexp_v1_readCTSS.RDS"))
