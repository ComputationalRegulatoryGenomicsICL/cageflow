#!/usr/bin/env Rscript

library(BSgenome)

args = commandArgs()
forge_seed = args[6]
seqs_srcdir = args[7]
genome_name = args[8]

bsgenome_destdir = paste0("BSgenome.custom.", genome_name)

dir.create(file.path(bsgenome_destdir))

forgeBSgenomeDataPkg(forge_seed, 
                     seqs_srcdir, 
                     bsgenome_destdir)
