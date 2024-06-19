#!/usr/bin/env Rscript

library(BSgenome)

args = commandArgs()
forge_seed = args[6]
seqs_srcdir = args[7]
genome_name = args[8]

forgeBSgenomeDataPkg(forge_seed, 
                     seqs_srcdir)
