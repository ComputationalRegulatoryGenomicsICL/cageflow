#!/usr/bin/env bash

# Download a test dataset with paired-end reads

# Use CAGE scan data from the FANTOM5 project:
# https://fantom.gsc.riken.jp/5/datafiles/latest/basic/
# The data is publicly available and can be used for the testing of the CAGEflow pipeline.
# More information on CAGE scan can be found in Plessy et al., 2010 (https://pubmed.gov/20543846), Kratz et al., 2014 (https://pubmed.gov/24904046) and Bertin et al., 2017 (https://pubmed.gov/28972578).


mkdir CAGEscan
cd CAGEscan
mkdir human_cell_line
cd human_cell_line
wget -r -l 1 -nd -np -A fq.gz -e robots=off https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.CAGEScan/
cd ..
mkdir human_primary_cell
cd human_primary_cell
wget -r -l 1 -nd -np -A fq.gz -e robots=off https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.primary_cell.CAGEScan/
cd ..
mkdir human_tissue
cd human_tissue
wget -r -l 1 -nd -np -A fq.gz -e robots=off https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.CAGEScan/
