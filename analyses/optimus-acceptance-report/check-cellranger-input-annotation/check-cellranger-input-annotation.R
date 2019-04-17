## This script is for checking the gtf input to cellranger and finding how many
## distinct genes are present in it (after the filtration that cellranger imposes).
## This is used to confirm that the reason that genes are missing from the cellranger
## output is the input filtering step.

## This is meant to be run interactively

## Installation of packages
## install.packages('data.table')
## source("https://bioconductor.org/biocLite.R")
## biocLite("rtracklayer")

## Load libraries
library(rtracklayer)

## Path for gtf input to cellranger
gtf.path <- '/home/nbarkas/disk3/cellranger/GRCh38/genes/genes.gtf'

## Load the gtf
gtf <- rtracklayer::import(gtf.path)

## Check  and extract unique genes
head(gtf)
genes <- unique(sort(as.character(gtf$gene_name)))

## Count unique genes
length(genes)
