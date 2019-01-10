## Name: checkRawMatrices.R
## Description: Analysis script for the raw matrices
## Author: Nikolas Barkas
## Date: December 2018

## Load libraries
library(fastSave)
library(Matrix)

## Load helper functions
source('functions.R')

## Read in the V2 and V3 matrices
raw.V2V2 <- "/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV2_annotV2/pbmc_4k/outs/raw_gene_bc_matrices"
raw.V3V3 <- "/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV3/pbmc_4k/outs/raw_feature_bc_matrix"

v2.raw <- read10xMatrix(raw.V2V2,version='V2')
v3.raw <- read10xMatrix(raw.V3V3, version='V3')

## Check dimensions of the matrices
dim(v2.raw)
dim(v3.raw)

## Get the cells reported by V2 and V3
bc.v2 <- colnames(v2.raw)
bc.v3 <- colnames(v3.raw)

## Check number of cells reported
length(bc.v2)
length(bc.v3)

## Get number of cells with at least one read for each chemistry
table(colSums(v2.raw) >= 1)
table(colSums(v3.raw) >=1)

## Get the names of the used barcodex
usedBCv2 <- colnames(v2.raw)[colSums(v2.raw) >= 1]
usedBCv3 <- colnames(v3.raw)[colSums(v3.raw) >= 1]

## Check how many are used
length(usedBCv2)
length(usedBCv3)

## Check the overlaps of used barcodes
length(setdiff(usedBCv3,usedBCv2))
length(setdiff(usedBCv2,usedBCv3))

## Save the session
preserve.state()

