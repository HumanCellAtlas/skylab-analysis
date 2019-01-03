#!/bin/bash

## File: test_emptyDropsWrapper.sh
## Author: Nikolas Barkas
## Date: January 3, 2019
## Description: Unit test for emptyDropsWrapper script

## Download some sample data from 10X
printf "Downloading pbmc4k data..."
wget -q http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
printf "done\n"

## Decompress dataset
printf "Decompressing archive..."
tar xzf pbmc4k_raw_gene_bc_matrices.tar.gz
printf "done\n"

## Prepare RDS files for reading from the the 10X matrix
printf "Converting mtx to rds..."
./prepRDS.R
printf "done\n"

## Run empty drops
printf "Running emptyDrops..."
../emptyDropsWrapper.R -i pbmc4k.rds -o pbmc4k_emptyDrops.csv
printf "done\n"

## Check the output md5 checksum
## Note that empty drops is based on MC simulation and therefore
## the output is not guaranteed to be deterministic, in our
## test it was and therefore the following checksum always
## validates
printf "Verifying checksum..."
md5out=`md5sum pbmc4k_emptyDrops.csv | cut -f 1 -d ' '`
if [ "$md5out" = "c2b0c8f24b9f8383f7de710774c9bbab" ];
then
    echo 'PASSED'
else
    echo 'FAIL'
fi
printf "done\n"

## Cleanup
rm -r pbmc4k_emptyDrops.csv pbmc4k_raw_gene_bc_matrices.tar.gz raw_gene_bc_matrices pbmc4k.rds
