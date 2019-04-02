# Optimus Acceptance Report

This folder contains the scripts used for the generate of figures for the Optimus acceptace report. It is structured as follows:

##check-cellranger-input-annotatinon
Script to check the number of input genes in the filtered cellranger annotation. It loads the annotation gtf from cellranger (after the filtering) and counts unique genes. This is used to check that the number of genes that the cellranger output has is limited due to the gtf input.

##manual-cellranger-runs
A list of commands issued to perform processing of the test datasets with cellranger version 3. These are manually run using cellranger.

## matrix-comparison
Scripts that compare matrices and generate the main figures for the report. These scripts compare the output between cellranger and optimus.

## optimus-reproducibility
Scripts that compare the reproducibility of optimus accross identical runs. They are based on three identical runs of the t_4k dataset from cellranger.

## umi-check
Script that checks the effect of applying or not umi correction to the Optimus output. The standard umi-corrected output of optimus is checked against an alternative counting of the dataset that is performed using the uncorrected umi tags.