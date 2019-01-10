# Cellranger V2 vs V3 analysis walkthrough

This document provides instructions for reproducing the comparison between version 2 and 3 of cellranger.

Software requirements to run the following are a working installation of docker and git. You will also need a machine capable of running cellranger and
at least 400GB of disk space.

The following instructions have been generated for linux systems but should also work on windows systems with minor modifications.

```
############################
## Setup the docker image
############################

## First of all define a location where the analysis will be stored in the host
## machine by setting the following environment variable
mkdir v2v3analysisData
export data=${PWD}/v2v3analysisData

## Obtain a copy of the git repository and export the required branch
## TODO: Update the branch once its merged
git clone https://github.com/HumanCellAtlas/skylab-analysis.git
cd skylab-analysis
export repos_root=${PWD}
git checkout nb-add-v2-v3-report-code
cd analyses/cellranger_compare_v2_v3/docker/

## Build the docker image. This step take a bit less than an hour to complete
## Alternative the built image is availabel on quay.io and can be pulled
## using the tag 'quay.io/humancellatlas/secondary-analysis-cellranger-comparison-v2-v3'
docker build -t quay.io/humancellatlas/secondary-analysis-cellranger-comparison-v2-v3 .

## Run the docker container mounting the checked out repository and the
## data folder specified above
docker run \
	--mount type=bind,source=${data},destination=/data/ \
	--mount type=bind,source=${repos_root},destination=/repos/ \
	-it quay.io/humancellatlas/secondary-analysis-cellranger-comparison-v2-v3 \
	/bin/bash

## Make a copy of the skeleton analysis in data 
cp -r /repos/analyses/cellranger_compare_v2_v3/run_skeleton /data/
```

At this point you need to download cellranger version 2.2.0 and 3.0.0 manually from teh 10X website.
Navigate to https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/ and obtain
links to the tar files required.

Then continue with the instructions below

```
cd /data/
mkdir software

## Download the files the links to which you obtained in this directory

## Extract the files
tar xzf cellranger-2.2.0.tar.gz
tar xzf cellranger-3.0.1.tar.gz

## Set the global cellranger executable paths
export cellranger_v2_exec=/data/software/cellranger-2.2.0/cellranger
export cellranger_v3_exec=/data/software/cellranger-3.0.1/cellranger


#########################
## Build the reference
#########################

## Download the fasta and gtf files required for building the references
cd  /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/reffiles
./getfiles

## Set the enviroment variable specifying the location of the reference files
export refdir_GRCh38_rel93=/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/reffiles

## Build V2 reference
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/v2
./buildrefV2.sh

## Build V3 reference
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/v3
./buildrefV3.sh

## Export paths to references
export referenceV2=/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/v2/GRCh38
export referenceV3=/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/v3/GRCh38

########################
## Download data
########################

## 10X 4k PBMC
cd /data/run_skeleton/fastqs/pbmc4k
./getData.sh
tar xf pbmc4k_fastqs.tar
export data_pbmc4k=/data/run_skeleton/fastqs/pbmc4k/fastqs

## HCA BM1
```
This step need to be performed manually.
Navigate to https://preview.data.humancellatlas.org/ adn obtain the
samples names MantonBM1_1_* (this should include two lanes of data
(lane 7 and lane8) into the current directory.
```
## Export the location of the hcabm1 data
export data_hcabm1=/data/run_skeleton/fastqs/hcaBM1

################################################
## Analysis PBMC4k with GRCh38.rel-93 reference
################################################

## Run counting with cellranger
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV2_annotV2
./count_execV2_annotV2.sh
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV2
./count_execV3_annotV2.sh
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV3
./count_execV3_annotV3.sh

## Do analysis in R
cd /data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/03-examine
Rscript checkFilteredMatrices.R
Rscript checkRawMatrices.R
Rscript seurat_analysis.R

######################################################
## Analysis of HCA_BM1 with GRCh38.rel-93
######################################################

## Run counting with cellranger
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93/execV2_annotV2
./count_execV2_annotV2.sh

cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93/execV3_annotV3
./count_execV3_annotV3.sh

## Run analysis
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93/compare
./analysis_seurat.R

###########################################################
## Analysis HCA_BM1 with Gencode used by the SS2 pipeline
###########################################################

# Get the GencodeV27 reference, this is a pre-build cellranger reference
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_Gencode/reference
gsutil cp gs://hca-dcp-mint-test-data/reference/GRCh38_Gencode/GRCh38_GencodeV27_Primary_CellRanger.tar .
tar xf GRCh38_GencodeV27_Primary_CellRanger.tar

## Export the location of the reference
export referenceGencode=/data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_Gencode/reference/GRCh38

## Run V2
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_Gencode/cellranger/cr_v2
./count.sh

## Run V3
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_Gencode/cellranger/cr_v3
./count.sh

## Run the comparative analysis
cd /data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_Gencode/analysis
Rscript hca_bm1_cr_v2_vs_v3_seraut.R


```