# Optimus Acceptance Report

This folder contains the scripts used for the generate of figures for the Optimus acceptace report. It is structured as 
\follows:

## Analysis Walkthrough

## Datasets used
The following datasets are used in this analysis:

* 4k_panT
* 4k_pbmc ()
* 8k_pbmc ()
* MantonBM1_1 
* ischaemia_1

### 4k_panT
This is the 10X "4k Pan T Cells from a Healthy Donor" Dataset. It can be downloaded from 
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_4k).

### 4k_pbmc
This is the 10X "4k PBMCs from a Healthy Donor" dataset. It can be downloaded from
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k).

## 8k_pbmc 
This is the 10X "8k PBMCs from a Healthy Donor" dataset. It can be downloaded from
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k).

## MantonBM1_1
This is one sample from the "Immune Cell Atlas Dataset" from the Human Cell Atlas Project.
As of this writing it can downloaded from 
[https://prod.data.humancellatlas.org/](https://prod.data.humancellatlas.org/).
Specifically the inputs files comprise: 

```
MantonBM1_HiSeq_1_S1_L007_I1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_I1_001.fastq.gz
MantonBM1_HiSeq_1_S1_L007_R1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R1_001.fastq.gz
MantonBM1_HiSeq_1_S1_L007_R2_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R2_001.fastq.gz
```

## Ischaemia_1
This is one sample from the "Ischaemic sensitivity of human tissue by single cell RNA seq".
As of this writing it can downloaded from 
[https://prod.data.humancellatlas.org/](https://prod.data.humancellatlas.org/).
It comprises of the following files:

```
HCATisStabAug177376561_S1_L001_I1_001.fastq.gz	
HCATisStabAug177376561_S1_L001_R2_001.fastq.gz
HCATisStabAug177376561_S1_L001_R1_001.fastq.gz
```

## Replicating the primary data processing
In order to reproduce the analysis the five following datasets must be run through both Cell Ranger and Optimus:

The Cell Ranger version used in this analysis is 3.0.2. The Optimus version is v1.0.0.
The reference used is GRCh38 with the Gencode V27 annotation, filtered for the requirements of Cell Ranger. 

The reference used can be downloaded using: 
```
gsutil cp gs://hca-dcp-mint-test-data/reference/GRCh38_Gencode/GRCh38_GencodeV27_Primary.tar .
```
###check-cellranger-input-annotatinon
Script to check the number of input genes in the filtered cellranger annotation. It loads the annotation gtf from 
cellranger (after the filtering) and counts unique genes. This is used to check that the number of genes that the 
cellranger output has is limited due to the gtf input.

###manual-cellranger-runs
A list of commands issued to perform processing of the test datasets with cellranger version 3. These are manually run 
using cellranger.

###matrix-comparison
Scripts that compare matrices and generate the main figures for the report. These scripts compare the output between 
cellranger and optimus.

###optimus-reproducibility
Scripts that compare the reproducibility of optimus accross identical runs. They are based on three identical runs of
the t_4k dataset from cellranger.

###umi-check
Script that checks the effect of applying or not umi correction to the Optimus output. The standard umi-corrected 
output of optimus is checked against an alternative counting of the dataset that is performed using the uncorrected umi 
tags.

