# Optimus Acceptance Report

This folder contains the scripts used for the generate of figures for the Optimus acceptace report. It is structured as 
\follows:

## Analysis Walkthrough

## Datasets used
The following datasets are used in this analysis:

* 4k_panT
* 4k_pbmc
* 8k_pbmc
* MantonBM1_1 
* ischaemia_1

### 4k_panT
This is the 10X "4k Pan T Cells from a Healthy Donor" Dataset. It can be downloaded from 
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_4k).

### 4k_pbmc
This is the 10X "4k PBMCs from a Healthy Donor" dataset. It can be downloaded from
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k).

### 8k_pbmc 
This is the 10X "8k PBMCs from a Healthy Donor" dataset. It can be downloaded from
[here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k).

### MantonBM1_1
This is one sample from the "Immune Cell Atlas Dataset" from the Human Cell Atlas Project.
As of this writing it can downloaded from 
[https://prod.data.humancellatlas.org/](https://prod.data.humancellatlas.org/).
Specifically the inputs files comprise: 
```
MantonBM1_HiSeq_1_S1_L007_I1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_I1_001.fastq.gz
MantonBM1_HiSeq_1_S1_L007_R1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R1_001.fastq.gz
MantonBM1_HiSeq_1_S1_L007_R2_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R2_001.fastq.gz
```

### Ischaemia_1
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
In order to reproduce the analysis the above datasets needs to run through both Optimus and Cell Ranger. For the 
reproducibility testing the 4k_panT dataset needs to be run through Optimus 3 times with caching turned off.

The Cell Ranger version used in this analysis is 3.0.2. The Optimus version is v1.0.0.
The reference used is GRCh38 with the Gencode V27 annotation, filtered for the requirements of Cell Ranger. 

The reference used can be downloaded using: 
```
gsutil cp gs://hca-dcp-mint-test-data/reference/GRCh38_Gencode/GRCh38_GencodeV27_Primary.tar .
```

The commands for performing the Cell Ranger runs can be found here in the ```manual-cellranger-runs``` directory.

## Replicating the downstream data processing

### Converting the data for input into R
The count output of Optimus can be converted into an 'rds' file that contains a compressed sparse matrix format
object that can be read into R. This utility can be found in the skylab repository 
[here](https://github.com/HumanCellAtlas/skylab/blob/master/docker/emptydrops/npz2rds/npz2txt.py).

### Main matrix comparison and data inspection
The majority of the code for inspecting the datasets can be found under the ```matrix-comparison/``` directory. These 
are R scripts that accept as input the location of the files generated from the primary analysis and generate plots
and pagoda2 apps.

### Comparison of Optimus Matrices with and without UMI correction
Script that checks the effect of applying or not umi correction to the Optimus output can be found under the 
```umi-check``` directory. The standard umi-corrected output of optimus is checked against an alternative counting of 
the dataset that is performed using the uncorrected umi 
tags.

### Check Cell Ranger input annotation
The Script to check the number of input genes in the filtered Cell Ranger annotation can be foung under
```check-cellranger-input-annotation```. It loads the annotation gtf from 
Cell Ranger (after the filtering) and counts unique genes. This is used to check that the number of genes that the 
Cell Ranger output has is limited due to the GTF input (and that the numbers match).

### Checking reproducibility
The scripts that check the reproducibility of the output can be found under ``optimus-reproducibility``.
they are based on three identical runs of the t_4k dataset from cellranger.
