# README for cellranger_compare_v2_v3

The following outlines the directory structure and describes the use of the script files.

docker/
	This directory contains the scripts required to build the docker image on which
	the R analysis is run on. The docker image is based on ubuntu 16.04 and R version 3.5.1.
	The docker file will install dependencies for Seurat and pagoda2 and install specific
	commits (current as of the time of generation of the image) of these two packages. The R
	requirements for pagoda and seurat are in the pagoda2.deps.r and seurat.deps.r files respectively.

	To build this image run:

	docker build -t quay.io/humancellatlas/secondary-analysis-cellranger-comparison-v2-v3

	The build image can also be found on quay with the above tag

README.md

        This file explaining the structure of the directories and steps required to reproduce the analysis

runs

	This folder contains all the scripts required to reproduce the analysis. This includes both the cellranger
	and the R (Custom, Seurat and Pagoda2 parts). The folder is organised by cellranger run name (dataset and
	reference used). In order to run the cellranger counting the v2 and v3 versions of cell ranger need to be
	downloaded independly. The downloads are available from the 10X website (https://www.10xgenomics.com).
	The versions of cellranger used for this analysis were: v2: v2.2.0 nd v3: v3.0.0.

	IMPORTANT NOTE: The run environment for these scripts is provided by the docker image above:
	quay.io/humancellatlas/secondary-analysis-cellranger-comparison-v2-v3

    cmpCellRanger2vs3_GRCh38.rel-93

	This folder contains the scripts for running cellranger (both versions 2 and 3) reference generation
	counting and them examining the resulting data count matrices in R.

	The 10X 4k PBMC samples are used for this analysis. These can be downloaded from the 10X website.

	In 01-mkref/, reffiles/getfile.sh
	downloads a specific version of the reference fasta and gtf annotation. The scripts v2/buildrefV2.sh and
	v3/buildrefV3.sh generate cellranger annotations using the respective cellranger version.

	In 02-count there are 3 run directories: execV2_annotV2, execV3_annotV2 and execV3_annotV3. These
	contain the counting runs for cellranger for different versions of cellranger and annotation as denoted
	in the name. For example execV2_annotV2 uses v2 cellranger generated annotation and the v2 of cellranger
	to do the counting.

	In 03-examine we examine the resulting matrices from the counting and generate the figures for the
	report. The main analysis is in checkFilteredMatrices.R. checkRawMatrices.R contains  a smaller
	version of the analysis performed on the raw matrices. functions.R contains auxilary functions
	used in the analysis. seurat_analsis.R contains the script that performs analysis of the dataset using
	Seurat.

    cmpCellRanger2vs3_HCA_BM1_Gencode

	This directory has a structure very similar to that of the directory cmpCellRanger2vs3_GRCh38.rel-93/
	the main difference is that it does not contain the reference generation step. Instead for this run
	the gencode v27 annotation is used. This is the annotation version that is used by the optimus pipeline
	and is used here for direct comparison. The reference used for this analysis can be downloaded from:

	gs:///hca-dcp-mint-test-data/reference/GRCh38_Gencode/GRCh38_GencodeV27_Primary_CellRanger.tar.gz

	Also, for this analysis the HCA_BM1 sample is used for this analysis. Specifically the following files
	were used:

	MantonBM1_HiSeq_1_S1_L007_I1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_I1_001.fastq.gz
	MantonBM1_HiSeq_1_S1_L007_R1_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R1_001.fastq.gz
	MantonBM1_HiSeq_1_S1_L007_R2_001.fastq.gz  MantonBM1_HiSeq_1_S1_L008_R2_001.fastq.gz

	These files can be dowloaded from the HCA preview website as part of the 'Census of Immune Cells'
	dataset. ( https://preview.data.humancellatlas.org/ ). 

    cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93

	This directory contains an analysis similar to that found in cmpCellRanger2vs3_HCA_BM1_Gencode/
	but with the same annotation used in cmpCellRanger2vs3_GRCh38.rel-93