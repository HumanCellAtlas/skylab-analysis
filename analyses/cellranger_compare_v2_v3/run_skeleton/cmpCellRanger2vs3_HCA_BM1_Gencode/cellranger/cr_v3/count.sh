#!/bin/bash

reference=/home/nbarkas/storage3/nbarkas/references/GRCh38_GencodeV27_Primary_CellRanger/GRCh38/
cellrangerExec=/home/nbarkas/software/cellranger-3.0.0/cellranger
fastqs=/home/nbarkas/storage3/nbarkas/work/cmpCellrangerVsOptimus_Gencode/hca_bm_s1/input/

$cellrangerExec count --id hca_bm_s1_cr3 \
		--transcriptome $reference \
		--fastqs $fastqs \
		--expect-cells 5500 \
		--localcores 8 --localmem 24 \
		--nosecondary
	    
