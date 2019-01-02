#!/bin/bash

reference=/home/nbarkas/storage2/nbarkas/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/v2/GRCh38
cellrangerExec=/home/nbarkas/software/cellranger-2.2.0/cellranger
fastqs=/home/nbarkas/storage2/nbarkas/cmpCellRanger2vs3_GRCh38.rel-93/data/pbmc4k/fastqs

$cellrangerExec count --id pbmc_4k \
		--transcriptome $reference \
		--fastqs $fastqs \
		--expect-cells 5000 \
		--localcores 8 --localmem 24 \
		--nosecondary
	    
