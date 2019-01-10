#!/bin/bash

$cellranger_v3_exec count --id pbmc_4k \
		--transcriptome $referenceV2 \
		--fastqs $data_pbmc4k \
		--expect-cells 5000 \
		--localcores 16 --localmem 60 \
		--nosecondary
	    
