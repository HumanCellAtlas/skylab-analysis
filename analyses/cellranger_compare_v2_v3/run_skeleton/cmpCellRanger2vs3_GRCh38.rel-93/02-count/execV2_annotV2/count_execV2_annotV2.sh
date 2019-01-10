#!/bin/bash

${cellranger_v2_exec} count --id pbmc_4k \
		--transcriptome $referenceV2 \
		--fastqs $data_pbmc4k \
		--expect-cells 5000 \
		--localcores 8 --localmem 24 \
		--nosecondary
	    
