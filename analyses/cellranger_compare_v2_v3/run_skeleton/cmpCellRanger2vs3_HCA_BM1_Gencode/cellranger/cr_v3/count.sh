#!/bin/bash

$cellranger_v3_exec count --id hca_bm_s1_cr3 \
		--transcriptome $referenceGencode \
		--fastqs $data_hcabm1 \
		--expect-cells 5500 \
		--localcores 8 --localmem 24 \
		--nosecondary
	    
