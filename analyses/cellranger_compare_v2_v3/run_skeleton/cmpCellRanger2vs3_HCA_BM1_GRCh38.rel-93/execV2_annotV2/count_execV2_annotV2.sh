#!/bin/bash

$cellranger_v2_exec count --id hca_bm_s1 \
		--transcriptome $referenceV2 \
		--fastqs $data_hcabm1 \
		--expect-cells 5000 \
		--localcores 8 --localmem 24 \
		--nosecondary
	    
