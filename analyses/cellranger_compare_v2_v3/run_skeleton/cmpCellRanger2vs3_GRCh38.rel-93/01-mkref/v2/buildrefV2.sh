#!/bin/bash

${cellranger_v2_exec} mkgtf ${refdir_GRCh38_rel93}/Homo_sapiens.GRCh38.93.gtf Homo_sapiens.GRCh38.93.filtered.gtf \
	   --attribute=gene_biotype:protein_coding \
	   --attribute=gene_biotype:lincRNA \
	   --attribute=gene_biotype:antisense \
	   --attribute=gene_biotype:IG_LV_gene \
	   --attribute=gene_biotype:IG_V_gene \
	   --attribute=gene_biotype:IG_V_pseudogene \
	   --attribute=gene_biotype:IG_D_gene \
	   --attribute=gene_biotype:IG_J_gene \
	   --attribute=gene_biotype:IG_J_pseudogene \
	   --attribute=gene_biotype:IG_C_gene \
	   --attribute=gene_biotype:IG_C_pseudogene \
	   --attribute=gene_biotype:TR_V_gene \
	   --attribute=gene_biotype:TR_V_pseudogene \
	   --attribute=gene_biotype:TR_D_gene \
	   --attribute=gene_biotype:TR_J_gene \
	   --attribute=gene_biotype:TR_J_pseudogene \
	   --attribute=gene_biotype:TR_C_gene

${cellranger_v2_exec} mkref --genome=GRCh38 \
	   --fasta=${refdir_GRCh38_rel93}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	   --genes=Homo_sapiens.GRCh38.93.filtered.gtf \
	   --nthreads 8 --memgb 24 \
	   --ref-version=1.2.0
