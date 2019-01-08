#!/bin/bash

reffilesdir=/home/nbarkas/storage2/nbarkas/cmpCellRanger2vs3_GRCh38.rel-93/01-mkref/reffiles
cellrangerExec=/home/nbarkas/software/cellranger-3.0.0/cellranger

${cellrangerExec} mkgtf ${reffilesdir}/Homo_sapiens.GRCh38.93.gtf Homo_sapiens.GRCh38.93.filtered.gtf \
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


${cellrangerExec} mkref --genome=GRCh38 \
	   --fasta=${reffilesdir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	   --genes=Homo_sapiens.GRCh38.93.filtered.gtf \
	   --nthreads 2 --memgb 24 \
	   --ref-version=3.0.0
