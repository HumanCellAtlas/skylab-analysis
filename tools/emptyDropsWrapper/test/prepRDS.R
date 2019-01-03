#!/usr/bin/env Rscript
library(nbHelpers)
d <- read10xMatrix('raw_gene_bc_matrices/GRCh38', version='V2')
saveRDS(d, 'pbmc4k.rds')
