## Author: Nikolas Barkas
## Date: December 2018
## Description: Analysis of HCA BM1_1 dataset with Seurat and
##   comparison between V2 and V3

## Load libraries
library(nbHelpers)
library(Seurat)

## Paths to load the filtered matrices from
v2.path <- '/home/nbarkas/storage3/nbarkas/work/cmpCellrangerVsOptimus_Gencode/hca_bm_s1/cr_v2/hca_bm_s1_cr2/outs/filtered_gene_bc_matrices/GRCh38'
v3.path <- '/home/nbarkas/storage3/nbarkas/work/cmpCellrangerVsOptimus_Gencode/hca_bm_s1/cr_v3/hca_bm_s1_cr3/outs/filtered_feature_bc_matrix'

## Load the matrices
v2.mat <- read10xMatrix(v2.path,'V2')
v3.mat <- read10xMatrix(v3.path,'V3')

## Check size of matrices
dim(v2.mat)
dim(v3.mat)

## find the cells in V3 that are not in V2
extracells <- setdiff(colnames(v3.mat), colnames(v2.mat))
length(extracells)

## Make the gene names unique to allow subsequent processing
rownames(v3.mat) <- make.unique(rownames(v3.mat))


## Process with Seurat -- Standard Analysis
hca_bm1_v3 <- CreateSeuratObject(raw.data = v3.mat, min.cells=0, min.genes=0, project='HCA_BM_V3')
hca_bm1_v3 <- NormalizeData(object = hca_bm1_v3, normalization.method='LogNormalize', scale.factor=10000)
hca_bm1_v3 <- FilterCells(object = hca_bm1_v3, subset.names=c("nGene"),low.thresholds = c(200), high.thresholds = Inf)
hca_bm1_v3 <- FindVariableGenes(object = hca_bm1_v3, mean.function= ExpMean, dispersion.function = LogVMR, x.low.cutoff=0.125, x.high.cutoff = 3, y.cutoff = 0.5)
hca_bm1_v3 <- ScaleData(object = hca_bm1_v3, vars.to.regress = c('nUMI'))
hca_bm1_v3 <- RunPCA(object = hca_bm1_v3, pc.genes = hca_bm1_v3@var.genes, do.print=TRUE, pcs.print=1:5,genes.print=5)
hca_bm1_v3 <- FindClusters(object=hca_bm1_v3, reduction.type='pca', dims.use=1:10, resolution=0.6, print.output=0, save.SNN=TRUE)
hca_bm1_v3 <- RunTSNE(object = hca_bm1_v3, dims.use=1:10, do.fast=TRUE)

## Extract the tSNE coordinates from the Seurat object
tsne.emb <- as.data.frame(hca_bm1_v3@dr$tsne@cell.embeddings)

## Extract the clusters from the seurat cluster
cl <- hca_bm1_v3@meta.data$res
names(cl) <- rownames(hca_bm1_v3@meta.data)

## Check that the order is identical
all(rownames(tsne.emb) == names(cl))

## append the clustering
tsne.emb$cl <- cl
tsne.emb$extra <- rownames(tsne.emb) %in% extracells

## Plot the tSNE with the clusters
p1 <- ggplot(tsne.emb, aes(x=tSNE_1, y=tSNE_2, color=as.factor(cl))) + geom_point()
ggsave('cluster.png',w=7,h=7,plot=p1)
p1

## Plot the tSNE annotating new cells
p2 <- ggplot(tsne.emb, aes(x=tSNE_1, y=tSNE_2, color=as.factor(extra))) + geom_point()
ggsave('extra.png',w=7,h=7,plot=p2)
p2

## Save final state
fastSave::preserve.state()

