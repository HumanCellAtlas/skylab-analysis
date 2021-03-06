## Description: Analysis of V2 and V3 output of the dataset with Seurat
## Author: Nikolas Barkas
## Date: December 2018

library(Seurat)
library(nbHelpers)

## Input file paths
v2.path <- '/data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93/execV2_annotV2/hca_bm_s1/outs/filtered_gene_bc_matrices/GRCh38'
v3.path <- '/data/run_skeleton/cmpCellRanger2vs3_HCA_BM1_GRCh38.rel-93/execV3_annotV3/hca_bm_s1/outs/filtered_feature_bc_matrix'


## Read in the matrices
v2.mat <- read10xMatrix(v2.path,'V2')
v3.mat <- read10xMatrix(v3.path,'V3')

## Check the matrix dimenstions
dim(v2.mat)
dim(v3.mat)

## Get the extra numbers of cells
extracells <- setdiff(colnames(v3.mat), colnames(v2.mat))
length(extracells)

## Make gene names unique
rownames(v3.mat) <- make.unique(rownames(v3.mat))

## Process with Seurat
hca_bm1_v3 <- CreateSeuratObject(raw.data = v3.mat, min.cells=0, min.genes=0, project='HCA_BM_V3')
hca_bm1_v3 <- NormalizeData(object = hca_bm1_v3, normalization.method='LogNormalize', scale.factor=10000)
hca_bm1_v3 <- FilterCells(object = hca_bm1_v3, subset.names=c("nGene"),
                          low.thresholds = c(200), high.thresholds = Inf)
hca_bm1_v3 <- FindVariableGenes(object = hca_bm1_v3, mean.function= ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff=0.125, x.high.cutoff = 3, y.cutoff = 0.5)
hca_bm1_v3 <- ScaleData(object = hca_bm1_v3, vars.to.regress = c('nUMI'))
hca_bm1_v3 <- RunPCA(object = hca_bm1_v3, pc.genes = hca_bm1_v3@var.genes, do.print=TRUE, pcs.print=1:5,
                     genes.print=5)
hca_bm1_v3 <- FindClusters(object=hca_bm1_v3, reduction.type='pca', dims.use=1:10,
                           resolution=0.6, print.output=0, save.SNN=TRUE)
hca_bm1_v3 <- RunTSNE(object = hca_bm1_v3, dims.use=1:10, do.fast=TRUE)

## Extract tSNE
tsne.emb <- as.data.frame(hca_bm1_v3@dr$tsne@cell.embeddings)

## Get the clusters
cl <- hca_bm1_v3@meta.data$res
names(cl) <- rownames(hca_bm1_v3@meta.data)

## Check that they are in the same order
all(rownames(tsne.emb) == names(cl))

## append the clustering
tsne.emb$cl <- cl
tsne.emb$extra <- rownames(tsne.emb) %in% extracells

## Plot the tSNE wby cluser
p1 <- ggplot(tsne.emb, aes(x=tSNE_1, y=tSNE_2, color=as.factor(cl))) + geom_point()
ggsave('cluster.png',w=7,h=7,plot=p1)
p1

## Plot the tsen annotating extra cells
p2 <- ggplot(tsne.emb, aes(x=tSNE_1, y=tSNE_2, color=as.factor(extra))) + geom_point()
ggsave('extra.png',w=7,h=7,plot=p2)
p2

## Plot the tSNE with clusters and labels
pl1 <- TSNEPlot(hca_bm1_v3,do.label=TRUE)
ggsave('tsne.cl.labeled.png',pl1)


## Find markers for all clusters
lapply(namedLevels(hca_bm1_v3@ident), function(id) {
    mrks <- FindMarkers(hca_bm1_v3,ident.1=id, ident.2=NULL, only.pos=TRUE)
    write.csv(mrks, paste0('cl.markers/',id,'.csv'))
})

## Make a plot of the number of cells in each cluster that are only identifed in V3
cl <- hca_bm1_v3@meta.data$res
names(cl) <- rownames(hca_bm1_v3@meta.data)
tmp.df1 <- data.frame(
    cellid=names(cl),
    cl
)
tmp.df1$extra <- rownames(tmp.df1) %in% extracells
aggr <- aggregate(tmp.df1$extra, by=list(tmp.df1$cl), function(x) {sum(x)/length(x)})

## plot
ggplot(aggr, aes(x=Group.1,y=x)) + geom_bar(stat='identity') +
    scale_y_continuous(lim=c(0,1), name='Percent of Cluster only identified by V3')  +
    scale_x_discrete(name='Cluster identifier')
ggsave('new.in.cluster.png')
