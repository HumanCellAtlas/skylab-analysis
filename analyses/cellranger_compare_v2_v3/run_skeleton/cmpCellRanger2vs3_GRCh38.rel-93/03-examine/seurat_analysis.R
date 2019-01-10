

## Load the libraries
library('Seurat')
library('nbHelpers')
library(ggplot2)
library('fastSave')

## Load helper funcgtions
ggsource('functions.R')

## The V2 and V3 filtered input files
v2filt.file <- '/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV2_annotV2/pbmc_4k/outs/filtered_feature_bc_matrix'
v3filt.file <- '/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV3/pbmc_4k/outs/filtered_feature_bc_matrix'


#install.packages('Seurat') # remember to install hdf5 when dockerizing

## Read in the input data
v2.mat <- read10xMatrix(v2filt.file, 'V2')
v3.mat <- read10xMatrix(v3filt.file, 'V3')

## Make a note of the extra cells that pass the filtering
extracells <- setdiff(colnames(v3.mat), colnames(v2.mat))

## Get the rownames
rownames(v3.mat) <- make.unique(rownames(v3.mat))

## Set the display and open an X window
setDisplay('12.0')
X11()

## Process with Seurat
pbmc <- CreateSeuratObject(raw.data = v3.mat, min.cells=0,min.genes=0, project='10X_PBMC_4k')

## Normalize the data
pbmc <- NormalizeData(object = pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)

## Filter the cells 
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene"),
                    low.thresholds = c(200), high.thresholds = c(2500))

## Find variable genes
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

## Only regress out UMIs
pbmc <- ScaleData(object = pbmc, vars.to.regress = c('nUMI'))

## Run PCA
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)

## Find the clusters
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10,
                         resolution = 0.6, print.output = 0, save.SNN = TRUE)

## Run tSNE
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

## Get the embedding form the Seurat obj
tsne.emb <- as.data.frame(pbmc@dr$tsne@cell.embeddings)

## Get the clusters called by Seurat above
cl <-  pbmc@meta.data$res
names(cl) <- rownames(pbmc@meta.data)

head(tsne.emb)
## Verify all in same order
all(rownames(tsne.emb) == names(cl))

## Append the clustering
tsne.emb$cl <- cl


## Annotate the extra cells
tsne.emb$extra <- rownames(tsne.emb) %in% extracells

## Plot the tSNEs
p1 <- ggplot(tsne.emb, aes(x=tSNE_1, y= tSNE_2, color=as.factor(cl))) + geom_point()
ggsave('outs/tsne.clusters.png',p1)

p2 <- ggplot(tsne.emb, aes(x=tSNE_1, y= tSNE_2, color=as.factor(extra))) + geom_point()
ggsave('outs/tsne.extracells.png',p2)

preserve.state()
##  "savepoint_2018-11-29_22:00:59_29599.RDataFS"

load.lbzip2("savepoint_2018-11-29_22:00:59_29599.RDataFS")

## Perform differential Expression of the New Clusters

X11()

library(Seurat)

pl1 <- TSNEPlot(pbmc,do.label=TRUE)
class(pl1)

ggsave('tsne_annot.png',pl1)

getwd()

cl9.markers <- FindMarkers(pbmc, ident.1='9', ident.2=NULL, only.pos=TRUE)
cl9.markers <- subset(cl9.markers, p_val_adj < 0.05)
head(cl9.markers)
write.csv(cl9.markers, 'cl9.markers.csv')
dim(cl9.markers)

## Perform Go analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GOstats", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")

library('GOstats')
library('org.Hs.eg.db')
library('reshape2')

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
entrez2symbol <- melt(xx)
rm(x, xx)
names(entrez2symbol)

rownames(cl9.markers) -> selectedGenes
rownames(v3.mat) -> universeGenes
head(universeGenes)

selectedIDs <- entrez2symbol[match(selectedGenes, entrez2symbol$value),]$L1
universeIDs <- entrez2symbol[match(universeGenes, entrez2symbol$value),]$L1

head(selectedIDs)
head(universeIDs)

selectedIDs <- selectedIDs[!is.na(selectedIDs)]
universeIDs <- universeIDs[!is.na(universeIDs)]

type <- c('BP','MF')
names(type) <- type
type <- as.list(type)

length(selectedIDs)
length(universeIDs)

goResults <- lapply(type, function(t) {
    param <- new("GOHyperGParams",
        geneIds = selectedIDs,
        universeGeneIds = universeIDs,
        annotation= 'org.Hs.eg',
        ontology=t,
        pvalueCutoff=0.05,
        conditional = TRUE,
        testDirection='over')
    hyperGTest(param)
})


write.csv(as.data.frame(summary(goResults$BP)),'cl9.gobp.csv')
write.csv(as.data.frame(summary(goResults$MF)),'cl9.gomf.csv')



preserve.state()
## "savepoint_2018-12-19_16:08:05_3547.RDataFS"
