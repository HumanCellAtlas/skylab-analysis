## Comparison analysis of the cellranger and optimus output matrices
## for the 8k PBMC dataset. This is meant to be run interactively

## Load required libraries
library(Matrix)
library(nbHelpers)
library(pheatmap)
library(pagoda2)
library(fastSave)

## Load the data

## This is the direct 10X output matrix
cr.path <- "/home/nbarkas/disk3/8k_pbmc/run/pbmc8k/outs/raw_feature_bc_matrix"
mat.10x <- read10xMatrix(cr.path,version='V3')

## This is the optimus matrix from the same run with converted to rds format                                                                    
## with npz2rds utility that can be found under the skylab repository under                                                                     
## /docker/emptydrops/npz2rds/   
umi.on.path <- '/home/nbarkas/disk2/umi-check/8k_pbmc/MergeCountFiles/call-MergeCountFiles/convert/out.rds'
cm.on <- readRDS(umi.on.path)

## Fix the 10X names that have an -1 appended at the end of the cell 
colnames(mat.10x) <- strpart(colnames(mat.10x), '-', 1)

## Compare the dimensions of the matrix 
dim(mat.10x)
dim(cm.on)

## Find the non-empty barcodes
non.emtpy.10x.bc <- (colnames(mat.10x)[colSums(mat.10x) > 0])

## Non-empty cells in 10X match the cells that are returned by Optimus
bc.10x.absent.cr <- setdiff(non.emtpy.10x.bc,rownames(cm.on))
length(bc.10x.absent.cr)
summary(colSums(mat.10x[,bc.10x.absent.cr]))

## Barcodes in optimus that are not present in the 10x output after removing all 0 cells
bc.cr.absent.cr <- setdiff(rownames(cm.on), non.emtpy.10x.bc)
length(bc.cr.absent.cr)
summary(rowSums(cm.on[bc.cr.absent.cr,]))

## Check the Matrix dimensions and the extra genes that Optimus returns
## This is supporting the main text section discussing differences between the
## matrix sizes
dim(mat.10x)
dim(cm.on)
head(rownames(mat.10x))
head(colnames(cm.on))
length(setdiff(rownames(mat.10x),colnames(cm.on)))
length(setdiff(colnames(cm.on), rownames(mat.10x)))
optimus.extra.genes <- setdiff(colnames(cm.on), rownames(mat.10x))
cm.on.extra.genes <- cm.on[,optimus.extra.genes]
cm.non.extra.genes <- cm.on[,!colnames(cm.on) %in% optimus.extra.genes]
sum(cm.on.extra.genes)/sum(cm.non.extra.genes)

## What percent of the total genes are extra?
length(optimus.extra.genes) / ncol(cm.on)

## Load the emptyDrops output from Optimus
ed.path <- "/home/nbarkas/disk2/emptydrops-check/data/optimus/8k_pbmc/RunEmptyDrops/empty_drops_result.csv"
ed.table <- read.csv(ed.path)
head(ed.table)

## Apply the updated cutoff
ed.table$IsCell <- ed.table$FDR < 0.01 & ed.table$Total > 100
keep.cells <- as.character(ed.table$CellId[ed.table$IsCell])
head(keep.cells)
length(keep.cells)

## Inspect the matrices
mat.10x[1:3,1:3]
cm.on[1:3,1:3]

## Subset both to keep cells from emptyDrops
mat.10x <- mat.10x[,keep.cells]
cm.on <- cm.on[keep.cells,]

## Transpose the optimus matrix
cm.on <- t(as(cm.on, 'dgCMatrix'))

## Compare the dimensions of the matrices
dim(mat.10x)
dim(cm.on)

## find and inspect the common genes
common.genes <- intersect(rownames(mat.10x), rownames(cm.on))
length(common.genes)
head(common.genes)

## Subset both matrices to the common genes
cm.on <- cm.on[common.genes,]
mat.10x <- mat.10x[common.genes,]

## Check the matrix dimensions and inspect the matrices
dim(mat.10x)
dim(cm.on)
mat.10x[1:10,1:10]
cm.on[1:10,1:10]

## Generate plots for the report
png('mat.img/genesum.cmp.png')
plot(log10(rowSums(mat.10x) + 1), log10(rowSums(cm.on) +1), xlab='10x', ylab='optimus',
     main='Sum of per gene counts (8k PBMC) in log scale')
abline(0,1,col='red')
dev.off()

## Calculate correlation coefficient for the above plot
df1 <- data.frame(cr=log10(rowSums(mat.10x) + 1), opt=log10(rowSums(cm.on) +1))
df1 <- df1[df1$cr > 1e-6& df1$opt > 1e-6,]
cor(df1$cr ,df1$opt)

png('mat.img/genesum.linear.cmp.png')
plot((rowSums(mat.10x) + 1), (rowSums(cm.on) +1), xlab='10x', ylab='optimus',
     main='Sum of per gene counts (8k PBMC)')
abline(0,1,col='red')
dev.off()

png('mat.img/perGeneAbsDeltaHist.png')
hist(log10(rowSums(abs(mat.10x - cm.on))+1),breaks=100,
     main='Histogram of per gene absolute matrix differences',
     xlab='log10(perGene sum of differences)')
dev.off()

## Get raw count cross correlation coeffient for the two matrices
cells.xcor <- sample(colnames(mat.10x),1000)
xcor <- cor(as(mat.10x[,cells.xcor],'matrix'),as(cm.on[,cells.xcor],'matrix'),method='spearman')

## inspect
xcor[1:10,1:10]

## Calculate summary stats of cross correlation for non identical cells
xcor.2 <- xcor
diag(xcor.2) <- c(0)
summary(xcor.2)

## Plot the distribution of same cell cross correlation coefficients between the two matrices
png('mat.img/hist.same.cell.xcor.png')
hist(diag(xcor),breaks=50,main='',xlab='Same cell cross correlation')
dev.off()

## Confirm that each cell cross-correlates maximally with itself
sum(apply(xcor, 1, which.max) == seq_along(rownames(xcor)))

## Plot the unordered heatmap of the cell cross-correlation
pheatmap(xcor,cluster_rows=FALSE, cluster_cols=FALSE, file='mat.img/xcor.heatmap.png', show_rownames=FALSE, show_colnames=FALSE)

## Do the same with columns and row clustering to confirm that we see the expected population structure
pheatmap(xcor,cluster_rows=T, cluster_cols=T, file='mat.img/xcor.heatmap.cluster.png', show_rownames=FALSE, show_colnames=FALSE)

## Normalize the matrices and recreate the heatmap
matA <- as(mat.10x[,cells.xcor],'matrix')
matB <- as(cm.on[,cells.xcor],'matrix')
matA <- sweep(matA, 2, apply(matA, 2, sum), '/')
matB <- sweep(matB, 2, apply(matA, 2, sum), '/')
xcor.v2 <- cor(matA,matB,method='spearman')
pheatmap(xcor.v2,cluster_rows=T, cluster_cols=T, file='mat.img/xcor.heatmap.cluster.v2.png', show_rownames=FALSE, show_colnames=FALSE)

## Generate pagoda2 apps
p2.cr <- basicP2proc(mat.10x[,colSums(mat.10x) > 10],n.cores=16)
p2.cr.web <- basicP2web(p2.cr,app.title='10x',n.cores=16)
p2.opt <- basicP2proc(cm.on[,colSums(cm.on) > 10],n.cores=16)
p2.opt.web <- basicP2web(p2.opt, app.title='opt',n.cores=16)

## Inspect the pagoda2 apps manually -- note this will start a rook server
show.app(p2.cr.web,browse=FALSE,name='cr')
show.app(p2.opt.web,browse=FALSE,name='opt')

## Get the delta matrix masking out entries that are 0 in both matrices
mask <- mat.10x + cm.on == 0
delta.mat <- as.matrix(mat.10x - cm.on)
mask <- as.matrix(mask)
delta.mat[mask] <- c(NA)
rs <- log10(rowSums(abs(delta.mat),na.rm=TRUE)+1)
summary(rs)

## This is the per cell mean change
10^0.7688/8201
dim(delta.mat)

## Plot Per Gene Sum of absolute per element differences between Optimus and Cellranger
png('mat.img/delta.hist.3.png')
hist(rs,breaks=50,xlab="Per Gene Sum of absolute per element differences between Optimus and Cellranger")
dev.off()

## Optionally preserver the final state
preserve.state(n.cores=16)
