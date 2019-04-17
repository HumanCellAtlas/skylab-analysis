## Comparison analysis of the cellranger and optimus output matrices
## for the 4k PBMC dataset. This is meant to be run interactively.

## Load required libraries
library(Matrix)
library(nbHelpers)
library(pheatmap)
library(pagoda2)
library(fastSave)

## Load the data

## Paths

## Cellranger output matrix
cr.path <- "/home/nbarkas/disk3/4k_pbmc/run/pbmc4k/outs/raw_feature_bc_matrix/"

## This is the optimus matrix from the same run with converted to rds format
## with npz2rds utility that can be found under the skylab repository under
## /docker/emptydrops/npz2rds/
umi.on.path <- '/home/nbarkas/disk2/umi-check/4k_pbmc/data/convert/out.rds'

## This is the direct 10X output matrix
mat.10x <- read10xMatrix(cr.path,version='V3')
cm.on <- readRDS(umi.on.path)

## Fix the 10X names that have an -1 appended at the end of the cell
colnames(mat.10x) <- strpart(colnames(mat.10x), '-', 1)

## Compare the dimensions of the matrix
dim(mat.10x)
dim(cm.on)

## Load the emptyDrops output from optimus
ed.path <- "/home/nbarkas/disk2/check-empty-drops/4k_pbmc/emptyDrops/empty_drops_result.csv"
ed.table <- read.csv(ed.path)
head(ed.table)

## Impose updated cutoff and get a vector of the cells that we want to keep
ed.table$IsCell <- ed.table$FDR < 0.01 & ed.table$Total > 100
keep.cells <- as.character(ed.table$CellId[ed.table$IsCell])

## Inspect the kept cells
head(keep.cells)
length(keep.cells)

## Check the arrays and their orientations (they are not the same)
mat.10x[1:3,1:3]
cm.on[1:3,1:3]

## Subset both matrices to keep cells from emptyDrops
mat.10x <- mat.10x[,keep.cells]
cm.on <- cm.on[keep.cells,]

## Transpose the matrix from optimus
cm.on <- t(as(cm.on, 'dgCMatrix'))

## Confirm the dimentions and inspect the matrices
dim(mat.10x)
dim(cm.on)
cm.on[1:3,1:3]
mat.10x[1:3,1:3]

## Get the common genes in both matrices
common.genes <- intersect(rownames(mat.10x), rownames(cm.on))

## Inspect the common genes
length(common.genes)
head(common.genes)

## Subset both matrices to the common genes
cm.on <- cm.on[common.genes,]
mat.10x <- mat.10x[common.genes,]

## Check dimensions and inspect
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

## Calculate matrix cross correlation for 1000 random cells 
cells.xcor <- sample(colnames(mat.10x),1000)
xcor <- cor(as(mat.10x[,cells.xcor],'matrix'),as(cm.on[,cells.xcor],'matrix'),method='spearman')

## Check if the maximal matrix cross-correlation always occurs for the same cells
sum(apply(xcor, 1, which.max) == seq_along(rownames(xcor)))

## Generate a non-rearranged heatmap for these cross-correlations
pheatmap(xcor,cluster_rows=FALSE, cluster_cols=FALSE, file='mat.img/xcor.heatmap.png', show_rownames=FALSE, show_colnames=FALSE)

## Generate pagoda2 apps for these samples and examine manually
p2.cr <- basicP2proc(mat.10x[,colSums(mat.10x) > 10],n.cores=16)
p2.cr.web <- basicP2web(p2.cr,app.title='10x',n.cores=16)
p2.opt <- basicP2proc(cm.on[,colSums(cm.on) > 10],n.cores=16)
p2.opt.web <- basicP2web(p2.opt, app.title='opt',n.cores=16)

## Show the apps (this will lauch a Rook server instance)
show.app(p2.cr.web,browse=FALSE,name='cr')
show.app(p2.opt.web,browse=FALSE,name='opt')

## Get the count correlation excluding 0s
df1 <- data.frame(cr=log10(rowSums(mat.10x) + 1), opt=log10(rowSums(cm.on) +1))
df1 <- df1[df1$cr > 1e-6& df1$opt > 1e-6,]
cor(df1$cr ,df1$opt)

## Get the delta matrix masking out elements that are 0s in both matrices
mask <- mat.10x + cm.on == 0
delta.mat <- as.matrix(mat.10x - cm.on)
mask <- as.matrix(mask)
delta.mat[mask] <- c(NA)
rs <- log10(rowSums(abs(delta.mat),na.rm=TRUE)+1)
summary(rs)
dim(delta.mat)
## This is per element mean change (calculated from above summary output)
10^0.6021/4091

## Per Gene Sum of absolute per element differences between Optimus and Cellrange
png('mat.img/delta.hist.3.png')
hist(rs,breaks=50,xlab="Per Gene Sum of absolute per element differences between Optimus and Cellranger")
dev.off()

## Optionally save the final image
preserve.state(n.cores=16)
