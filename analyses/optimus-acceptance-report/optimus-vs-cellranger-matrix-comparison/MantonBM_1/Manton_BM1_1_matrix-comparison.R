## Comparison analysis of the cellranger and optimus output matrices
## for the Ischaemia 1 dataset
library(Matrix)
library(nbHelpers)
library(pheatmap)
library(pagoda2)
library(fastSave)

## Load the data

## this is the direct 10X output matrix
cr.path <- "/home/nbarkas/disk3/MantonBM1_1/run/MantonBM1_HiSeq_1/outs/raw_feature_bc_matrix/"
mat.10x <- read10xMatrix(cr.path,version='V3')

## This is the optimus matrix from the same run with converted to rds format                                                                    
## with npz2rds utility that can be found under the skylab repository under                                                                     
## /docker/emptydrops/npz2rds/   
umi.on.path <- '/home/nbarkas/disk2/umi-check/MantonBM1_HiSeq_1/MergeCountFiles/convert/out.rds'
cm.on <- readRDS(umi.on.path)

## Fix the 10X names that have an -1 appended at the end of the cell
colnames(mat.10x) <- strpart(colnames(mat.10x), '-', 1)

## Load the emptyDrops output on the optimus
ed.path <- "/home/nbarkas/disk2/check-empty-drops/MantonBM1_1/RunEmptyDrops/empty_drops_result.csv"
ed.table <- read.csv(ed.path)
head(ed.table)

## Update the cutoffs 
ed.table$IsCell <- ed.table$FDR < 0.01 & ed.table$Total > 100
keep.cells <- as.character(ed.table$CellId[ed.table$IsCell])
head(keep.cells)
length(keep.cells)
dim(mat.10x)
dim(cm.on)

## Inspect
mat.10x[1:3,1:3]
cm.on[1:3,1:3]

## Subset both to keep cells from emptyDrops
mat.10x <- mat.10x[,keep.cells]
cm.on <- cm.on[keep.cells,]

## Tranpose the optimus matrix
cm.on <- t(as(cm.on, 'dgCMatrix'))

## Check dims and inspect
dim(mat.10x)
dim(cm.on)
cm.on[1:3,1:3]
mat.10x[1:3,1:3]

## Get the genes common in both matrices
common.genes <- intersect(rownames(mat.10x), rownames(cm.on))
length(common.genes)
head(common.genes)

## Subset to the genes common in both matrices
cm.on <- cm.on[common.genes,]
mat.10x <- mat.10x[common.genes,]

## Check dims and inspect
dim(mat.10x)
dim(cm.on)
mat.10x[1:10,1:10]
cm.on[1:10,1:10]

## Generate the plots
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

## Get the gens that are not expressed at all in both matrices
allzerogenes <- rowSums(mat.10x) == 0 & rowSums(cm.on) ==0
table(allzerogenes)

## Make a plot exclusing these genes
png('mat.img/perGeneAbsDeltaHist.filter0.png')
hist(log10(rowSums(abs(mat.10x - cm.on)[!allzerogenes,])+1),breaks=100,
     main='Histogram of per gene absolute matrix differences',
     xlab='log10(perGene sum of differences)')
dev.off()

## Calculte the cross-correlation for 1000 random cells
cells.xcor <- sample(colnames(mat.10x),1000)
xcor <- cor(as(mat.10x[,cells.xcor],'matrix'),as(cm.on[,cells.xcor],'matrix'),method='spearman')

## Check that every cell maximally x-correlates with itself
sum(apply(xcor, 1, which.max) == seq_along(rownames(xcor)))

## Plot a heatmap with out col or row reordering
pheatmap(xcor,cluster_rows=FALSE, cluster_cols=FALSE, file='mat.img/xcor.heatmap.png', show_rownames=FALSE, show_colnames=FALSE)

## Generate pagoda2 apps for manual inspection of hte dataset
p2.cr <- basicP2proc(mat.10x[,colSums(mat.10x) > 10],n.cores=16)
p2.cr.web <- basicP2web(p2.cr,app.title='10x',n.cores=16)

p2.opt <- basicP2proc(cm.on[,colSums(cm.on) > 10],n.cores=16)
p2.opt.web <- basicP2web(p2.opt, app.title='opt',n.cores=16)

## display the apps -- note that this will start a Rook web server
show.app(p2.cr.web,browse=FALSE,name='cr')
show.app(p2.opt.web,browse=FALSE,name='opt')

## Histogram removing genes that are 0 in both 
df1 <- data.frame(cr=log10(rowSums(mat.10x) + 1), opt=log10(rowSums(cm.on) +1))
df1 <- df1[df1$cr > 1e-6& df1$opt > 1e-6,]

## Calculate the correlation for the above plot
cor(df1$cr ,df1$opt)

## get delta mat masking out elements that were 0s in both cases
mask <- mat.10x + cm.on == 0
delta.mat <- as.matrix(mat.10x - cm.on)
mask <- as.matrix(mask)
delta.mat[mask] <- c(NA)
rs <- log10(rowSums(abs(delta.mat),na.rm=TRUE)+1)
summary(rs)
dim(delta.mat)
## the difference on a per cell basis
10^0.5077/4766

## Plot Per Gene Sum of absolute per element differences between Optimus and Cellranger
png('mat.img/delta.hist.3.png')
hist(rs,breaks=50,xlab="Per Gene Sum of absolute per element differences between Optimus and Cellranger")
dev.off()

## Count all zero gens
allzerogenes <- rowSums(mat.10x) == 0 & rowSums(cm.on) ==0
table(allzerogenes)

## Same plot with all zero genes removed
png('mat.img/perGeneAbsDeltaHist.filter0.png')
hist(log10(rowSums(abs(mat.10x - cm.on)[!allzerogenes,])+1),breaks=100,
     main='Histogram of per gene absolute matrix differences',
     xlab='log10(perGene sum of differences)')
dev.off()

## correlation coefficient for scatter above
df1 <- data.frame(cr=log10(rowSums(mat.10x) + 1), opt=log10(rowSums(cm.on) +1))
df1 <- df1[df1$cr > 1e-6& df1$opt > 1e-6,]
cor(df1$cr ,df1$opt)

## Optionally save the final state
preserve.state(n.cores=16)
