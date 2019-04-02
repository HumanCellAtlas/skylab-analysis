
## Load the libraries used
library(pagoda2)
library(fastSave)
library(nbHelpers)
library(reshape2)
library(ggplot2)

## Load the matrices generated with and without UMI correction
umi.on.path <- "/home/nbarkas/disk2/umi-check/4k_pbmc/data/convert/out.rds"
umi.off.path <- "/home/nbarkas/disk2/umi-check/4k_pbmc/data_bam/uncorrect-umis/manualcount/convert/out.rds"
cm.on <- readRDS(umi.on.path)
cm.off <- readRDS(umi.off.path)

## Check the dimentions of the matrix
dim(cm.on)
dim(cm.off)
nrow(cm.on) - nrow(cm.off)

## Common cells
common_row <- intersect(rownames(cm.on), rownames(cm.off))
cm.on.cr <- cm.on[common_row,]
cm.off.cr <- cm.off[common_row,]

## Keep cells with over 1000 reads in the umi corrected matrix
cells <- rownames(cm.on.cr)[rowSums(cm.on.cr) > 1000]
length(cells)
cm.on.cr <- cm.on.cr[cells,]
cm.off.cr <- cm.off.cr[cells,]

## Check dimentions
dim(cm.on.cr)
dim(cm.off.cr)

## Check that the rownames are identically ordered
all(rownames(cm.on.cr) == rownames(cm.off.cr))

## Get the delat matrix
delta.mat <- cm.on.cr - cm.off.cr

## Get histogram of eh delta matrix values
png('delta_hist.png')
hist(as.vector(as.matrix(delta.mat)),breaks=100)
dev.off()

## Get historgram of the delta matrix values with the zero entries removed
png('delta_hist_nz.png')
tmp <- as.vector(as.matrix(delta.mat))
tmp <- tmp[tmp!=0]
hist(tmp,breaks=100)
dev.off()

## Generate pagoda2 apps for both umi corrected an non-corrected matrices
cm.on.cr <- as(t(cm.on.cr), 'dgCMatrix')
p2.on <- basicP2proc(cm.on.cr,n.cores=8)
p2.on.web <- basicP2web(p2.on)

cm.off.cr <- as(t(cm.off.cr), 'dgCMatrix')
p2.off <- basicP2proc(cm.off.cr,n.cores=8)
p2.off.web <- basicP2web(p2.off)

## View the matrices, this will start the rook web server
show.app(p2.on.web,name='on',browse=FALSE)
show.app(p2.off.web,name='off',browse=FALSE)

## Serialize the apps to binary files for future viewing
p2.on.web$serializeToStaticFast('4k_pbmc.umi-on.bin')
p2.off.web$serializeToStaticFast('4k_pbmc.umi-off.bin')

## Load the cellranger matrix for comparison
cr.path <- '/home/nbarkas/disk3/4k_pbmc/run/pbmc4k/outs/raw_feature_bc_matrix/'
cm.cr <- read10xMatrix(cr.path, version='V3')
## fix the names
colnames(cm.cr) <- strpart(colnames(cm.cr), '-', 1)
cells.keep <- cells

## Keep the cells we kept above and do an additional filtering to
## remove some small cells, that cause problems downstream (<10 cells)
cm.cr.f <- cm.cr[,cells.keep]
cm.cr.f <- cm.cr.f[rowSums(cm.cr.f) > 10, colSums(cm.cr.f) > 10]

## Generate a pagoda2 app for the cr analysis
p2.cr <- basicP2proc(cm.cr.f,n.cores=16)
p2.w.cr <- basicP2web(p2.cr,n.cores=16)
p2.w.cr$serializeToStaticFast('4kPBMC_1_cr.bin')

## Show all the aps
show.app(p2.on.web, name='on')
show.app(p2.off.web, name='off')
show.app(p2.w.cr, name='cr')

## Here we want to assess the strability of clustering to subsetting
## the fraction of cells to keep
f <- 0.7
cl.diff <- lapply(1:100, function(i) {
    cells.sample <- sample(colnames(cm.on.cr),floor(ncol(cm.on.cr) * f))
    p2.on.sample <- basicP2proc(cm.on.cr[,cells.sample],n.cores=16, get.tsne=F,make.geneknn=F,get.largevis=F)
    p2.off.sample <- basicP2proc(cm.off.cr[,cells.sample],n.cores=16, get.tsne=F,make.geneknn=F,get.largevis=F)
    nlevels(p2.on.sample$clusters$PCA$multilevel) - nlevels(p2.off.sample$clusters$PCA$multilevel)
})

## Melt and inspect
x <- melt(cl.diff)
head(x)
table(x$value)

## Plot the histogram of the number of clusters found
png('delta.clusters.png')
hist(x$value)
dev.off()

 ## Plot with ggplot
p <- ggplot(x, aes(x=value)) + geom_histogram(binwidth=1,color='black') + theme_bw() + scale_x_continuous(name='difference in number of clusters')
ggsave('delta.clusters.png',p)

## Look at the per-gene variance before and after the correction
## this a complex relationship and the way it varies is also complex
df.on <- pagoda2:::colMeanVarS(p2.on$counts, NULL, 16)
df.off <- pagoda2:::colMeanVarS(p2.off$counts, NULL, 16)
all(colnames(p2.on$counts) == colnames(p2.off$counts))

## Plot the log10 per gene variancs in scatter plot
psc <- 1e-8 # pseudocount
png('img/var.umi.on.off.png')
plot(log10(df.on$v+psc), log10(df.off$v+psc),xlab='log10(Per Gene Variance with UMI Correction)',
     ylab='log10(Per Gene Variance without UMI Correction')
dev.off()

## Optionally save the final state
preserve.state()
