library(pagoda2);
library(Matrix)

umi.on.path <- '/home/nbarkas/disk2/umi-check/8k_pbmc/MergeCountFiles/call-MergeCountFiles/convert/out.rds'
umi.off.path <- '/home/nbarkas/disk2/umi-check/8k_pbmc/MergeSorted/manualcount/convert/out.rds'

cm.on <- readRDS(umi.on.path)
cm.off <- readRDS(umi.off.path)

dim(cm.on)
dim(cm.off)

nrow(cm.on) - nrow(cm.off)

common_row <- intersect(rownames(cm.on), rownames(cm.off))

cm.on.cr <- cm.on[common_row,]
cm.off.cr <- cm.off[common_row,]

cells <- rownames(cm.on.cr)[rowSums(cm.on.cr) > 1000]
length(cells)

cm.on.cr <- cm.on.cr[cells,]
cm.off.cr <- cm.off.cr[cells,]

all(rownames(cm.on.cr) == rownames(cm.off.cr))

delta.mat <- cm.on.cr - cm.off.cr

png('delta_hist.png')
hist(as.vector(as.matrix(delta.mat)),breaks=100)
dev.off()

png('delta_hist_nz.png')
tmp <- as.vector(as.matrix(delta.mat))
tmp <- tmp[tmp!=0]
hist(tmp,breaks=100)
dev.off()


library(pagoda2)
cm.on.cr <- as(t(cm.on.cr), 'dgCMatrix')
p2.on <- basicP2proc(cm.on.cr,n.cores=24)
p2.on.web <- basicP2web(p2.on)

cm.off.cr <- as(t(cm.off.cr), 'dgCMatrix')
p2.off <- basicP2proc(cm.off.cr,n.cores=12)
p2.off.web <- basicP2web(p2.off)

p2.on.web$serializeToStaticFast('8k_pbmc_1_onUMI.bin')
p2.off.web$serializeToStaticFast('8k_pbmc_1_offUMI.bin')

show.app(p2.on.web,'on',browse=F)
show.app(p2.off.web,'off',browse=F)

## The Cellranger processed version of the datasets
library(nbHelpers)
cr.path <- "/home/nbarkas/disk2/umi-check/8k_pbmc/10Xproc/raw_gene_bc_matrices/GRCh38"
mat.10x <- read10xMatrix(cr.path)
mat.10x[1:3,1:3]
colnames(mat.10x) <- strpart(colnames(mat.10x),'-',1)
mat.10x[1:3,1:3]
mat.10x <- mat.10x[,cells]
dim(mat.10x)
mat.10x <- mat.10x[rowSums(mat.10x) > 0,]
## one cell doesn't have counts in 10x
mat.10x <- mat.10x[,cells]
mat.10x <- mat.10x[,colSums(mat.10x) > 0]
dim(mat.10x)

p2.10x <- basicP2proc(mat.10x,n.cores=12)
p2.10x.web <- basicP2web(p2.10x, n.cores=10)
show.app(p2.10x.web,name='CR',browse=F)
p2.10x.web$serializeToStaticFast('8k_pbmb_cellranger.bin')

cm.on.cr[1:3,1:3]

df.on <- pagoda2:::colMeanVarS(p2.on$counts, NULL, 16)
df.off <- pagoda2:::colMeanVarS(p2.off$counts, NULL, 16)

all(colnames(p2.on$counts) == colnames(p2.off$counts))

psc <- 1e-8

png('img/var.umi.on.off.png')
plot(log10(df.on$v+psc), log10(df.off$v+psc),xlab='log10(Per Gene Variance with UMI Correction)',
     ylab='log10(Per Gene Variance without UMI Correction')
dev.off()

