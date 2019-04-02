library(pagoda2)
library(nbHelpers)

umi.on.path <- "/home/nbarkas/disk2/umi-check/MantonBM1_HiSeq_1/MergeCountFiles/convert/out.rds"
umi.off.path <- "/home/nbarkas/disk2/umi-check/MantonBM1_HiSeq_1/MergeSorted/manualcount/convert/out.rds"

cm.on <- readRDS(umi.on.path)
cm.off <- readRDS(umi.off.path)

dim(cm.on)
dim(cm.off)

cells.keep <- rownames(cm.on)[rowSums(cm.on) > 1000]

cm.on.f <- cm.on[cells.keep,]
cm.off.f <- cm.off[cells.keep,]

cm.on.f.t <- as(t(as(cm.on.f,'matrix')),'dgCMatrix')
cm.off.f.t <- as(t(as(cm.off.f,'matrix')),'dgCMatrix')

p2.on <- basicP2proc(cm.on.f.t,n.cores=12)
p2.off <- basicP2proc(cm.off.f.t,n.cores=12)

p2.w.on <- basicP2web(p2.on,n.cores=4)
p2.w.off <- basicP2web(p2.off, n.cores=4)

p2.w.on$serializeToStaticFast('BM1_1_onUMI.bin')
p2.w.off$serializeToStaticFast('BM1_1_offUMI.bin')

show.app(p2.w.on, 'on', browse=FALSE)
show.app(p2.w.off, 'off', browse=FALSE)


##

cr.path <- "/home/nbarkas/disk3/MantonBM1_1/run/MantonBM1_HiSeq_1/outs/raw_feature_bc_matrix/"
cm.cr <- read10xMatrix(cr.path, version='V3')

head(colnames(cm.cr))

colnames(cm.cr) <- strpart(colnames(cm.cr), '-', 1)

cm.cr.f <- cm.cr[,cells.keep]
cm.cr.f <- cm.cr.f[rowSums(cm.cr.f) > 10, colSums(cm.cr.f) > 10]
class(cm.cr.f)

p2.cr <- basicP2proc(cm.cr.f)
p2.w.cr <- basicP2web(p2.cr)
p2.w.cr$serializeToStaticFast('BM1_1_cr.bin')

show.app(p2.w.cr, 'cellranger', browse=FALSE)

str1(p2.on)

str1(as.list(p2.on@.xData))

##

df.on <- pagoda2:::colMeanVarS(p2.on$counts, NULL, 16)
df.off <- pagoda2:::colMeanVarS(p2.off$counts, NULL, 16)

all(colnames(p2.on$counts) == colnames(p2.off$counts))

psc <- 1e-8
png('img/var.umi.on.off.png')
plot(log10(df.on$v+psc), log10(df.off$v+psc),xlab='log10(Per Gene Variance with UMI Correction)',
     ylab='log10(Per Gene Variance without UMI Correction')
dev.off()


## Cluster to cluster plot to show stability
cl.on <- p2.on$clusters$PCA$multilevel
cl.off <- p2.off$clusters$PCA$multilevel

all(names(cl.on) == names(cl.off))

df1 <- data.frame(as.character(cl.on), as.character(cl.off))
head(df1)

library(reshape2)
install.packages('pheatmap')
library(pheatmap)

x <- acast(df1, cl.on ~ cl.off,fun.aggregate=length)
x
x <- sweep(x, 2, apply(x, 2, sum), FUN='/')

pheatmap(log10(x+1), filename='img/cluster.correspondance.png')

## library(devtools)
## install_github('barkasn/fastSave')

library('fastSave')
preserve.state()

load.lbzip2('savepoint_2019-03-19_15:40:10_2447.RDataFS')


###
png('img/var.umi.on.off.png')
plot(log10(df.on$v+psc), log10(df.off$v+psc),xlab='log10(Per Gene Variance with UMI Correction)',
     ylab='log10(Per Gene Variance without UMI Correction')
dev.off()


keep.rows <- df.on$v > 0 & df.off$v >0
table(keep.rows)

delta.var <- (log10(df.on$v+psc) - log10(df.off$v+psc))[keep.rows]

head(df.on)

png('img/deltavar.vs.exp.png')
plot( log10((df.on$m + df.off$m)/2+psc),(log10(df.on$v+psc) - log10(df.off$v+psc)))
abline(0,0,col='red')
dev.off()

png('img/hist.delta.var.png')
hist(delta.var,breaks=50,xlab='Var(UMI_ON)-Var(UMI_OFF)')
dev.off()

getwd()

table((log10(df.on$v+psc) - log10(df.off$v+psc) > 0.3))

x <- colnames(p2.on$counts)[log10(df.on$v+psc) - log10(df.off$v+psc) > 0.3]
paste0('^',paste(x,collapse='$|^'),'$')


show.app(p2.w.cr,name='cr',browse=FALSE)

head(x)
e(x %in% colnames(p2.on$misc$rawCounts))

summary(colSums(p2.on$misc$rawCounts[,x]))
