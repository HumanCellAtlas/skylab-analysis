umi.on.path <- '/home/nbarkas/disk2/umi-check/ischaemia_1/MergeCountFiles/convert/out.rds'
umi.off.path <- '/home/nbarkas/disk2/umi-check/ischaemia_1/MergeSorted/manualcount/convert/out.rds'

cm.on <- readRDS(umi.on.path)
cm.off <- readRDS(umi.off.path)

dim(cm.on)
dim(cm.off)

nrow(cm.on) - nrow(cm.off)

common <- row <- intersect(rownames(cm.on), rownames(cm.off))

cm.on.cr <- cm.on[common <- row,]
cm.off.cr <- cm.off[common <- row,]

table(as.vector(cm.on.cr - cm.off.cr) == 0)

cells <- rownames(cm.on.cr)[rowSums(cm.on.cr) > 1000]
length(cells)

cm.on.cr <- cm.on.cr[cells,]
cm.off.cr <- cm.off.cr[cells,]

dim(cm.on.cr)
dim(cm.off.cr)

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
p2.on <- basicP2proc(cm.on.cr,n.cores=8)
p2.on.web <- basicP2web(p2.on)

cm.off.cr <- as(t(cm.off.cr), 'dgCMatrix')
p2.off <- basicP2proc(cm.off.cr,n.cores=8)
p2.off.web <- basicP2web(p2.off)

show.app(p2.on.web,name='on',browse=FALSE)
show.app(p2.off.web,name='off',browse=FALSE)

p2.on.web$serializeToStaticFast('4k_pbmc.umi-on.bin')
p2.off.web$serializeToStaticFast('4k_pbmc.umi-off.bin')

library(fastSave)
preserve.state()

df.on <- pagoda2:::colMeanVarS(p2.on$counts, NULL, 16)
df.off <- pagoda2:::colMeanVarS(p2.off$counts, NULL, 16)

all(colnames(p2.on$counts) == colnames(p2.off$counts))

psc <- 1e-8
png('img/var.umi.on.off.png')
plot(log10(df.on$v+psc), log10(df.off$v+psc),xlab='log10(Per Gene Variance with UMI Correction)',
     ylab='log10(Per Gene Variance without UMI Correction')
dev.off()
