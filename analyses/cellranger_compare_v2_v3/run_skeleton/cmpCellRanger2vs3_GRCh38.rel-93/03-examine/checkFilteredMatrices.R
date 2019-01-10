## Name: checkFilteredMatrices.R
## Author: Nikolas Barkas
## Date: December 2018
## Description: Generate the plots for the main part of the report by comparing the V2V2 and V3V3 matrices

## Load required libraries
library(nbHelpers)
library(fastSave)
library(Matrix)
library('VennDiagram')
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pagoda2)

## Load the auxilary functions
source('functions.R')

## The V2 and V3 filtered input files
v2filt.file <- '/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV2_annotV2/pbmc_4k/outs/filtered_gene_bc_matrices/GRCh38'
v3filt.file <- '/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV3/pbmc_4k/outs/filtered_feature_bc_matrix'

## mixed V3 exec with V2 build ref
v2v3filt.file <- '/data/run_skeleton/cmpCellRanger2vs3_GRCh38.rel-93/02-count/execV3_annotV2/pbmc_4k/outs/filtered_feature_bc_matrix'

## Read in the input data
v2.mat <- read10xMatrix(v2filt.file, 'V2')
v3.mat <- read10xMatrix(v3filt.file, 'V3')

## Check the V2 V3 vs the other two
v3v2.mat <- read10xMatrix(v2v3filt.file, 'V3')
dim(v2.mat)
dim(v3.mat)
dim(v3v2.mat)
v3_v3v2.diff <- drop0(v3.mat - v3v2.mat)

## Check total number of UMIs
sum(v2.mat)
sum(v3.mat)
sum(v3v2.mat)

## Quickly checking for non zero differences
## Any different entries after drop0() will be in the @x vector
v3_v3v2.diff@x
str(v3_v3v2.diff)
## x is empty, therefore identical to v3

## Genes and order are identical
table(rownames(v2.mat) == rownames(v3.mat))

## There are differences in the number of cells
ncol(v2.mat) ## 4325
ncol(v3.mat) ## 4600

## Let's focus on the cells in common
common.cells <- intersect(colnames(v2.mat), colnames(v3.mat))
length(common.cells) ## All the cells from V2
## V3 are a superset of V3

## Let's make a Venn of the cell overlap
pal <- colorspace::rainbow_hcl(n=2)

## Plot the venn
png('outs/cell.overlap.venn.png', width=600, height=600)
draw.pairwise.venn(ncol(v2.mat),ncol(v3.mat),length(common.cells),euler.d=FALSE,scaled=FALSE, fill=pal,
                   cex=3)
dev.off()

## Subset matrices to the relavent cells only
v2.mat.c <- v2.mat[,common.cells]
v3.mat.c <- v3.mat[,common.cells]

## Double check the dimentions are identical
all(rownames(v2.mat.c) == rownames(v3.mat.c))
all(colnames(v2.mat.c) == colnames(v3.mat.c))

## Make a matrix of the differences
diff.mat <- as.matrix(v2.mat.c - v3.mat.c)

## How many entries of the matrices are identical?
identical.count <- sum(as.matrix(diff.mat) == 0)
identical.percentage <- identical.count / length(diff.mat)
identical.percentage ## 0.996404

## How many increase and how many decrease?
sum(diff.mat < 0) ## Values are only lower in V3, never higher!!!
sum(diff.mat > 0) ## 524,869 (these are higher in v2), this is a very small percentage of the entire matrix
## Clearly they are just discarding reads in V3, which is consistend with the extra filtering step
## on chimeric reads

X11()

## Let's look at the matrix now from the gene perspective

## How many entries are different per gene
diff.entries.per.gene <- apply(diff.mat > 0, 1, sum)

## Mean expression per gene in v2 vs v3
mean.diff.per.gene <- apply(diff.mat, 1, mean)

## mean expression per gene in V2 and V3
v2.mean.expr <- apply(v2.mat.c,1,mean)
v3.mean.expr <- apply(v3.mat.c,1, mean)

## expression variance for V2 and V3
v2.expr.var <- apply(v2.mat.c, 1, var)
v3.expr.var <- apply(v3.mat.c, 1, var)

## Mean Gene difference vs Mean Gene Expression
## in a single data.frame
gene.diff.summary <- data.frame(
    gene.id=names(v2.expr.var),
    v2.mean.expr,
    v3.mean.expr,
    mean.diff.per.gene,
    diff.entries.per.gene,
    v2.expr.var,
    v3.expr.var
)

## Histogram of gene sums
png('outs/genesum_dist.png',width=480*2,height=480)
par(mfrow=c(1,2))
breaks <- 25
hist(log10(apply(v2.mat.c,1,sum)+1),
     breaks=breaks,
     main='V2V2',
     cex=2,
     xlab='log10(genesum + 1) for V2V2')
hist(log10(apply(v3.mat.c,1,sum)+1),
     breaks=breaks,
     main='V3V3',
     cex=2,
     xlab='log10(genesum + 1) for V3V3')
dev.off()

## Test with KS test
set.seed(2018)
ks.test(
    sample(log10(apply(v2.mat.c,1,sum)+1),5000),
    sample(log10(apply(v3.mat.c,1,sum)+1),5000)
)

ks.test(
    log10(apply(v2.mat.c,1,sum)+1),
    log10(apply(v3.mat.c,1,sum)+1),
)

## Histogram of gene sums
png('outs/cellsum_dist.png',width=480*2,height=480)
par(mfrow=c(1,2))
breaks <- 25
hist(log10(apply(v2.mat.c,2,sum)+1),
     breaks=breaks,
     main='V2V2',
     cex=2,
     xlab='log10(cellsum + 1) for V2V2')
hist(log10(apply(v3.mat.c,2,sum)+1),
     breaks=breaks,
     main='V3V3',
     cex=2,
     xlab='log10(cellsum + 1) for V3V3')
dev.off()


ks.test(
    log10(apply(v2.mat.c,2,sum)+1),
    log10(apply(v3.mat.c,2,sum)+1),
)


## Check how many genes are affected and how many are not
gene.diff.summary$has.change <- apply(diff.mat > 0, 1, sum) > 0

table(gene.diff.summary$has.change)

custom_theme <- theme_bw() + theme(axis.text=element_text(size=16,color='black'),
                                   axis.title=element_text(size=16,face="bold",color='black'),
                                   title =element_text(size=12, face='bold'))

## Plot the number of genes that have at least one value changed
ggplot(gene.diff.summary, aes(x=has.change)) + geom_bar() +
    ggtitle('Number of Genes with at least one change') +
    custom_theme + scale_x_discrete(name='Gene has > 1 change') +
    scale_y_continuous(name='Number of genes')
ggsave('outs/num.genes.with.change.png',w=7,h=7)

## Plot the number of changes per gene
ggplot(gene.diff.summary, aes(x=diff.entries.per.gene)) + geom_bar() + scale_y_log10() +
    ggtitle('Number of Entries that different per gene') +
    custom_theme
ggsave('outs/num.entries.diff.pergene.png',w=7,h=7)

## Plot the number of changes per gene
ggplot(gene.diff.summary, aes(y=log10(diff.entries.per.gene+1),x="All genes")) + geom_boxplot() +
    ggtitle('Number of Entries different per gene') +
    scale_y_continuous(name=expression(log[10]("Diff Entry Count" + 1))) + 
    custom_theme + scale_x_discrete(name='')
ggsave('outs/hist.num.genes.with.change.png',w=7,h=7)

## V2 vs V3 expr
ggplot(gene.diff.summary, aes(x=v2.mean.expr, y=v3.mean.expr)) + geom_point() +
    ggtitle('Mean Expression in V2V2 vs V3V3') + scale_x_continuous(name='V2V2 mean expression') +
    scale_y_continuous(name='V3V3 mean expression') + 
    geom_abline(slope=1,intercept=0, color='red') + coord_fixed() +
    custom_theme
ggsave('outs/v2.vs.v3.expression.linear.png',w=7,h=7)

## Correlation of log10 values
with(gene.diff.summary,
     cor(log10(v2.mean.expr+1), log10(v3.mean.expr+1)))
## 0.9994629

with(gene.diff.summary,
     cor(log10(v2.mean.expr+1), log10(v3.mean.expr+1), method='spearman'))
## 0.9994629

## Log10 V2 vs V3 expr
ggplot(gene.diff.summary, aes(x=log10(v2.mean.expr+1), y=log10(v3.mean.expr+1))) + geom_point() +
    ggtitle('Mean Expression in V2V2 vs V3V3') +
    scale_x_continuous(name=expression(log10("V2V2 mean expression" + 1))) +
    scale_y_continuous(name=expression(log10("V3V3 mean expression" + 1))) +
    geom_abline(slope=1,intercept=0, color='red') + coord_fixed() +
    custom_theme
ggsave('outs/v2.vs.v3.expression.log10.png',w=7,h=7)

## Rank V2 vs V3 expr
gene.diff.summary$v2.rank <- rank(v2.mean.expr)
gene.diff.summary$v3.rank <- rank(v3.mean.expr)

## Rank correlation
with(gene.diff.summary,
          cor(v2.rank, v3.rank))
## 0.9957007

## a temp dataframe for the datapoint to label
tmp1 <- subset(gene.diff.summary, #(v2.rank - v3.rank) > 100 &
                                  ((v2.rank - v3.rank) / ((v2.rank + v3.rank)/2)) > 0.05 &
                                  v2.rank > 20000 & v3.rank > 20000)
nrow(tmp1)




p <- ggplot(gene.diff.summary, aes(x=v2.rank, y=v3.rank)) + geom_point() +
    ggtitle('Rank comparison of Mean Expression in V2V2 vs V3V3') +
    scale_x_continuous(name='Rank of V2V2 mean expression') +
    scale_y_continuous(name='Rank of V3V3 mean expression') +
    geom_label_repel(aes(x=v2.rank,y=v3.rank,label=gene.id, color='red'), data=tmp1) +
    custom_theme
ggsave('outs/rank.comparison.labelled.png', p, width=7, height=7)
rm(p)
    
## log10 mean diff per gene vs gene rank in V2
#ggplot(gene.diff.summary, aes(x=rank(v2.mean.expr), y=log10(mean.diff.per.gene+1))) + geom_point()

## Mean Variance for V2 and V3
p.v2 <- ggplot(gene.diff.summary, aes(v2.mean.expr, v2.expr.var)) + geom_point() + 
    ggtitle('V2V2 mean vs variance') + scale_x_continuous(name='Mean Expression') +
    scale_y_continuous(name='Variance') + custom_theme
p.v3 <- ggplot(gene.diff.summary, aes(v3.mean.expr, v3.expr.var)) + geom_point() + 
    ggtitle('V3V3 mean vs variance') + scale_x_continuous(name='Mean Expression') +
    scale_y_continuous(name='Variance')  + custom_theme
## install.packages('gridExtra')

p <- grid.arrange(p.v2,p.v3)
ggsave('outs/mean.vs.var.cmp.png',p)
rm(p.v2, p.v3, p)

## log(expr) vs log(var) for V2 and V3
p.v2 <- ggplot(gene.diff.summary, aes(
                              x=log10(v2.expr.var+1),
                              y=log10(v2.mean.expr+1)
                              )) + geom_point() +
    scale_x_continuous(name=expression(log[10](expr +1))) +
    scale_y_continuous(name=expression(log[10](Var(expr)+1))) +
    ggtitle('A   V2V2') + custom_theme
p.v3 <- ggplot(gene.diff.summary, aes(
                              x=log10(v3.expr.var+1),
                              y=log10(v3.mean.expr+1)
                              )) + geom_point() +
    scale_x_continuous(name=expression(log[10](expr +1))) +
    scale_y_continuous(name=expression(log[10](Var(expr)+1))) +
    ggtitle('B   V3V3') + custom_theme
p.v3vsv2 <- ggplot(gene.diff.summary, aes(
                              x=log10(v3.expr.var+1) - log10(v2.expr.var+1), 
                              y=log10(v3.mean.expr+1) - log10(v2.mean.expr+1)
                              )) + geom_point() +
    scale_x_continuous(name=expression(log[10](expr +1))) +
    scale_y_continuous(name=expression(log[10](Var(expr)+1))) +
    ggtitle('C Difference between V3V3 and V2V2') + custom_theme
p <- grid.arrange(grobs=list(p.v2, p.v3, p.v3vsv2), ncol=2, nrow=2)
ggsave('outs/mean.vs.var.panel1.png',p)
rm(p.v2, p.v3, p.v3vsv2)

## Relationship between gene rank and change
p <- ggplot(gene.diff.summary, aes(y=rank(v3.mean.expr), x=has.change)) + geom_violin() +
    ggtitle('Ranks of genes Vs changes') +
    scale_x_discrete(name='Gene as at least one change') +
    scale_y_continuous(name='Rank of Gene Expression in V3V3') +
    custom_theme
ggsave('outs/changed.vs.rank.png',p,w=7,h=7)

## Look at the changed values directly

## Get matrix entries as a 1d vector
direct.value.compare <- cbind(v2.val=as.vector(as.matrix(v2.mat.c)), v3.val=as.vector(as.matrix(v3.mat.c)))
 
## Keep only entries that are non-zero in both
nz.rows <- rowSums(direct.value.compare) > 0
nz.direct.value.compare <- direct.value.compare[nz.rows,]

## Calculate the difference for these values
diff <- nz.direct.value.compare[,1] - nz.direct.value.compare[,2]

## Keep only values that have changed
changed <- nz.direct.value.compare[diff !=0, ]

## Calculate difference between the two values against their mean
rel.diff <- (changed[,1] - changed[,2])/rowMeans(changed)

## Plot relative difference vs mean expression
tmp1 <- cbind(log10(rowMeans(changed)+1), rel.diff)

png('outs/value.relative.diff.png')
plot(tmp1, main="Relative difference of Changed Values vs gene mean",
     xlab=expression(log[10](mean)),cex=2,cex.axis=2)
dev.off()

## Focus on genes with expresion over 100
## The relative difference is < 10%
png('outs/value.relative.diff.over100.png')
plot(tmp1[tmp1[,1] > 2,], main="Relative difference of Changed Values vs gene mean",
     xlab=expression(log[10](mean)),cex=2,cex.axis=2)
dev.off()

#########################
## Per gene summary
png('outs/genes.detected.VS.UMI_A.png',width=480,height=480)
plot(
    log10(apply(v2.mat.c, 2, sum)+1),
    log10(apply(v2.mat.c >0, 2, sum)+1),
    main='V2V2 cellranger',
    xlab=expression(log[10](UMI)),
    ylab=expression(log[10](Genes)),
    cex=2,
    pch='.'
)
dev.off()

png('outs/genes.detected.VS.UMI_B.png',width=480,height=480)
plot(
    log10(apply(v3.mat.c, 2, sum)+1),
    log10(apply(v3.mat.c >0, 2, sum)+1),
    main='V3V3 cellranger',
    pch='.',
    xlab=expression(log[10](UMI)),
    ylab=expression(log[10](Genes)),
    cex=2
)
dev.off()

png('outs/genes.detected.VS.UMI_C.png',width=480,height=480)
plot(
    log10(apply(v3.mat.c, 2, sum)+1) - log10(apply(v2.mat.c, 2, sum)+1),
    log10(apply(v3.mat.c >0, 2, sum)+1) - log10(apply(v2.mat.c >0, 2, sum)+1),
    main='V3V3 cellranger',
    pch='.',
    xlab=expression(log[10](UMI)),
    ylab=expression(log[10](Genes)),
    cex=2
)
dev.off()

## Let's process both to see what the data actually look like

## V2
p2.v2 <- basicP2proc(v2.mat.c,n.cores =4)
web.v2 <- basicP2web(p2.v2)
web.v2$serializeToStaticFast('pbmc_4k_cellrangerV2V2.bin')
## V3
p2.v3 <- basicP2proc(v3.mat.c,n.cores=4)
web.v3 <- basicP2web(p2.v3)
web.v3$serializeToStaticFast('pbmc_4k_cellrangerV3V3.bin')

## % change in number of genes detected
pc.bc.diff <- (colSums(v3.mat.c > 0) - colSums(v2.mat.c > 0))/colSums(v2.mat.c >0) * 100
setDisplay('11.0')
X11()

png('outs/pc.change.num.genes.png')
hist(pc.bc.diff,breaks=40,main='Histogram of % delta(number of genes) between V2V2 and V3V3',
     xlab='% change in number of detected genes')
dev.off()

## Matrix with the cells only in V3
v3.mat.extra <- v3.mat[,setdiff(colnames(v3.mat), colnames(v2.mat))]

## Compare gene dist of additional cells
png('outs/cmpDisNoGenesOfadditional248_Orig.png',width=480,height=480)
hist(colSums(v3.mat.c >= 1),breaks=seq(0,5000,by=250),xlim=c(0,5000),
     main='Histogram of number of genes detected in V3V3 and V2V2',
     xlab='#genes / cell')
dev.off()

png('outs/cmpDisNoGenesOfadditional248_Extra.png',width=480,height=480)
hist(colSums(v3.mat.extra >= 1),breaks=seq(0,5000,by=250),xlim=c(0,5000),
     main='Histogram of number of genes detected for additional cells',
     xlab='#genes / cell')
dev.off()

## Venn of overlap of barcodes with at lease one read in V2V2 and V3V3

## Let's make a Venn of the cell overlap
pal <- colorspace::rainbow_hcl(n=2)

## Manually draw the Venn
png('outs/raw.cell.overlap.venn.png', width=600, height=600)
draw.pairwise.venn(
    272315,
    279136,
    271858,
    euler.d=FALSE,scaled=FALSE, fill=pal, cex=3)
dev.off()

