## This is a template script for generating pagoda2 web apps that were used to inspect the
## datasets and generate some of the figures

## Load the
install.github('hms-dbmi/pagoda2')
library(pagoda2)
library('Matrix')

## This is the name of the app (e.g. "4k_pbmc")
appname <- 'appname'

## Load the count matrix, this has been converted from the Optimus output
## to the an rds file using the npz2rds utility (can be found in skylab github
## repository).
cm.raw <- readRDS('../convRDS/counts.rds')

## Load emptydrops output table to identify
## which cells to retain
ed.out <- read.csv('../RunEmptyDrops/empty_drops_result.csv')

## Get only cells called as true cells
called.cells.def <- ed.out$CellId[ed.out$IsCell]
length(called.cells.def)

## Subset matrix so that it can be processed by
## pagoda2.

## Keep only cells identified by emptyDrops
cm <- cm.raw[called.cells.def,]

## Keep only genes with at least 10 reads
cm <- cm[,colSums(cm) > 10]

## Filter any cells with fewer than 1000 reads
cm <- cm[rowSums(cm) > 1000,]

## Transpose, while retaining correct class
cm <- as(t(as(cm,'matrix')),'dgCMatrix')

## Finally, process with pagoda2
p2 <- basicP2proc(cm, n.cores=8)
p2w <- basicP2web(p2,app.title=appname,n.cores=8)

## Save as a binary file that can be inspected later
p2w$serializeToStaticFast(paste0(appname,'.bin'))

## Load the app in a web server
##show.app(p2w, name=appname, browse=FALSE)                          
