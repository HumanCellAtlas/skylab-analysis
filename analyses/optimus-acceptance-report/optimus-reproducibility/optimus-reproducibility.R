## This script compares three runs from Optimus to check the reproducibility of the pipeline
## For the comparison to be valid the three runs must be performed on identical inputs and caching
## must be turned off

## Load libraries
library(Matrix)

## Helper Functions
NA2VALUE <- function(x, val) {
    x[is.na(x)] <- c(val); x
}

NA2FALSE <- function(x) {
    NA2VALUE(x,FALSE)
}

## Load the datasets

## run1
run1.cm.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run1/MergeCountFiles/convert/out.rds'
run1.emptyDrops.table.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run1/RunEmptyDrops/empty_drops_result.csv'
run1.cm <- readRDS(run1.cm.path)
dim(run1.cm)
run1.ed <- read.csv(run1.emptyDrops.table.path)
head(run1.ed)

## run2
run2.cm.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run2/MergeCountFiles/convert/out.rds'
run2.emptyDrops.table.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run2/runEmptyDrops/empty_drops_result.csv'
run2.cm <- readRDS(run2.cm.path)
dim(run2.cm)
run2.ed <- read.csv(run2.emptyDrops.table.path)
head(run2.ed)

## run3
run3.cm.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run3/MergeCountFiles/convert/out.rds'
run3.emptyDrops.table.path <- '/home/nbarkas/disk2/optimus-reproducibility/data/run3/RunEmptyDrops/empty_drops_result.csv'
run3.cm <- readRDS(run3.cm.path)
dim(run3.cm)
run3.ed <- read.csv(run3.emptyDrops.table.path)
head(run3.ed)

## Perform comparisons of the count matrices
## As the matrices are identical its easier to just compare the compressed sparse array
## representations directly

## compare 1 and 2
dim(run1.cm)
dim(run2.cm)
all(rownames(run1.cm) == rownames(run2.cm))
all(colnames(run1.cm) == colnames(run2.cm))
str(run1.cm)
all(run1.cm@x == run2.cm@x)
all(run1.cm@p == run2.cm@p)
all(run1.cm@j == run2.cm@j)

## compare 1 and 3
dim(run1.cm)
dim(run3.cm)
all(rownames(run1.cm) == rownames(run3.cm))
all(colnames(run1.cm) == colnames(run3.cm))
str(run1.cm)
all(run1.cm@x == run3.cm@x)
all(run1.cm@p == run3.cm@p)
all(run1.cm@j == run3.cm@j)


## Perform comparisons of the emptyDrops output matrices

## run1 vs run2
head(run1.ed)
all(run1.ed$CellId == run2.ed$CellId)
all(run1.ed$IsCell == run2.ed$IsCell)
table(run1.ed$IsCell == run2.ed$IsCell)

## run2 vs run3
all(run2.ed$CellId == run3.ed$CellId)
all(run2.ed$IsCell == run3.ed$IsCell)
table(run2.ed$IsCell == run3.ed$IsCell)

## run1 vs run3
all(run1.ed$CellId == run3.ed$CellId)
all(run1.ed$IsCell == run3.ed$IsCell)
table(run1.ed$IsCell == run3.ed$IsCell)

## Set the row names to be cell ids
rownames(run1.ed) <- run1.ed$CellId
rownames(run2.ed) <- run2.ed$CellId
rownames(run3.ed) <- run3.ed$CellId

## Run1 vs Run2 -- what are the cells that are different?
diff.cells <- as.character(run1.ed$CellId[NA2FALSE(run1.ed$IsCell) != NA2FALSE(run2.ed$IsCell)])
table(run1.ed[diff.cells,]$IsCell)
table(run2.ed[diff.cells,]$IsCell)
length(diff.cells)
summary(run1.ed[diff.cells,]$Total)

## Run2 vs Run3
diff.cells <- as.character(run2.ed$CellId[NA2FALSE(run2.ed$IsCell) != NA2FALSE(run3.ed$IsCell)])
length(diff.cells)
table(run2.ed[diff.cells,]$IsCell)
table(run3.ed[diff.cells,]$IsCell)
summary(run2.ed[diff.cells,]$Total)

## Run1 vs Run3
diff.cells <- as.character(run1.ed$CellId[NA2FALSE(run1.ed$IsCell) != NA2FALSE(run3.ed$IsCell)])
length(diff.cells)
table(run1.ed[diff.cells,]$IsCell)
table(run3.ed[diff.cells,]$IsCell)
summary(run1.ed[diff.cells,]$Total)

