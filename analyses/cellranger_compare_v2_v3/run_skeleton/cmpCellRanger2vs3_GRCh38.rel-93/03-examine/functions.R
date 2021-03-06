## Name: function.R
## Author: Nikolas Barkas
## Date: December 2018
## Description: Accessory functions for analysis. 


## Function to read 10X matrices (both V2 and V3)
## Derived from the barkasn/nbHelpers package function
## but enhanced to read V3 matrices
read10xMatrix <- function(path, version='V2') {
    if(version == 'V2') {
        unpackFunction <- I
        suffix <- ''
    } else if (version == 'V3') {
        unpackFunction <- gzfile
        suffix <- '.gz'
    } else {
        stop('Unknown file version!')
    }
    matrixFile <- paste0(path, '/matrix.mtx', suffix);
    if (version == 'V2') {
        genesFile <- paste0(path, '/genes.tsv', suffix);
    } else if (version == 'V3') {
        genesFile <- paste0(path, '/features.tsv', suffix);
    }
    barcodesFile <- paste0(path, '/barcodes.tsv', suffix);
    if (!file.exists(matrixFile)) { stop('Matrix file does not exist');  }
    if (!file.exists(genesFile)) { stop('Genes file does not exist'); }
    if (!file.exists(barcodesFile)) { stop('Barcodes file does not exist'); }
    x <- as(Matrix::readMM(unpackFunction(matrixFile)), 'dgCMatrix')
    genes <- read.table(unpackFunction(genesFile));
    rownames(x) <- genes[,2];
    barcodes <- read.table(unpackFunction(barcodesFile));
    colnames(x) <- barcodes[,1]
    invisible(x);
}
