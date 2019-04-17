# SCATTER PLOTS

# the list of samples 
samples=c("ischaemia_1", "pbmc4k", "pbmc8k", "t_4k")


# Scatter plots from the tables with each line as 
#  <corrlation> <tab> cell_ranger row count<tab> optimus row count 
# created by the python script compare_matrices.py 
# for the case where the data is not normalized 
for( sample in samples) {
    #read the table containing the correlations, and the total row-wise (barcode wise counts) from Cell Ranger and Optimus count matrices.
    tab1 <- read.csv(   paste(sample, "non-normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
     
    # remove the ones with 0 counts in any one of the matrices.
     x<-tab1[ tab1$optimus > 0,]
     x<-x[ x$cell_ranger > 0,]

    # for adj R square value
    fit<-lm(x[,2]~x[,3])

    # print the plot
    png(paste(sample, "x_y_non-normalized_50000_v2.png", sep="_"))
    plot( x[, 2], x[,3], log="xy",  xlab="Gene counts for barcodes in CellRanger", ylab="Gene counts for barcodes in Optimum", main=paste(sample, " the sum total gene counts in barcode by \n Cellranger and Optimus (non-normalized)", "\nadj. R square ", format(summary(fit)$adj.r.squared, digits=4), sep=":"))
    dev.off()
 }

# Scatter plots for the python script compare_matrices.py where the counts for each gene columns is normalized by the total counts
# specified with the options --normalize 
for( sample in samples) {
    tab1 <- read.csv(   paste(sample, "normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
    # remove the ones with 0 counts in any one of the matrices.
     x<-tab1[ tab1$optimus > 0,]
     x<-x[ x$cell_ranger >0,]
    fit<-lm(x[,2]~x[,3])
    png(paste(sample, "x_y_normalized_50000_v2.png", sep="_"))
    plot( x[, 2], x[,3], log="xy",  xlab="Gene counts for barcodes in CellRanger", ylab="Gene counts for barcodes in Optimum", main=paste(sample, "sum total gene counts in barcode by \n Cellranger and Optimus (normalized)", "\nadj. R square ", format(summary(fit)$adj.r.squared, digits=4), sep=":"))
    dev.off()
 }

# Scatter plots for the python script compare_matrices.py where the counts for each gene columns is normalized by the total counts
# specified with the options --normalize and also the columns of one of the matrices is randomly permuted, 
# achieved with the --shuffle option
for(sample in samples) {
    tab1 <- read.csv(   paste(sample, "normalized_shuffle_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
    # keep the ones with more than 10 total (across the genes) counts in any one of the matrices.
     x<-tab1[ tab1$optimus > 10,]
     x<-x[ x$cell_ranger > 10,]
    fit<-lm(x[,2]~x[,3])

    png(paste(sample, "x_y_min_10_normalized_50000_random_permutation_v2.png", sep="_"))
plot( x[, 2], x[,3], log="xy",  xlab="Gene counts for barcodes in CellRanger", ylab="Gene counts for barcodes in Optimum", main=paste(sample, "sum total gene counts (> 10) in barcode by \n Cellranger and Optimus (normalized)", "\nadj. R square and with one row vector randomly permuted ", format(summary(fit)$adj.r.squared, digits=4), sep=":"))
    dev.off()
 }

# Scatter plots for the python script compare_matrices.py where the counts for each gene columns is normalized by the total counts
# for the entire column # specified with the options --normalize and the pair of row vectors are randomly selected from the pair of matrices. Note that 
# each of the row vector of counts correspond to one barcode. The random pair selection is the done with the  --cross-correlation option
for(sample in samples) {
    tab1 <- read.csv(   paste(sample, "normalized_50000_cross_correl_v2.txt", sep="_"), header=TRUE, sep='\t')
    # keep the ones with more than 10 total (across the genes) counts in any one of the matrices.
     x<-tab1[ tab1$optimus > 10,]
     x<-x[ x$cell_ranger > 10,]
    fit<-lm(x[,2]~x[,3])

    png(paste(sample, "x_y_min_10_normalized_50000_cross_correlation_v2.png", sep="_"))
    plot( x[, 2], x[,3], log="xy",  xlab="Gene counts for barcodes in CellRanger", ylab="Gene counts for barcodes in Optimum", main=paste(sample, "for sum total gene counts (> 10) in barcode by \n Cellranger and Optimus (normalized), \nadj. R square for random barcode pairs", format(summary(fit)$adj.r.squared, digits=4), sep=":"))
    dev.off()
 }


# BOX PLOTS for samples "ischaemia_1", "pbmc4k", "pbmc8k", "t_4k"
tab1 <- read.csv(paste("ischaemia_1", "normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
tab1<-tab1[ tab1$optimus > 10,][, 1]

tab2 <- read.csv(paste("pbmc4k", "normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
tab2<-tab2[ tab2$optimus > 10,][, 1]
 
tab3 <- read.csv(paste("pbmc8k", "normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
tab3<-tab3[ tab3$optimus > 10,][, 1]
 
tab4 <- read.csv(paste("t_4k", "normalized_50000_v2.txt", sep="_"), header=TRUE, sep='\t')
tab4<-tab4[ tab4$optimus > 10,][, 1]

png("boxplot_correlation.v1.png")
boxplot(tab1, tab2, tab3, tab4,  medcol = "magenta", outcex=0.05, names=c("ischeamia_1", "PMBC 4K", "PBMC 8K", "Pan T 4K"), main="correlations for barcodes with gene counts of more than 10")
dev.off()

