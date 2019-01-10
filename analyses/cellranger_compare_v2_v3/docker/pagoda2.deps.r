install.packages(c('dendsort', 'RcppArmadillo', 'rjson', 'RMTstat', 'Rook', 'fastcluster'),
                 repos='https://cloud.r-project.org')
                 
install.packages("devtools", repos='https://cloud.r-project.org')
                                        # Install Bioconductor dependencies
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager", repos='https://cloud.r-project.org')
BiocManager::install(ask=F,update=F)
BiocManager::install("GO.db", version = "3.8",ask=F,update=F)
BiocManager::install("org.Hs.eg.db", version = "3.8",ask=F,update=F)
BiocManager::install("org.Mm.eg.db", version = "3.8",ask=F,update=F)
BiocManager::install("pcaMethods", version = "3.8",ask=F,update=F)
library(devtools)
install_github("igraph/rigraph") # Don't install with install.packages()
install_github("jkrijthe/Rtsne",ref="openmp")
install.packages(c("Cairo","urltools"), repos='https://cloud.r-project.org')

