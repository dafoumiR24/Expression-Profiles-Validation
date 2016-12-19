
##################################### R Scripts for RNA-seq/miRNA-seq/RNA microarray data analysis ##############################

A. REQUIREMENTS
------------------------------------------------------------------------------------------------------------------------------------------
 1. Windows 10
 2. R version 3.2.1 (2015-06-18)
 3. R packages: affy, simpleaffy, frma, RColorBrewer, affyPLM, inSilicoMerging, amap, genefilter, limma, hgu133plus2.db, annotate, gplots, devtools, easyGgplot2, edgeR, corrplot, reshape2, ggplot2, ggbiplot

B. INSTALLATION AND SET UP
-------------------------------------------------------------------------------------------------------------------------------------------

1.Into R session, Install R/BioConductor dependencies 
  download the BioC installation routines
 source("http://bioconductor.org/biocLite.R")
  install the core packages. It will take some time!!
  biocLite(c())
  install the Bioconductor packages
 biocLite(c("affy", "simpleaffy", "frma", "affyPLM", "inSilicoMerging", "genefilter", "limma", "hgu133plus2.db", "annotate", "edgeR"))
  install the CRAN packages
 install.packages(c("RColorBrewer", "amap", "gplots", "devtools", "corrplot", "reshape2", "ggplot2" ))
  install the Github packages
 library(devtools)
install_github(c("kassambara/easyGgplot2", "vqv/ggbiplot"))

2.Download 





