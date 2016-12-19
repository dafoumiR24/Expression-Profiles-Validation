

GSE_download <- function(GSE_acc, outputFileFolder ){
library(GEOquery)
dir.create(path = outputFileFolder, recursive = TRUE)
setwd(outputFileFolder)

getGEOSuppFiles(GSE_acc)

#Unpack the CEL files
setwd( GSE_acc)
untar(paste(GSE_acc,"_RAW.tar", sep = ""), exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")
assign(paste(GSE_acc,"_cels", sep = ""), cels, envir = .GlobalEnv)

setwd("..")
setwd("..")

}

