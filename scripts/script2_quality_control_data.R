

GSE_quality_control <- function(GSE_acc, dir_download_GSES ){
  
library(affy)
library(simpleaffy)
library(frma)
library(RColorBrewer)
library(affyPLM)

setwd(paste(dir_download_GSES,"/", GSE_acc, "/", "data", sep = ""))
celfiles<-read.affy(covdesc= "phenodata.txt")
sample_num = length(list.files( pattern = "CEL"))

celfiles_frma <- frma(celfiles)

# set colour palette
cols <- brewer.pal(8, "Set1")

dir.create(path = paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control", sep = ""), recursive = TRUE);


png(paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "Boxplots_before_norm.png",sep = "")  ,width = 1200, height = 1000); 
boxplot(celfiles, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM required to interrogate celfiles.gcrma
png( paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "Boxplots_after_norm.png",sep = "") ,width = 1200, height = 1000); 
boxplot(celfiles_frma, col=cols)
dev.off()

# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data
#png(paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "Histograms_before_norm.png",sep = ""),width = 1200, height = 1000); 
#hist(celfiles, col=cols)
#dev.off()

# Plot a density vs log intensity histogram for the normalised data
#png(paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "Histograms_after_norm.png",sep = "") ,width = 1200, height = 1000); 
#hist(celfiles_frma, col=cols)
#dev.off()


celfiles.qc <- fitPLM(celfiles)


for(i in 1:sample_num){
 # print(paste("PLM", i, sep = "") )
  png( paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "PLM",i, ".png",sep=""),width = 1200, height = 1000); 
  image(celfiles.qc, which=i, add.legend=TRUE)
  dev.off()
}

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.  GSM524665.CEL is an outlier
RLE_stats <- RLE(celfiles.qc,type="stats")
RLE_stats <- as.data.frame(RLE_stats)
write.table(RLE_stats, paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "RLE_stats.txt",sep = "") ,sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 

assign(paste(GSE_acc,"RLE_stats", sep = ""), RLE_stats, envir = .GlobalEnv)

png( paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "RLE.png",sep = "") ,width = 1200, height = 1000); 
RLE(celfiles.qc, main="RLE")
dev.off()

# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
NUSE_STATS<- NUSE(celfiles.qc, type="stats")
write.table(NUSE_STATS, paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "NUSE_STATS.txt",sep = ""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 

assign(paste(GSE_acc,"NUSE_STATS", sep = ""), NUSE_STATS, envir = .GlobalEnv)

png(  paste(dir_download_GSES,"/", GSE_acc, "/", "Quality_Control","/", "NUSE.png",sep = ""),width = 1200, height = 1000); 
NUSE(celfiles.qc, main="NUSE")
dev.off()
setwd(dir_download_GSES)
setwd('..')

}