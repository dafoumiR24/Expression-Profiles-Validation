

GSE_clear_low_quality_norm <- function(GSE_acc,  cels_REMOVE_AFTER_QC,dir_download_GSES ){
  library(affy)
  library(simpleaffy)
    library(frma)
  
  dir.create(path = paste(dir_download_GSES,"/", GSE_acc, "/", "data_after_QC_filter", sep = ""), recursive = TRUE);
  
  setwd(paste(dir_download_GSES,"/", GSE_acc, "/", "data", sep = ""))
  phenodata <- read.delim( "phenodata.txt")
  if(length(cels_REMOVE_AFTER_QC)>0){
  pheno_data_new <- phenodata[-which(phenodata[,1]%in%cels_REMOVE_AFTER_QC==TRUE),]
  } else {
    pheno_data_new <- phenodata
  }
    
  file.copy(from=list.files(pattern = "CEL"), to= paste(dir_download_GSES,"/", GSE_acc, "/", "data_after_QC_filter", sep = ""),copy.mode = TRUE)
  setwd(paste(dir_download_GSES,"/", GSE_acc, "/", "data_after_QC_filter", sep = ""))
  write.table(pheno_data_new,"phenodata.txt",sep='\t',quote=FALSE,row.names=F,col.names=T) 
  if(length(cels_REMOVE_AFTER_QC)>0){
    
  file.remove(cels_REMOVE_AFTER_QC)
  }
  celfiles_after_QC<-read.affy(covdesc="phenodata.txt")
  celfiles_after_QC_frma <- frma(celfiles_after_QC)
  
  
  save(file= paste(GSE_acc,"_celfiles_after_QC_frma.RData", sep = ""), celfiles_after_QC_frma)
       
  setwd(dir_download_GSES)
  setwd('..')
  
}