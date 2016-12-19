
PREPROCESS_AND_COMMON_MIRNA_MRNA <- function(RNA_seq_data, miRNA_seq_data,TCGA_source=TRUE, outputFileFolder){
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
 if(TCGA_source){
  Normal_samples_RNA_seq  <- as.vector(colnames(RNA_seq_data)[which(substr(colnames(RNA_seq_data),14,14) == "1")])
  Tumor_samples_RNA_seq <- as.vector(colnames(RNA_seq_data)[which(substr(colnames(RNA_seq_data),14,14) == "0")])
  
  Normal_samples_miRNA_seq  <- as.vector(colnames(miRNA_seq_data)[which(substr(colnames(miRNA_seq_data),14,14) == "1")])
  Tumor_samples_miRNA_seq <- as.vector(colnames(miRNA_seq_data)[which(substr(colnames(miRNA_seq_data),14,14) == "0")])
  
  RNA_seq_data_normal <-as.data.frame(t( RNA_seq_data[,colnames(RNA_seq_data)%in%Normal_samples_RNA_seq]))
  RNA_seq_data_tumor <- as.data.frame(t(RNA_seq_data[,colnames(RNA_seq_data)%in%Tumor_samples_RNA_seq]))
  
  miRNA_seq_data_normal <- as.data.frame(t(miRNA_seq_data[, colnames(miRNA_seq_data)%in%Normal_samples_miRNA_seq]))
  miRNA_seq_data_tumor <-as.data.frame(t( miRNA_seq_data[, colnames(miRNA_seq_data)%in%Tumor_samples_miRNA_seq]))
  
  RNA_seq_data_normal$BRC_SAMPLES <- as.vector(row.names(RNA_seq_data_normal))
  RNA_seq_data_tumor$BRC_SAMPLES <- as.vector(row.names(RNA_seq_data_tumor))

  miRNA_seq_data_normal$BRC_SAMPLES <- as.vector(row.names(miRNA_seq_data_normal))
  miRNA_seq_data_tumor$BRC_SAMPLES <- as.vector(row.names(miRNA_seq_data_tumor))
  
  
  
  
  RNA_seq_data_normal$BRC_PATIENTS <- substr(RNA_seq_data_normal$BRC_SAMPLES,1,12)
  RNA_seq_data_tumor$BRC_PATIENTS <- substr(RNA_seq_data_tumor$BRC_SAMPLES,1,12)
  
  miRNA_seq_data_normal$BRC_PATIENTS <- substr(miRNA_seq_data_normal$BRC_SAMPLES,1,12)
  miRNA_seq_data_tumor$BRC_PATIENTS <- substr(miRNA_seq_data_tumor$BRC_SAMPLES,1,12)
  
  RNA_seq_data_normal <- RNA_seq_data_normal[!duplicated(RNA_seq_data_normal$BRC_PATIENTS),]
  RNA_seq_data_tumor <- RNA_seq_data_tumor[!duplicated(RNA_seq_data_tumor$BRC_PATIENTS), ]
  
  miRNA_seq_data_normal <- miRNA_seq_data_normal[!duplicated(miRNA_seq_data_normal$BRC_PATIENTS),]
  miRNA_seq_data_tumor <- miRNA_seq_data_tumor[!duplicated(miRNA_seq_data_tumor$BRC_PATIENTS),]
  
  INTER_NORMAL <- intersect(RNA_seq_data_normal$BRC_PATIENTS, miRNA_seq_data_normal$BRC_PATIENTS)
  INTER_TUMOR <- intersect(RNA_seq_data_tumor$BRC_PATIENTS, miRNA_seq_data_tumor$BRC_PATIENTS)
  
  RNA_seq_data_normal <- RNA_seq_data_normal[RNA_seq_data_normal$BRC_PATIENTS%in%INTER_NORMAL,]
  miRNA_seq_data_normal <- miRNA_seq_data_normal[miRNA_seq_data_normal$BRC_PATIENTS%in%INTER_NORMAL,]
  RNA_seq_data_tumor <- RNA_seq_data_tumor[RNA_seq_data_tumor$BRC_PATIENTS%in%INTER_TUMOR,]
  miRNA_seq_data_tumor <- miRNA_seq_data_tumor[miRNA_seq_data_tumor$BRC_PATIENTS%in%INTER_TUMOR,]
  
  RNA_seq_data_normal <- RNA_seq_data_normal[order(RNA_seq_data_normal$BRC_PATIENTS),]
  miRNA_seq_data_normal <- miRNA_seq_data_normal[order(miRNA_seq_data_normal$BRC_PATIENTS),]
  RNA_seq_data_tumor <- RNA_seq_data_tumor[order(RNA_seq_data_tumor$BRC_PATIENTS),]
  miRNA_seq_data_tumor <- miRNA_seq_data_tumor[order(miRNA_seq_data_tumor$BRC_PATIENTS),]
  
  
  RNA_seq_data_normal$PATIENT_NUM_CODE <- paste("N.", 1:nrow(RNA_seq_data_normal), sep = "")
  miRNA_seq_data_normal$PATIENT_NUM_CODE <- paste("N.", 1:nrow(miRNA_seq_data_normal), sep = "")

  RNA_seq_data_tumor$PATIENT_NUM_CODE <- paste("T.", 1:nrow(RNA_seq_data_tumor), sep = "")
  miRNA_seq_data_tumor$PATIENT_NUM_CODE <- paste("T.", 1:nrow(miRNA_seq_data_tumor), sep = "")
  
  RNA_seq_data_normal$SAMPLE_TYPE <- rep("NORMAL", nrow(RNA_seq_data_normal))
  miRNA_seq_data_normal$SAMPLE_TYPE <- rep("NORMAL", nrow(miRNA_seq_data_normal))
  
  RNA_seq_data_tumor$SAMPLE_TYPE <- rep("TUMOR", nrow(RNA_seq_data_tumor))
  miRNA_seq_data_tumor$SAMPLE_TYPE <- rep("TUMOR", nrow(miRNA_seq_data_tumor))
  
  
  RNA_seq_data_new <- rbind(RNA_seq_data_tumor, RNA_seq_data_normal)
  miRNA_seq_data_new <- rbind(miRNA_seq_data_tumor, miRNA_seq_data_normal)
  

  cat("IDENTICAL SAMPLES BETWEEN NEW MRNA AND MIRNA DATAFRAMES = ", (identical(RNA_seq_data_new$BRC_PATIENTS, miRNA_seq_data_new$BRC_PATIENTS)), "\n")
  PHENO_DATA_MIRNA <- miRNA_seq_data_new[, colnames(miRNA_seq_data_new)%in%c("PATIENT_NUM_CODE", "BRC_PATIENTS","BRC_SAMPLES", "SAMPLE_TYPE" )]
  PHENO_DATA_MRNA <- RNA_seq_data_new[, colnames(RNA_seq_data_new)%in%c("PATIENT_NUM_CODE", "BRC_PATIENTS","BRC_SAMPLES", "SAMPLE_TYPE" )]
  row.names(PHENO_DATA_MRNA) <- PHENO_DATA_MRNA$PATIENT_NUM_CODE
  PHENO_DATA_MRNA$PATIENT_NUM_CODE <- NULL
  
  row.names(PHENO_DATA_MIRNA) <- PHENO_DATA_MIRNA$PATIENT_NUM_CODE
  PHENO_DATA_MIRNA$PATIENT_NUM_CODE <- NULL
  
  EXPR_DATA_MIRNA_NEW <- miRNA_seq_data_new
  row.names(EXPR_DATA_MIRNA_NEW) <- EXPR_DATA_MIRNA_NEW$PATIENT_NUM_CODE
  
  EXPR_DATA_MIRNA_NEW <- as.data.frame(t(EXPR_DATA_MIRNA_NEW[,-which(colnames(EXPR_DATA_MIRNA_NEW)%in%c("PATIENT_NUM_CODE", "BRC_PATIENTS","BRC_SAMPLES", "SAMPLE_TYPE" )==TRUE)]))
  

  EXPR_DATA_RNA_SEQ_NEW <- RNA_seq_data_new
  row.names(EXPR_DATA_RNA_SEQ_NEW) <- EXPR_DATA_RNA_SEQ_NEW$PATIENT_NUM_CODE
  
  EXPR_DATA_RNA_SEQ_NEW <- as.data.frame(t(EXPR_DATA_RNA_SEQ_NEW[,-which(colnames(EXPR_DATA_RNA_SEQ_NEW)%in%c("PATIENT_NUM_CODE", "BRC_PATIENTS","BRC_SAMPLES", "SAMPLE_TYPE" )==TRUE)]))
  assign("EXPR_DATA_RNA_SEQ_NEW" ,EXPR_DATA_RNA_SEQ_NEW, envir = .GlobalEnv )
  
  write.table(EXPR_DATA_RNA_SEQ_NEW,paste(outputFileFolder, "/",  "EXPR_DATA_RNA_SEQ_NEW.txt", sep = "") ,sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  
  assign("EXPR_DATA_MIRNA_NEW" ,EXPR_DATA_MIRNA_NEW, envir = .GlobalEnv )
  
  write.table(EXPR_DATA_MIRNA_NEW,paste(outputFileFolder, "/", "EXPR_DATA_MIRNA_NEW.txt", sep = ""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  assign("PHENO_DATA_MRNA" ,PHENO_DATA_MRNA, envir = .GlobalEnv )
  
  write.table(PHENO_DATA_MRNA, paste(outputFileFolder, "/",  "PHENO_DATA_MRNA.txt", sep = ""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  assign("PHENO_DATA_MIRNA" ,PHENO_DATA_MIRNA, envir = .GlobalEnv )
  
  write.table(PHENO_DATA_MIRNA,paste(outputFileFolder, "/", "PHENO_DATA_MIRNA.txt", sep = "") ,sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  
  
  
} else {
  
  INTER_SAMPLES <- intersect(colnames(RNA_seq_data), colnames(miRNA_seq_data))
  RNA_seq_data_new <- RNA_seq_data[, colnames(RNA_seq_data)%in%INTER_SAMPLES]
  miRNA_seq_data_new <- miRNA_seq_data[, colnames(miRNA_seq_data)%in%INTER_SAMPLES]
  RNA_seq_data_new <- RNA_seq_data_new[, order(colnames(RNA_seq_data_new))]
  miRNA_seq_data_new <- miRNA_seq_data_new[, order(colnames(miRNA_seq_data_new))]
  cat("IDENTICAL SAMPLES BETWEEN NEW MRNA AND MIRNA DATAFRAMES = ", identical(colnames(RNA_seq_data_new),colnames(miRNA_seq_data_new)), "\n")
  
  
  assign("EXPR_DATA_RNA_SEQ_NEW" ,RNA_seq_data_new, envir = .GlobalEnv )
  
  write.table(RNA_seq_data_new,paste(outputFileFolder, "/", "EXPR_DATA_RNA_SEQ_NEW.txt", sep = ""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  assign("EXPR_DATA_MIRNA_SEQ_NEW" ,miRNA_seq_data_new, envir = .GlobalEnv )
  
  write.table(miRNA_seq_data_new,paste(outputFileFolder, "/", "EXPR_DATA_MIRNA_SEQ_NEW.txt", sep = ""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  
  
}

}

  