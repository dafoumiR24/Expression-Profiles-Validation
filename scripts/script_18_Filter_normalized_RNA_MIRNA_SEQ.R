

Filter_normalized_RNA_MIRNA_SEQ <- function(Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type,data_type, outputFileFolder ){
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
  library(edgeR)
  library(limma)
  
  
  diff_exp_type_split <- as.character(strsplit(as.character(diff_exp_type),split="_vs_" )[[1]])
  
  
  Pheno_data_table_new <- Pheno_data_table[Pheno_data_table[, sample_type_colname]%in%diff_exp_type_split,]
  Exprs_table_new <- Exprs_table[,colnames(Exprs_table)%in%row.names(Pheno_data_table)]
  
  Exprs_table_new <- Exprs_table_new[,order(colnames(Exprs_table_new))]
  Pheno_data_table_new <- Pheno_data_table_new[order(row.names(Pheno_data_table_new)),]
  
  paste(cat("features before filtering"," = ", nrow(Exprs_table_new), "\n", sep = ""))
  
  
  Exprs_table_new$EXPR_PERC <- apply(Exprs_table_new, 1, function(x) length(which(x>=10))/ncol(Exprs_table_new)*100)
  
  Exprs_table_new <- as.data.frame(Exprs_table_new[Exprs_table_new$EXPR_PERC>=90,])
  Exprs_table_new$EXPR_PERC <- NULL
  
  paste(cat("features after filtering"," = ", nrow(Exprs_table_new), "\n", sep = ""))
  
  y <- DGEList(counts=as.matrix(Exprs_table_new))
  
  
  ##Apply normalization
  y <- calcNormFactors(y, method="TMM")
  
  Sample_type <- as.vector(Pheno_data_table_new[,sample_type_colname])
  
  design <- model.matrix(~0+Sample_type)
  
  rownames(design) <- colnames(y);
  
  v <- voom(y,design) 
  
 
  NORM_EXPR_VALUES <-  as.data.frame(v$E)
  PHENO_DATA <- data.frame(SAMPLES=as.vector(colnames(NORM_EXPR_VALUES)), Sample_type=Sample_type)
  row.names(PHENO_DATA) <- PHENO_DATA$SAMPLES
  colnames(PHENO_DATA)[2] <- sample_type_colname
  assign(paste("NORM_EXPR_VALUES_",diff_exp_type,"_", data_type, sep="") ,NORM_EXPR_VALUES, envir = .GlobalEnv )
  
  write.table(NORM_EXPR_VALUES,paste(outputFileFolder,"/", "NORM_EXPR_VALUES_",diff_exp_type, "_", data_type,".txt", sep=""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  
  assign(paste("PHENO_DATA_",diff_exp_type,"_", data_type, sep="") ,PHENO_DATA, envir = .GlobalEnv )
  
  write.table(PHENO_DATA,paste(outputFileFolder,"/", "PHENO_DATA_",diff_exp_type, "_", data_type,".txt", sep=""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
  
    
  
}