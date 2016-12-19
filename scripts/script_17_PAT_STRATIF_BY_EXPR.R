
PAT_STRATIF_BY_EXPR <- function(Exprs_table, Pheno_data_table, gene_mirnas, outputFileFolder){

  dir.create(path = outputFileFolder, recursive = TRUE)
  
  Exprs_table <- as.data.frame(t(Exprs_table))
  INTER <- intersect(row.names(Pheno_data_table), row.names(Exprs_table))
  
  Exprs_table <- Exprs_table[row.names(Exprs_table)%in%INTER,]
  Pheno_data_table <- Pheno_data_table[row.names(Pheno_data_table)%in%INTER,]
  
  Exprs_table <- Exprs_table[order(row.names(Exprs_table)),]
  
  Pheno_data_table <- Pheno_data_table[order(row.names(Pheno_data_table)),]
  
  for(i in 1:length(gene_mirnas)){
    Pheno_data_table$PATIENT_EXPR_GROUP <- rep(NA,dim(Pheno_data_table)[1])
    
  PATIENTS_LOW_EXPR <- row.names(Exprs_table) [which(Exprs_table[,gene_mirnas[i]]<= quantile(as.numeric(Exprs_table[,gene_mirnas[i]]),0.25))]
  
  PATIENTS_HIGH_EXPR <- row.names(Exprs_table) [which(Exprs_table[,gene_mirnas[i]]>= quantile(as.numeric(Exprs_table[,gene_mirnas[i]]),0.75))]
  PATIENTS_MEDIUM_EXPR <- row.names(Exprs_table)[which(Exprs_table[,gene_mirnas[i]]> quantile(as.numeric(Exprs_table[,gene_mirnas[i]]),0.25)&Exprs_table[,gene_mirnas[i]]< quantile(as.numeric(Exprs_table[,gene_mirnas[i]]),0.75))]
  
  Pheno_data_table$PATIENT_EXPR_GROUP[which(row.names(Pheno_data_table)%in%PATIENTS_LOW_EXPR==TRUE)] <- rep("LOW_EXPR", length(PATIENTS_LOW_EXPR))
  Pheno_data_table$PATIENT_EXPR_GROUP[which(row.names(Pheno_data_table)%in%PATIENTS_HIGH_EXPR==TRUE)] <- rep("HIGH_EXPR", length(PATIENTS_HIGH_EXPR))
  Pheno_data_table$PATIENT_EXPR_GROUP[which(row.names(Pheno_data_table)%in%PATIENTS_MEDIUM_EXPR==TRUE)] <- rep("MEDIUM_EXPR", length(PATIENTS_MEDIUM_EXPR))
  colnames(Pheno_data_table)[which(colnames(Pheno_data_table)=="PATIENT_EXPR_GROUP")] <- paste("PATIENTS_GROUPS_BY_", gene_mirnas[i], sep = "")
  
  }
  
  write.table(Pheno_data_table, file = paste(outputFileFolder,"/","PHENO_DATA_PAT_STRATIF_BY_EXPR.txt", sep=""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  assign("PHENO_DATA_PAT_STRATIF_BY_EXPR",Pheno_data_table,envir = .GlobalEnv)
  
}
  
  
  