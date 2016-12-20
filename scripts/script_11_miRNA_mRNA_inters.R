miRNA_mRNA_inters <- function (NORM_EXPR_VALUES_MRNA,NORM_EXPR_VALUES_MIRNA,gene_list, miRNA_list,method= "pearson",adjust_pvalue_method= NULL, miRNA_annot, cor_pvalue_cutoff,outputFileFolder) {
  
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  NORM_EXPR_VALUES_MRNA <- NORM_EXPR_VALUES_MRNA[row.names(NORM_EXPR_VALUES_MRNA)%in%gene_list,]
  NORM_EXPR_VALUES_MIRNA <- NORM_EXPR_VALUES_MIRNA[row.names(NORM_EXPR_VALUES_MIRNA)%in%miRNA_list,]
  
  NORM_EXPR_VALUES_MRNA <- as.matrix(NORM_EXPR_VALUES_MRNA[,order(colnames(NORM_EXPR_VALUES_MRNA))])
  NORM_EXPR_VALUES_MIRNA <- as.matrix(NORM_EXPR_VALUES_MIRNA[,order(colnames(NORM_EXPR_VALUES_MIRNA))])
  miRNA_annot[,1] <- as.character(miRNA_annot[,1])
  miRNA_annot[,2] <- as.character(miRNA_annot[,2])
  
  
  
correlation.matrix<-matrix(NA,nrow=nrow(NORM_EXPR_VALUES_MIRNA),ncol=nrow(NORM_EXPR_VALUES_MRNA))
colnames(correlation.matrix)<-rownames(NORM_EXPR_VALUES_MRNA)
rownames(correlation.matrix)<-rownames(NORM_EXPR_VALUES_MIRNA)
pval.matrix<-matrix(NA,nrow=nrow(NORM_EXPR_VALUES_MIRNA),ncol=nrow(NORM_EXPR_VALUES_MRNA))
colnames(pval.matrix)<-rownames(NORM_EXPR_VALUES_MRNA)
rownames(pval.matrix)<-rownames(NORM_EXPR_VALUES_MIRNA)

  
  for (i in 1:nrow(NORM_EXPR_VALUES_MIRNA)) {
    for (j in 1:nrow(NORM_EXPR_VALUES_MRNA)) {
      x<-NORM_EXPR_VALUES_MIRNA[i,]
      y<-NORM_EXPR_VALUES_MRNA[j,]
      cor.t<-cor.test(x,y,method=method, alternative = "less")
      correlation.matrix[i,j]<-cor.t$estimate
      pval.matrix[i,j]<-cor.t$p.value
    }
  }
  


library("reshape2")
cor_wide_table <- correlation.matrix
pval_wide_table <- pval.matrix

Correlation.melted_cor <- melt(cor_wide_table, value.name="CorrelationCoef")
colnames(Correlation.melted_cor) <- c("miRNA", "mRNA", "cor")

Correlation.melted_cor$miRNA_mRNA_inter <- paste(Correlation.melted_cor$miRNA, Correlation.melted_cor$mRNA, sep = ":")

Correlation.melted_pval <- melt(pval_wide_table, value.name="pval")
colnames(Correlation.melted_pval)[c(1,2)] <- c("miRNA", "mRNA")

Correlation.melted_pval$miRNA_mRNA_inter <- paste(Correlation.melted_pval$miRNA, Correlation.melted_pval$mRNA, sep = ":")

Correlation.melted_cor_pval <- merge(Correlation.melted_cor, Correlation.melted_pval[,c(3,4)], by ="miRNA_mRNA_inter" )


 MODIFIED_microCosm_v5_18 <- read.delim("MODIFIED_microCosm_v5_18.txt")
 MODIFIED_targetScan_v6.2_18 <- read.delim("MODIFIED_targetScan_v6.2_18.txt")
Correlation.melted_cor_pval$Target_site_info_MicroCosm <- rep(0, nrow(Correlation.melted_cor_pval))
Correlation.melted_cor_pval$Target_site_info_MicroCosm[which(Correlation.melted_cor_pval$miRNA_mRNA_inter%in%row.names(MODIFIED_microCosm_v5_18)==TRUE)]= 1

Correlation.melted_cor_pval$Target_site_info_Targetscan <- rep(0, nrow(Correlation.melted_cor_pval))
Correlation.melted_cor_pval$Target_site_info_Targetscan[which(Correlation.melted_cor_pval$miRNA_mRNA_inter%in%row.names(MODIFIED_targetScan_v6.2_18)==TRUE)]= 1
Correlation.melted_cor_pval$Target_site_info_sum <- Correlation.melted_cor_pval$Target_site_info_MicroCosm+Correlation.melted_cor_pval$Target_site_info_Targetscan

if (length(adjust_pvalue_method)>0){
  
  Correlation.melted_cor_pval$adj_pval <- p.adjust(Correlation.melted_cor_pval$pval,method=adjust_pvalue_method)
  
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- subset(Correlation.melted_cor_pval, Correlation.melted_cor_pval$cor < 0 &Correlation.melted_cor_pval$Target_site_info_sum>=1&Correlation.melted_cor_pval$adj_pval<=cor_pvalue_cutoff)
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- merge(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS, miRNA_annot, by.x = "miRNA", by.y = colnames(miRNA_annot)[1])
  
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <-  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[moveme(names(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS), paste(colnames(miRNA_annot)[2] ,"first",sep = " "))]
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[order(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$adj_pval),]
  
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[,3] <- NULL
  colnames(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS)[c(1,2,3)] <- c("MIR_SYMBOLS", "MIMAT_SYMBOLS", "MRNAS")
  
}  else {
    
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- subset(Correlation.melted_cor_pval, Correlation.melted_cor_pval$cor < 0 &Correlation.melted_cor_pval$Target_site_info_sum>=1&Correlation.melted_cor_pval$pval<=cor_pvalue_cutoff)
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- merge(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS, miRNA_annot, by.x = "miRNA", by.y = colnames(miRNA_annot)[1])
  
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <-  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[moveme(names(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS), paste(colnames(miRNA_annot)[2] ,"first",sep = " "))]
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS <- CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[order(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$pval),]
 
  CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[,3] <- NULL
  colnames(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS)[c(1,2,3)] <- c("MIR_SYMBOLS", "MIMAT_SYMBOLS", "MRNAS")
  

}
  
  
  
ALL_MIRNAS <- row.names(cor_wide_table)
ALL_MIRNAS_STAT_SEM <- unique(as.vector(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MIMAT_SYMBOLS))


mIrna_inters <- NULL
Freq <- NULL
for (i in 1:length(ALL_MIRNAS_STAT_SEM)){
  
  
  mIrna_inters <-c(mIrna_inters, paste(as.vector(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[which(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MIMAT_SYMBOLS==ALL_MIRNAS_STAT_SEM[i]),3]), collapse = " ,"))
  Freq <- c(Freq, length(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[which(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MIMAT_SYMBOLS==ALL_MIRNAS_STAT_SEM[i]),3]))
  
}

mIrna_targets <- data.frame(MIRNAS= ALL_MIRNAS_STAT_SEM, freq= Freq, names= mIrna_inters)


MIRNAS_ZERO_FREQ <- ALL_MIRNAS[which(ALL_MIRNAS%in%ALL_MIRNAS_STAT_SEM==FALSE)]

if(length(MIRNAS_ZERO_FREQ)>0) {
  mIrna_targets_zero_freq <- data.frame(MIRNAS= MIRNAS_ZERO_FREQ, freq= rep(0, length(MIRNAS_ZERO_FREQ)), names= rep(" ", length(MIRNAS_ZERO_FREQ)))
  mIrna_targets <- rbind(mIrna_targets, mIrna_targets_zero_freq)
}

mIrna_targets <- mIrna_targets[order(mIrna_targets$freq, decreasing = TRUE),]
row.names(mIrna_targets) <- mIrna_targets$MIRNAS
mIrna_targets$MIRNAS <- NULL
  
  
  
ALL_MRNAS <- colnames(cor_wide_table)
ALL_MRNAS_STAT_SEM <- unique(as.vector(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MRNAS))


mRna_inters <- NULL
Freq <- NULL
for (i in 1:length(ALL_MRNAS_STAT_SEM)){
  
  
  mRna_inters <-c(mRna_inters, paste(as.vector(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[which(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MRNAS==ALL_MRNAS_STAT_SEM[i]),2]), collapse = " ,"))
  Freq <- c(Freq, length(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS[which(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MRNAS==ALL_MRNAS_STAT_SEM[i]),2]))
  
}

mRna_targets <- data.frame(MRNAS= ALL_MRNAS_STAT_SEM, freq= Freq, names= mRna_inters)


MRNAS_ZERO_FREQ <- ALL_MRNAS[which(ALL_MRNAS%in%ALL_MRNAS_STAT_SEM==FALSE)]

if(length(MRNAS_ZERO_FREQ)>0) {
  mRna_targets_zero_freq <- data.frame(MRNAS= MRNAS_ZERO_FREQ, freq= rep(0, length(MRNAS_ZERO_FREQ)), names= rep(" ", length(MRNAS_ZERO_FREQ)))
  mRna_targets <- rbind(mRna_targets, mRna_targets_zero_freq)
}

mRna_targets <- mRna_targets[order(mRna_targets$freq, decreasing = TRUE),]
row.names(mRna_targets) <- mRna_targets$MRNAS
mRna_targets$MRNAS <- NULL





miRNA_Targets <- mIrna_targets
mir_accession <- as.vector(row.names(miRNA_Targets))
miRNA_Targets <- cbind(mir_accession, miRNA_Targets)
miRNA_Targets <- merge(miRNA_Targets, miRNA_annot, by.x = "mir_accession", by.y = colnames(miRNA_annot)[1])
miRNA_Targets <-  miRNA_Targets[moveme(names(miRNA_Targets), paste(colnames(miRNA_annot)[2],"first", sep= " "))]
miRNA_Targets <- miRNA_Targets[order(- miRNA_Targets$freq),]
colnames(miRNA_Targets)[c(1,2)] <- c("MIR_SYMBOLS","MIMAT_SYMBOLS")
colnames(miRNA_Targets)[4] <- "MRNAS"
row.names(miRNA_Targets) <- NULL

mRNA_Targets <- mRna_targets

GENES <- row.names(mRNA_Targets)
mRNA_Targets_new <- cbind(GENES, mRNA_Targets)
s <- strsplit(as.character(mRNA_Targets_new$names), split = " ,")
mRNA_Targets_new <- data.frame(GENES = rep(mRNA_Targets_new$GENES, sapply(s, length)),freq = rep(mRNA_Targets_new$freq, sapply(s, length)), names = unlist(s))
#returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
mRNA_Targets_new$names <- trim(mRNA_Targets_new$names)
mRNA_Targets_new <- mRNA_Targets_new[mRNA_Targets_new$freq>=1,]
MIR_SYMBOLS <- NULL
for (i in 1:dim(mRNA_Targets_new)[1]){
  MIMAT <- as.character(mRNA_Targets_new[i, which(colnames(mRNA_Targets_new)=="names")])
  MIR_SYMBOLS <- c(MIR_SYMBOLS,as.character(miRNA_annot[grep(MIMAT,as.character(miRNA_annot[,1]) ),2]))
  
}
mRNA_Targets_new <- cbind(mRNA_Targets_new,MIR_SYMBOLS)
IND_ZERO_FREQ <- which(mRNA_Targets$freq==0)  
mRNA_targets_zero_freq <- mRNA_Targets[IND_ZERO_FREQ,]
remove(GENES)
GENES <- row.names(mRNA_targets_zero_freq)
mRNA_targets_zero_freq <- cbind(GENES, mRNA_targets_zero_freq)
row.names(mRNA_targets_zero_freq) <- NULL
mRNA_targets_zero_freq["MIR_SYMBOLS"] <- " "
FINAL_mRNA_targets <- rbind(mRNA_Targets_new,mRNA_targets_zero_freq)

GENES_SYMBOLS <- c()
for (i in 1:length(levels(FINAL_mRNA_targets$GENES))){
  GENES_SYMBOLS <- c(GENES_SYMBOLS, as.character(levels(FINAL_mRNA_targets$GENES)[i]))
}
FREQUENCY <- c()
for (i in 1:length(levels(FINAL_mRNA_targets$GENES))){
  FREQUENCY <- c(FREQUENCY, as.character(FINAL_mRNA_targets[which(FINAL_mRNA_targets$GENES==levels(FINAL_mRNA_targets$GENES)[i])[1],which(colnames(FINAL_mRNA_targets)=="freq")]))
}

MIMAT_NAMES <- c()
for (i in 1:length(levels(FINAL_mRNA_targets$GENES))){
  MIMAT_NAMES <- c(MIMAT_NAMES, paste(as.character(FINAL_mRNA_targets[c(which(FINAL_mRNA_targets$GENES==levels(FINAL_mRNA_targets$GENES)[i])),which(colnames(FINAL_mRNA_targets)=="names")]),collapse = ", "))
}

MIRNA_SYMBOLS <- c()
for (i in 1:length(levels(FINAL_mRNA_targets$GENES))){
  MIRNA_SYMBOLS <- c(MIRNA_SYMBOLS, paste(as.character(FINAL_mRNA_targets[c(which(FINAL_mRNA_targets$GENES==levels(FINAL_mRNA_targets$GENES)[i])),which(colnames(FINAL_mRNA_targets)=="MIR_SYMBOLS")]),collapse = ", "))
}

mRNA_TARGETS_WITH_MIMAT <- data.frame(MRNAS =GENES_SYMBOLS, freq= FREQUENCY, MIMAT_SYMBOLS= MIMAT_NAMES, MIR_SYMBOLS= MIRNA_SYMBOLS)
mRNA_TARGETS_WITH_MIMAT$freq <- as.numeric(as.character(mRNA_TARGETS_WITH_MIMAT$freq))
mRNA_TARGETS_WITH_MIMAT <- mRNA_TARGETS_WITH_MIMAT[order( mRNA_TARGETS_WITH_MIMAT$freq, decreasing = TRUE),]
row.names(mRNA_TARGETS_WITH_MIMAT) <- NULL


NODE_ATTRIBUTES_MRNA <- data.frame(FEATURES = as.vector(unique(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MRNAS)), FEATURE_TYPE= rep("MRNA", length(as.vector(unique(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MRNAS)))) ) 
NODE_ATTRIBUTES_MIRNA <- data.frame(FEATURES = as.vector(unique(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MIR_SYMBOLS)), FEATURE_TYPE= rep("MIRNA", length(as.vector(unique(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$MIR_SYMBOLS)))) ) 
NODE_ATTRIBUTES <- rbind(NODE_ATTRIBUTES_MRNA, NODE_ATTRIBUTES_MIRNA)

CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS$COR_TYPE <- rep("NEG", nrow(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS))

write.table(CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS, file = paste(outputFileFolder,"/","CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(miRNA_Targets, file = paste(outputFileFolder,"/","miRNA_Targets.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(mRNA_TARGETS_WITH_MIMAT, file = paste(outputFileFolder,"/","mRNA_TARGETS_WITH_MIMAT.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(NODE_ATTRIBUTES, file = paste(outputFileFolder,"/","NODE_ATTRIBUTES_MIR_MRNA_COR_NET.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



#assign("CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS",CYTOSCAPE_MATRIX_FILTER_NEGATIVE_CORS,envir = .GlobalEnv)
#assign("miRNA_Targets",miRNA_Targets,envir = .GlobalEnv)
#assign("mRNA_TARGETS_WITH_MIMAT",mRNA_TARGETS_WITH_MIMAT,envir = .GlobalEnv)
#assign("NODE_ATTRIBUTES_MIR_MRNA_COR_NET",NODE_ATTRIBUTES,envir = .GlobalEnv)


}


