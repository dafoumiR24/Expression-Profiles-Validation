
Feature_cor_net <- function(NORM_EXPR_VALUES, feature_list,feature_type,cor_pval_cutoff, outputFileFolder){
  
 dir.create(path = outputFileFolder, recursive = TRUE)
 NORM_EXPR_VALUES <- NORM_EXPR_VALUES[row.names(NORM_EXPR_VALUES)%in%feature_list,]
 
 
 library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r")


cormat <- rquery.cormat(data.frame(t(NORM_EXPR_VALUES)),graph=FALSE, type = "FULL")

cor_wide_table <- cormat$r

pval_wide_table <- cormat$p

library("reshape2")
#library("Matrix")


LOWER.TRIANGLE_cor <- lower.tri(cor_wide_table,diag=F)
CORRELATIONSLOWER.TRIANGLE_cor <- cor_wide_table
CORRELATIONSLOWER.TRIANGLE_cor [!LOWER.TRIANGLE_cor] <- NA
Correlation.melted_cor <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_cor, value.name="CorrelationCoef"))
colnames(Correlation.melted_cor) <- c("Feature1", "Feature2", "cor")
    
Correlation.melted_cor$Feature1_Feature2_inter <- paste(Correlation.melted_cor$Feature1, Correlation.melted_cor$Feature2, sep = ":")

LOWER.TRIANGLE_pval <- lower.tri(pval_wide_table,diag=F)
CORRELATIONSLOWER.TRIANGLE_pval <- pval_wide_table
CORRELATIONSLOWER.TRIANGLE_pval [!LOWER.TRIANGLE_pval] <- NA
Correlation.melted_pval <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_pval, value.name="pval"))
colnames(Correlation.melted_pval)[c(1,2)] <- c("Feature1", "Feature2")

Correlation.melted_pval$Feature1_Feature2_inter <- paste(Correlation.melted_pval$Feature1, Correlation.melted_pval$Feature2, sep = ":")

Correlation.melted_cor_pval <- merge(Correlation.melted_cor, Correlation.melted_pval[,c(3,4)], by ="Feature1_Feature2_inter" )

Correlation.melted_cor_pval$Feature1_Feature2_inter <- NULL
Correlation.melted_cor_pval$Feature1 <- as.character(Correlation.melted_cor_pval$Feature1)
Correlation.melted_cor_pval$Feature2 <- as.character(Correlation.melted_cor_pval$Feature2)
ALL_FEATURES <- unique(as.vector(c(Correlation.melted_cor_pval$Feature1, Correlation.melted_cor_pval$Feature2)))


Correlation.melted_cor_pval_stat_sem <- Correlation.melted_cor_pval[Correlation.melted_cor_pval$pval<=cor_pval_cutoff,]
ALL_FEATURES_STAT_SEM <- unique(as.vector(c(Correlation.melted_cor_pval_stat_sem$Feature1, Correlation.melted_cor_pval_stat_sem$Feature2)))


  feature_inters <- NULL
Freq <- NULL
for (i in 1:length(ALL_FEATURES_STAT_SEM)){
  
  
  feature_inters <-c(feature_inters, paste(as.vector(c(Correlation.melted_cor_pval_stat_sem[which(Correlation.melted_cor_pval_stat_sem$Feature1==ALL_FEATURES_STAT_SEM[i]),2], Correlation.melted_cor_pval_stat_sem[which(Correlation.melted_cor_pval_stat_sem$Feature2==ALL_FEATURES_STAT_SEM[i]),1] )), collapse = " ,"))
  Freq <- c(Freq, length(c(Correlation.melted_cor_pval_stat_sem[which(Correlation.melted_cor_pval_stat_sem$Feature1==ALL_FEATURES_STAT_SEM[i]),2], Correlation.melted_cor_pval_stat_sem[which(Correlation.melted_cor_pval_stat_sem$Feature2==ALL_FEATURES_STAT_SEM[i]),1] )))

}

feature_targets <- data.frame(FEATURES= ALL_FEATURES_STAT_SEM, freq= Freq, FEATURE_INTERS= feature_inters)


FEATURES_ZERO_FREQ <- ALL_FEATURES[which(ALL_FEATURES%in%ALL_FEATURES_STAT_SEM==FALSE)]

if(length(FEATURES_ZERO_FREQ)>0) {
  feature_targets_zero_freq <- data.frame(FEATURES= FEATURES_ZERO_FREQ, freq= rep(0, length(FEATURES_ZERO_FREQ)), FEATURE_INTERS= rep(" ", length(FEATURES_ZERO_FREQ)))
  feature_targets <- rbind(feature_targets, feature_targets_zero_freq)
}

feature_targets <- feature_targets[order(feature_targets$freq, decreasing = TRUE),]
  row.names(feature_targets) <- NULL
  

node_attribute <- data.frame(FEATURES= ALL_FEATURES_STAT_SEM, FEATURE_TYPE= rep(feature_type, length(ALL_FEATURES_STAT_SEM)))

Correlation.melted_cor_pval_stat_sem$COR_TYPE <- rep(NA, nrow(Correlation.melted_cor_pval_stat_sem))

Correlation.melted_cor_pval_stat_sem$COR_TYPE [which(Correlation.melted_cor_pval_stat_sem$cor<=0)]= rep("NEG", length(which(Correlation.melted_cor_pval_stat_sem$cor<=0)))
Correlation.melted_cor_pval_stat_sem$COR_TYPE [which(Correlation.melted_cor_pval_stat_sem$cor>0)]= rep("POS", length(which(Correlation.melted_cor_pval_stat_sem$cor>0)))


write.table(node_attribute, file = paste(outputFileFolder,"/","NODE_ATTRIBUTES_", feature_type  ,".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(feature_targets, file = paste(outputFileFolder,"/",  feature_type ,"_targets.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(Correlation.melted_cor_pval_stat_sem, file = paste(outputFileFolder,"/",  "CYTOSCAPE_",feature_type, "_COR_NET" ,".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



assign(paste("NODE_ATTRIBUTES_", feature_type, sep = ""),node_attribute,envir = .GlobalEnv)
assign(paste(feature_type, "_targets", sep = ""),feature_targets,envir = .GlobalEnv)
assign(paste("CYTOSCAPE_",feature_type, "_COR_NET", sep = ""),Correlation.melted_cor_pval_stat_sem,envir = .GlobalEnv)



}
