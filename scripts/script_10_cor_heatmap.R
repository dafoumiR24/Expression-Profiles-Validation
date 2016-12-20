  

cor_heatmap <- function(NORM_EXPR_VALUES, features_of_interest=NULL, outputFileFolder){
  
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
  if(length(features_of_interest)>0){
    for(i in 1: length(features_of_interest)){ 
      
  
  library(corrplot)
  source("http://www.sthda.com/upload/rquery_cormat.r")
  
  
  cormat <- rquery.cormat(data.frame(t(NORM_EXPR_VALUES)),graph=FALSE, type = "FULL")
  
  cor_wide_table <- cormat$r
  
  pval_wide_table <- cormat$p
  cor_wide_table_NEG <- cor_wide_table[unique(c(as.vector(which(cor_wide_table[features_of_interest[i],]<=0)), which(row.names(cor_wide_table)==features_of_interest[i]))), unique(c(as.vector(which(cor_wide_table[,features_of_interest[i]]<=0)), which(row.names(cor_wide_table)==features_of_interest[i])))]
  cor_wide_table_POS <- cor_wide_table[unique(c(as.vector(which(cor_wide_table[features_of_interest[i],]>0)), which(row.names(cor_wide_table)==features_of_interest[i]))), unique(c(as.vector(which(cor_wide_table[,features_of_interest[i]]>0)), which(row.names(cor_wide_table)==features_of_interest[i])))]
  
  pval_wide_table_NEG <- pval_wide_table[row.names(pval_wide_table)%in%row.names(cor_wide_table_NEG), colnames(pval_wide_table)%in%colnames(cor_wide_table_NEG)]
  pval_wide_table_POS <- pval_wide_table[row.names(pval_wide_table)%in%row.names(cor_wide_table_POS), colnames(pval_wide_table)%in%colnames(cor_wide_table_POS)]
  
  library("reshape2")
  
  
  LOWER.TRIANGLE_cor_NEG <- lower.tri(cor_wide_table_NEG,diag=F)
  CORRELATIONSLOWER.TRIANGLE_cor_NEG <- cor_wide_table_NEG
  CORRELATIONSLOWER.TRIANGLE_cor_NEG [!LOWER.TRIANGLE_cor_NEG] <- NA
  Correlation.melted_cor_NEG <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_cor_NEG, value.name="CorrelationCoef"))
  colnames(Correlation.melted_cor_NEG) <- c("Feature1", "Feature2", "cor")
  
  Correlation.melted_cor_NEG$Feature1_Feature2_inter <- paste(Correlation.melted_cor_NEG$Feature1, Correlation.melted_cor_NEG$Feature2, sep = ":")
  
  LOWER.TRIANGLE_pval_NEG <- lower.tri(pval_wide_table_NEG,diag=F)
  CORRELATIONSLOWER.TRIANGLE_pval_NEG <- pval_wide_table_NEG
  CORRELATIONSLOWER.TRIANGLE_pval_NEG [!LOWER.TRIANGLE_pval_NEG] <- NA
  Correlation.melted_pval_NEG <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_pval_NEG, value.name="pval"))
  colnames(Correlation.melted_pval_NEG)[c(1,2)] <- c("Feature1", "Feature2")
  
  Correlation.melted_pval_NEG$Feature1_Feature2_inter <- paste(Correlation.melted_pval_NEG$Feature1, Correlation.melted_pval_NEG$Feature2, sep = ":")
  
  Correlation.melted_cor_pval_NEG <- merge(Correlation.melted_cor_NEG, Correlation.melted_pval_NEG[,c(3,4)], by ="Feature1_Feature2_inter" )
  
  Correlation.melted_cor_pval_NEG$Feature1_Feature2_inter <- NULL
  Correlation.melted_cor_pval_NEG$Feature1 <- as.character(Correlation.melted_cor_pval_NEG$Feature1)
  Correlation.melted_cor_pval_NEG$Feature2 <- as.character(Correlation.melted_cor_pval_NEG$Feature2)
  
  

  
  
  
  
  
  LOWER.TRIANGLE_cor_POS <- lower.tri(cor_wide_table_POS,diag=F)
  CORRELATIONSLOWER.TRIANGLE_cor_POS <- cor_wide_table_POS
  CORRELATIONSLOWER.TRIANGLE_cor_POS [!LOWER.TRIANGLE_cor_POS] <- NA
  Correlation.melted_cor_POS <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_cor_POS, value.name="CorrelationCoef"))
  colnames(Correlation.melted_cor_POS) <- c("Feature1", "Feature2", "cor")
  
  Correlation.melted_cor_POS$Feature1_Feature2_inter <- paste(Correlation.melted_cor_POS$Feature1, Correlation.melted_cor_POS$Feature2, sep = ":")
  
  LOWER.TRIANGLE_pval_POS <- lower.tri(pval_wide_table_POS,diag=F)
  CORRELATIONSLOWER.TRIANGLE_pval_POS <- pval_wide_table_POS
  CORRELATIONSLOWER.TRIANGLE_pval_POS [!LOWER.TRIANGLE_pval_POS] <- NA
  Correlation.melted_pval_POS <- na.omit(melt(CORRELATIONSLOWER.TRIANGLE_pval_POS, value.name="pval"))
  colnames(Correlation.melted_pval_POS)[c(1,2)] <- c("Feature1", "Feature2")
  
  Correlation.melted_pval_POS$Feature1_Feature2_inter <- paste(Correlation.melted_pval_POS$Feature1, Correlation.melted_pval_POS$Feature2, sep = ":")
  
  Correlation.melted_cor_pval_POS <- merge(Correlation.melted_cor_POS, Correlation.melted_pval_POS[,c(3,4)], by ="Feature1_Feature2_inter" )
  
  Correlation.melted_cor_pval_POS$Feature1_Feature2_inter <- NULL
  Correlation.melted_cor_pval_POS$Feature1 <- as.character(Correlation.melted_cor_pval_POS$Feature1)
  Correlation.melted_cor_pval_POS$Feature2 <- as.character(Correlation.melted_cor_pval_POS$Feature2)
  
  
  
    
  Correlation.melted_cor_NEG_4_PLOT <- melt(cor_wide_table_NEG, value.name="CorrelationCoef")
  
library(ggplot2)


png(paste(outputFileFolder, "/",features_of_interest[i], "_NEG_COR_HEATMAP_",".png",sep=""),width = 900, height = 600); 

print(ggplot(data =Correlation.melted_cor_NEG_4_PLOT, aes(Var2, Var1, fill = CorrelationCoef))+
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "green", mid = "black",
                       midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed())
dev.off()


Correlation.melted_cor_POS_4_PLOT <- melt(cor_wide_table_POS, value.name="CorrelationCoef")


png(paste(outputFileFolder, "/",features_of_interest[i], "_POS_COR_HEATMAP_",".png",sep=""),width = 900, height = 600); 

print(ggplot(data =Correlation.melted_cor_POS_4_PLOT, aes(Var2, Var1, fill = CorrelationCoef))+
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "green", mid = "black",
                       midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed())
dev.off()


write.table(Correlation.melted_cor_pval_POS, file = paste(outputFileFolder,"/",features_of_interest[i] , "_POS.CORS_PVALS.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#assign(paste(features_of_interest[i] , "_POS.CORS_PVALS", sep = ""),Correlation.melted_cor_pval_POS,envir = .GlobalEnv)


write.table(Correlation.melted_cor_pval_NEG, file = paste(outputFileFolder,"/",features_of_interest[i] , "_NEG.CORS_PVALS.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#assign(paste(features_of_interest[i] , "_NEG.CORS_PVALS", sep = ""),Correlation.melted_cor_pval_NEG,envir = .GlobalEnv)

}

} else {
  
  
  
  library(corrplot)
  source("http://www.sthda.com/upload/rquery_cormat.r")
  
  
  cormat <- rquery.cormat(data.frame(t(NORM_EXPR_VALUES)),graph=FALSE, type = "FULL")
  
  cor_wide_table <- cormat$r
  
  pval_wide_table <- cormat$p
  
  library("reshape2")
  
  
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
  
  
  
  Correlation.melted_cor_4_PLOT <- melt(cor_wide_table, value.name="CorrelationCoef")
  
  
  png(paste(outputFileFolder, "/", "COR_HEATMAP_",".png",sep=""),width = 900, height = 600); 
  
  print(ggplot(data =Correlation.melted_cor_4_PLOT, aes(Var2, Var1, fill = CorrelationCoef))+
    geom_tile()+
    scale_fill_gradient2(low = "red", high = "green", mid = "black",
                         midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed())
  dev.off()
  
  
  write.table(Correlation.melted_cor_pval, file = paste(outputFileFolder,"/", "CORS_PVALS.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  #assign( "CORS_PVALS",Correlation.melted_cor_pval,envir = .GlobalEnv)
  
}
}
