
PCA <- function(Exprs_table, Pheno_data_table, pheno_var2study,PCs= c(1,2), labels = TRUE,circle= FALSE, ellipse= FALSE,var.axes= FALSE, outputFileFolder){ 
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
  Exprs_table <- as.data.frame(t(Exprs_table))
  Exprs_table <- Exprs_table[order(row.names(Exprs_table)),]
  
  Pheno_data_table <- Pheno_data_table[order(row.names(Pheno_data_table)),]
  
    
pca <- prcomp( Exprs_table , scale = TRUE, center = TRUE)

png(paste(outputFileFolder,"/", "Variances_vs_PCs_Explaining_",  pheno_var2study ,".png",sep=""),width = 900, height = 600)
A <- plot(pca, type = "l")
print(A)
dev.off()



PCA_SUMMARY<- summary(pca)

save (PCA_SUMMARY, file= paste("PCA_SUMMARY_Explaining_", pheno_var2study, ".RData", sep = "")) 


library(ggbiplot)

png(paste(outputFileFolder,"/", "PCA_Explaining_", pheno_var2study,".png",sep=""),width = 900, height = 600)


if(labels){
  
g <- ggbiplot(pca,choices = PCs, obs.scale = 1, var.scale = 1, 
              groups = as.vector(Pheno_data_table[,pheno_var2study]), ellipse = ellipse, 
              circle = circle, labels = as.vector(row.names(Exprs_table)), var.axes= var.axes)
} else {
  
  g <- ggbiplot(pca,choices = PCs, obs.scale = 1, var.scale = 1, 
                groups = as.vector(Pheno_data_table[,pheno_var2study]), ellipse = ellipse, 
                circle = circle, labels = NULL, var.axes= var.axes )
}

g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

dev.off()


}



