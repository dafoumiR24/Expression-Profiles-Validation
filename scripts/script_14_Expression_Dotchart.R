
Expression_Dotchart <- function(Exprs_table, Pheno_data_table, var2study,group_names, genes_or_mirnas,  x_label, y_label, outputFileFolder){
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
  Exprs_table <- as.data.frame(t(Exprs_table))
  
  Pheno_data_table <- Pheno_data_table[Pheno_data_table[,var2study]%in%group_names,]
  
  Exprs_table <- Exprs_table[row.names(Exprs_table)%in%row.names(Pheno_data_table),]
  Exprs_table$SAMPLES <- as.vector(row.names(Exprs_table))
  
  Pheno_data_table$SAMPLES <- as.vector(row.names(Pheno_data_table))
  
  for(i in 1:length(genes_or_mirnas)){
  
    
  x <- merge(Exprs_table[,colnames(Exprs_table)%in%c("SAMPLES",genes_or_mirnas[i] )], Pheno_data_table[,colnames(Pheno_data_table)%in%c("SAMPLES",var2study )], by = "SAMPLES" )
  
  x[,var2study] <- as.factor(x[,var2study])
    
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
 
  x$colors <- rep(NA, dim(x)[1])
  for(j in 1: length(group_names)){
   
  x$colors[x[,var2study]==group_names[j]]= col_vector[j]

 }

x <- x[order(x[,genes_or_mirnas[i]]),]


colors_4_legend <- NULL
for(k in 1:length(group_names)){
  colors_4_legend <- c(colors_4_legend, as.character(x[which(x[,var2study]==group_names[k])[1], which(colnames(x)=="colors")]))

}

png(paste(outputFileFolder, "/", "Dotchart_of_", genes_or_mirnas[i],".png",sep=""),width = 900, height = 600); 

dotchart(x[,genes_or_mirnas[i]],labels=row.names(x),groups= x[,var2study],main=paste("Dotchart_of_", genes_or_mirnas[i], sep = ""),xlab=x_label,ylab= y_label, color=x$colors, lcolor = "white")
legend("bottomright",legend=group_names  ,fill=colors_4_legend, border=FALSE, bty="n")
dev.off()


}

}
