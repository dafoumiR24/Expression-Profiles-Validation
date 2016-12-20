 Boxplots <- function(Exprs_table, Pheno_data_table, var2study,group_names, genes_or_mirnas,addDot=TRUE,  x_label, y_label,  outputFileFolder){
  
   dir.create(path = outputFileFolder, recursive = TRUE)
   
  groups = c()
  for (i in 1: length(group_names)){
    
    
    groups = c(groups, replicate (length(which(Pheno_data_table[,var2study]==group_names[i])), group_names[i]))
  }
  
  for (i in 1:length(genes_or_mirnas)){
    values <- NULL
    
    
    for (j in 1: length(group_names)){
      
      values <- c( values,as.numeric(Exprs_table[which(row.names(Exprs_table)==genes_or_mirnas[i]),colnames(Exprs_table)%in%row.names(Pheno_data_table[Pheno_data_table[,var2study]==group_names[j],])]))
    }
    
    df <- data.frame(group =groups, value = values, stringsAsFactors = TRUE)
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    library(easyGgplot2)
    png( paste(outputFileFolder,"/",  "Boxplot_",genes_or_mirnas[i], ".png",sep=""),width = 900, height = 600); 
    
    
    print(ggplot2.boxplot(data=df, xName='group',yName='value', groupName='group',
                  showLegend=FALSE,  groupColors= sample(col_vector, length(group_names)),
                  backgroundColor="white", xtitle=x_label, ytitle=y_label,  
                  mainTitle=genes_or_mirnas[i],
                  addDot=addDot, dotSize=1.7,
                  removePanelGrid=TRUE,removePanelBorder=FALSE,
                  axisLine=c(0.5, "solid", "black"),
                  dotPosition= "jitter", jitter=0.2))
    dev.off()      
  }
 }
  