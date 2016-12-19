
KM_PLOT_GENERATION <- function(Exprs_table, time_event_data, genes_or_mirnas,xlabel= "time (months)", ylabel= "survival probality", outputFileFolder ){
  
dir.create(path = outputFileFolder, recursive = TRUE)

inter <- intersect(row.names(time_event_data), colnames(Exprs_table))
time_event_data <- time_event_data[row.names(time_event_data)%in%inter,]
Exprs_table <- Exprs_table[, colnames(Exprs_table)%in%inter]
time_event_data <- time_event_data[order(row.names(time_event_data)),]
Exprs_table <- Exprs_table[,order(colnames(Exprs_table))]




for (i in 1:length(genes_or_mirnas)){
  
  PATIENTS_High_expression <-  colnames(Exprs_table)[which(Exprs_table[genes_or_mirnas[i],]>=median(as.numeric(Exprs_table[genes_or_mirnas[i] ,])))]
  
  PATIENTS_LOW_expression <-  colnames(Exprs_table)[which(Exprs_table[genes_or_mirnas[i],]<median(as.numeric(Exprs_table[genes_or_mirnas[i] ,])))]
  x1= length(PATIENTS_LOW_expression)
  x2= length(PATIENTS_High_expression)
  
  time_event_data$feature_expr_group <- rep(NA, dim(time_event_data)[1])
  time_event_data$feature_expr_group [which(row.names(time_event_data)%in%PATIENTS_LOW_expression==TRUE)] <- rep("A",length(PATIENTS_LOW_expression) )
  
  time_event_data$feature_expr_group [which(row.names(time_event_data)%in%PATIENTS_High_expression==TRUE)] <- rep("B",length(PATIENTS_High_expression)) 
  
  library(survival)
  
    
  surv <- Surv(time_event_data$time,time_event_data$event)
  kmfit <-  survfit(surv ~ time_event_data$feature_expr_group, conf.type = "log-log")
  
   
  model2 <- coxph(surv ~ time_event_data$feature_expr_group, method="breslow")
  sumcph <- summary(model2)
  pval <- sumcph$sctest[[3]]
  HAZ_RAT <- sumcph$coefficients[[2]]
 
  pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
  HAZ_RAT_txt <- paste("HR =", signif(HAZ_RAT, 3))
  
  png(paste(outputFileFolder,"/", genes_or_mirnas[i],"_", "kmplot.png",sep=""),width = 900, height = 600)
  plot1 <-   plot(kmfit,mark.time= TRUE, lty = 1:2, col= 2:3, xlab = xlabel, ylab = ylabel)
  
  legend1 <-  legend("bottomleft", c(paste("LOW EXPRESSION = ",x1,sep=""), paste("HIGH EXPRESSION = ",x2,sep="")), lty = 1:2,col= 2:3)
  title1 <- title(paste("KM-PLOT ", genes_or_mirnas[i],sep="" ))
  legend2 <- legend("bottomright",paste(pvaltxt, "\n",HAZ_RAT_txt, sep = ""), plot=T)
  
  dev.off() 
  
  
  
  time_event_data$feature_expr_group <- NULL
}
}