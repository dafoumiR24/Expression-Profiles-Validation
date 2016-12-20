cluster_heatmap <- function(Exprs_table, Pheno_data_table, vars2study, varType= c("cat"),dist_method,clust_method,main_title, outputFileFolder){
library("gplots")
library("devtools")
#Load latest version of heatmap.3 function
#Set a working directory for output files
dir.create(path = outputFileFolder, recursive = TRUE)

Exprs_table <- Exprs_table[, order(colnames(Exprs_table))]

Pheno_data_table <- Pheno_data_table[order(row.names(Pheno_data_table)), ]

 library(RColorBrewer)
 #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
 #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_vector = c("red","green3","blue","cyan","magenta","yellow","gray","black",
           "orange","darkred","green","darkblue","darkcyan","darkmagenta",
           "darkorchid1","darkgoldenrod3","aquamarine","antiquewhite",
           "darkolivegreen3")
colfunc <- colorRampPalette(c("green","yellow", "red"))


cat_var_elements <- NULL
cat_var_colors <- NULL
cat_var_elements_4_legend <- NULL
cat_var_colors_4_legend <- NULL
for(i in 1:length(vars2study)){
  if(varType[i]=="cat"){
    cat_var_elements <- c(cat_var_elements, unique(as.vector(Pheno_data_table[,vars2study[i]])))
    cat_var_colors <- c(cat_var_colors, col_vector[length(cat_var_colors)+1:length(cat_var_colors)+length(unique(as.vector(Pheno_data_table[,vars2study[i]])))])
    cat_var_elements_4_legend <- c(cat_var_elements_4_legend, c("", unique(as.vector(Pheno_data_table[,vars2study[i]]))))
    cat_var_colors_4_legend <- c(cat_var_colors_4_legend, c("white",col_vector[length(cat_var_colors_4_legend)+1:length(cat_var_colors_4_legend)+length(unique(as.vector(Pheno_data_table[,vars2study[i]])))]))
    
    
  }
  
  if(varType[i]=="ord"){
    cat_var_elements <- c(cat_var_elements, sort(unique(as.vector(Pheno_data_table[,vars2study[i]]))))
    cat_var_colors <- c(cat_var_colors, colfunc(length(unique(as.vector(Pheno_data_table[,vars2study[i]])))))
    cat_var_elements_4_legend <- c(cat_var_elements_4_legend, c("", sort(unique(as.vector(Pheno_data_table[,vars2study[i]])))))
    cat_var_colors_4_legend <- c(cat_var_colors_4_legend, c("white", colfunc(length(unique(as.vector(Pheno_data_table[,vars2study[i]]))))))
    
  }
  
}
cat_var_elements_colors <- cbind(cat_var_elements, cat_var_colors)
cat_var_elements_4_legend <- cat_var_elements_4_legend[-1]
cat_var_colors_4_legend <- cat_var_colors_4_legend[-1]
VAR_COLOR_DATAFRAME <- data.frame(ID = 1:nrow(Pheno_data_table), PATIENTS= row.names(Pheno_data_table))
    
for(i in 1:length(vars2study)){
  V <- rep(NA, nrow(VAR_COLOR_DATAFRAME))
  
  V_pheno_elements <- unique(as.vector(Pheno_data_table[, vars2study[i]]))
  
  for(k in 1:length(V_pheno_elements)){
  V [as.vector(Pheno_data_table[, vars2study[i]])==V_pheno_elements[k]] <- as.vector(cat_var_elements_colors[which(cat_var_elements_colors[,1]==V_pheno_elements[k]),2])[1]
  
  }
  VAR_COLOR_DATAFRAME <- cbind(VAR_COLOR_DATAFRAME, V)
  colnames(VAR_COLOR_DATAFRAME)[2+i] <- vars2study[i]
}
 
clab= as.matrix(VAR_COLOR_DATAFRAME[, 3:ncol(VAR_COLOR_DATAFRAME)])



#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method=dist_method)}
myclust=function(c) {hclust(c,method=clust_method)}

#Create heatmap using custom heatmap.3 source code loaded above
png(paste(outputFileFolder, "/", "Heatmap_", main_title,".png",sep=""),width = 1200, height = 1000); 
#pdf(paste("Heatmap_", main_title,".pdf",sep="")) 

heatmap.3(Exprs_table, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row", dendrogram="both", Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE,density.info="none", trace="none", main=main_title, labCol=FALSE,  col=redgreen)
legend("topright",legend=cat_var_elements_4_legend ,fill=cat_var_colors_4_legend, border=FALSE, bty="n")
dev.off()




}
