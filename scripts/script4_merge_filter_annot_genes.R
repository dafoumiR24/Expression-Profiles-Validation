GSE_merge_filter_annot_genes <- function(GSE_acc_4_merge, targetAnnot, batchAnnot,  transform_method= "COMBAT", dir_download_GSES ){
  
  
  COLORS = c("red","green3","blue","cyan","magenta","yellow","gray","black",
             "orange","darkred","green","darkblue","darkcyan","darkmagenta",
             "darkorchid1","darkgoldenrod3","aquamarine","antiquewhite",
             "darkolivegreen3");
  
  # create color map
  #-------------------------------------------------------------------------------
  makeColorMap = function(eset, label)
  {
    colMap = list();
    vec = unique(as.vector(pData(eset)[,label]));
    for(i in 1:length(vec)) { colMap[[ vec[i] ]] = COLORS[i]; }
    return(colMap);
  }  
  
  # create color vector for plots
  #-------------------------------------------------------------------------------
  makeColorVec = function(eset, label, colMap, o=TRUE)
  {
    if(o)
    {
      labels = as.vector(pData(eset)[colnames(exprs(eset)),label]);
    }
    else
    {
      labels = as.vector(pData(eset)[,label]);
    }
    return(as.vector(unlist(sapply(labels, function(x) { colMap[x]; }))));
  }
  
  # plot colored MDS plot
  #-------------------------------------------------------------------------------
  plotMDS = function(eset, batchAnnot="Batch", targetAnnot=NULL, legend=TRUE, 
                     file=NULL, ...)
  {
    if(is.null(targetAnnot)) { targetAnnot = batchAnnot; }
    
    if(! is.null(file)) { pdf(file, width=12, height=7); }
    
    mds = cmdscale(dist(t(exprs(eset))), eig=TRUE);
    
    colMap = makeColorMap(eset, targetAnnot); 
    colVec = makeColorVec(eset, targetAnnot, colMap);
    
    #--Add margin to the right for the legend 
    tmp = par()$mar;
    if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
    
    range_x = range(mds$points[,1]);
    range_y = range(mds$points[,2]);
    
    plot(mds$points,
         col=colVec,
         pch=as.numeric(pData(eset)[,batchAnnot]),
         panel.first={ U = par("usr");
                       rect(U[1],U[3],U[2],U[4],
                            col="azure2",
                            border="black",
                            lwd=3)},
         lwd=2,
         xlab="",
         ylab="",
         xlim=range_x,
         ylim=range_y,
         ...);
    
    if(legend)
    {
      x = range_x[2] + (range_x[2]-range_x[1])*0.1;
      y = range_y[2] - (range_y[2]-range_y[1])*0.1;
      
      syms = unique(pData(eset)[,batchAnnot])
      legend(x,y,
             legend = syms,
             pt.lwd=2,
             pch = as.numeric(syms),
             box.lwd=3,
             bg="azure2");
      
      legend(x,y-(length(syms)*(range_y[2]-range_y[1])*0.1),
             legend = names(colMap),
             pt.lwd=2,
             pch=19,
             col = unlist(colMap),
             box.lwd=3,
             bg="azure2");
      
      #--Reset margin
      par(xpd=F, mar=tmp)
    }
    
    if(! is.null(file)) { dev.off(); }
  }
  
  
  # plot dendrogram
  #-------------------------------------------------------------------------------
  plotDendrogram = function(eset, batchAnnot="Batch", legend=TRUE, file=NULL, ...)
  {
    if(! is.null(file)) { pdf(file, width=12, height=7); }
    
    #--Add margin to the right for the legend 
    tmp = par()$mar;
    if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
    
    #--Create clustering
    hc = hcluster(t(exprs(eset)), method="pearson");
    hc$labels = as.numeric(pData(eset)[,batchAnnot]);
    
    #--Plot tree
    plot(hc,
         panel.first={U = par("usr");
                      rect(U[1],U[3],U[2],U[4],
                           col="azure2",
                           border="black",
                           lwd=3);},
         xlab="", ylab="", sub="", axes=FALSE,
         ...);
    
    #--Reset margin
    if(legend) { par(xpd=F, mar=tmp) }
    
    if(! is.null(file)) { dev.off(); }
  }
  
  
  dir.create(paste(dir_download_GSES,"/", "merged_datasets", sep = ""), recursive = TRUE);
  
esets <- list()
for(i in 1:length(GSE_acc_4_merge)){
  setwd(paste(dir_download_GSES,"/", GSE_acc_4_merge[i], "/",  "data_after_QC_filter", sep = ""))
  load(paste(GSE_acc_4_merge[i], "_celfiles_after_QC_frma.RData", sep = ""))
  esets[[length(esets)+1]]<-celfiles_after_QC_frma
    remove(celfiles_after_QC_frma)
}
  
setwd(paste(dir_download_GSES,"/", "merged_datasets", sep = ""))

library(inSilicoMerging);
eset_NONE = merge(esets, method="NONE")
png( "PLOTMDS_NO_TRANS.png",width = 1200, height = 1000); 
plotMDS(eset_NONE, targetAnnot= targetAnnot, batchAnnot= batchAnnot , main="NONE (No Transformation)");
dev.off()


eset_TRANSF = merge(esets, method=transform_method);
png(paste( "PLOTMDS_", transform_method, ".png"),width = 1200, height = 1000); 
plotMDS(eset_TRANSF, targetAnnot= targetAnnot,batchAnnot= batchAnnot , main=transform_method);
dev.off()


# par(mfrow=c(1,2))
library(amap)
png( "plotDendrogram_NO_TRANS.png",width = 1200, height = 1000)
plotDendrogram(eset_NONE, batchAnnot= batchAnnot, legend=FALSE, main="NONE");
dev.off()

png(paste( "plotDendrogram_", transform_method, ".png"),width = 1200, height = 1000)
plotDendrogram(eset_TRANSF, batchAnnot= batchAnnot,legend=FALSE,main=transform_method)
dev.off()

save(file= "eset_NONE.RData", eset_NONE)
save(file= "eset_TRANSF.RData", eset_TRANSF)
assign("eset_NONE", eset_NONE, envir = .GlobalEnv)
assign("eset_TRANSF", eset_TRANSF, envir = .GlobalEnv)


library(genefilter)
eset_TRANSF_filtered <- nsFilter(eset_TRANSF, require.entrez=FALSE, remove.dupEntrez=FALSE)
eset_TRANSF_filtered <-eset_TRANSF_filtered$eset
TRANS_FILTER_Exprs_table <- as.data.frame(exprs(eset_TRANSF_filtered))
TRANS_FILTER_Pheno_data_table <- pData(eset_TRANSF_filtered)
assign("TRANS_FILTER_Exprs_table", TRANS_FILTER_Exprs_table, envir = .GlobalEnv)
assign("TRANS_FILTER_Pheno_data_table", TRANS_FILTER_Pheno_data_table, envir = .GlobalEnv)
write.table(TRANS_FILTER_Exprs_table,"TRANS_FILTER_Exprs_table.txt",sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 
write.table(TRANS_FILTER_Pheno_data_table,"TRANS_FILTER_Pheno_data_table.txt",sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 


eset_UNTRANSF_filtered <- nsFilter(eset_NONE, require.entrez=FALSE, remove.dupEntrez=FALSE)
eset_UNTRANSF_filtered <-eset_UNTRANSF_filtered$eset
UNTRANS_FILTER_Exprs_table <- as.data.frame(exprs(eset_UNTRANSF_filtered))
UNTRANS_FILTER_Pheno_data_table <- pData(eset_UNTRANSF_filtered)
assign("UNTRANS_FILTER_Exprs_table", UNTRANS_FILTER_Exprs_table, envir = .GlobalEnv)
assign("UNTRANS_FILTER_Pheno_data_table", UNTRANS_FILTER_Pheno_data_table, envir = .GlobalEnv)
write.table(UNTRANS_FILTER_Exprs_table,"UNTRANS_FILTER_Exprs_table.txt",sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 
write.table(UNTRANS_FILTER_Pheno_data_table,"UNTRANS_FILTER_Pheno_data_table.txt",sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE) 


setwd(dir_download_GSES)
setwd("..")


}

