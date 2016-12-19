

DEA_MICROARRAYS <- function(Exprs_table, Pheno_data_table,sample_type_colname, diff_exp_type,adjust_pvalue_method="BH", LogFC_cuttoff=NULL, Adj_p_value_cuttoff=NULL,  annotation_microarray_probesets=FALSE,  basename= "hgu133plus2", outputFileFolder ){
 
library(limma)
dir.create(path = outputFileFolder, recursive = TRUE)

diff_exp_type_split <- as.character(strsplit(as.character(diff_exp_type),split="_vs_" )[[1]])


Pheno_data_table_new <- Pheno_data_table[Pheno_data_table[, sample_type_colname]%in%diff_exp_type_split,]
Exprs_table_new <- Exprs_table[,colnames(Exprs_table)%in%row.names(Pheno_data_table)]
Exprs_table_new <- Exprs_table_new[,order(colnames(Exprs_table_new))]
Pheno_data_table_new <- Pheno_data_table_new[order(row.names(Pheno_data_table_new)),]


Sample_type <- as.vector(Pheno_data_table_new[,sample_type_colname])

design <- model.matrix(~0+Sample_type)

fit <- lmFit(Exprs_table_new,design)

contrast.matrix <- makeContrasts(contrasts=paste("Sample_type", diff_exp_type_split[1], "-","Sample_type", diff_exp_type_split[2], sep = "" ), levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
TOPTABLE <- toptable(fit2,num="all", adjust=adjust_pvalue_method)

TOPTABLE <- TOPTABLE[order(TOPTABLE$adj.P.Val),]
if(!is.null(LogFC_cuttoff)){
  TOPTABLE <- TOPTABLE[TOPTABLE$logFC>=LogFC_cuttoff, ]
}

if(!is.null(Adj_p_value_cuttoff)){
  TOPTABLE <- TOPTABLE[TOPTABLE$adj.P.Val<=Adj_p_value_cuttoff, ]
}

TOPTABLE <- TOPTABLE[order(TOPTABLE$adj.P.Val),]

if(annotation_microarray_probesets){ 
  library("hgu95av2.db")
  library("hgu133plus2.db")
  
 library(annotate)
TOPTABLE$gene.symbols <- getSYMBOL(row.names(TOPTABLE), basename)

}

assign(paste("TOPTABLE_",diff_exp_type, sep="") ,TOPTABLE, envir = .GlobalEnv )

write.table(TOPTABLE,paste(outputFileFolder, "/",  "TOPTABLE_",diff_exp_type,".txt", sep=""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)


png(paste(outputFileFolder, "/", "volcanoplot_", diff_exp_type,".png",sep=""),width = 900, height = 600); 
volcanoplot(fit2,coef=1)
dev.off()


}