
  
Prepr_DEA_RNA_MIRNA_SEQ <- function(Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type,data_type, adjust_pvalue_method="BH", LogFC_cuttoff=NULL, Adj_p_value_cuttoff=NULL, outputFileFolder){
  
  dir.create(path = outputFileFolder, recursive = TRUE)
  
library(edgeR)
library(limma)


diff_exp_type_split <- as.character(strsplit(as.character(diff_exp_type),split="_vs_" )[[1]])


Pheno_data_table_new <- Pheno_data_table[Pheno_data_table[, sample_type_colname]%in%diff_exp_type_split,]
Exprs_table_new <- Exprs_table[,colnames(Exprs_table)%in%row.names(Pheno_data_table)]

Exprs_table_new <- Exprs_table_new[,order(colnames(Exprs_table_new))]
Pheno_data_table_new <- Pheno_data_table_new[order(row.names(Pheno_data_table_new)),]

paste(cat("features before filtering"," = ", nrow(Exprs_table_new), "\n", sep = ""))


Exprs_table_new$EXPR_PERC <- apply(Exprs_table_new, 1, function(x) length(which(x>=10))/ncol(Exprs_table_new)*100)

Exprs_table_new <- as.data.frame(Exprs_table_new[Exprs_table_new$EXPR_PERC>=90,])
Exprs_table_new$EXPR_PERC <- NULL

paste(cat("features after filtering"," = ", nrow(Exprs_table_new), "\n", sep = ""))

y <- DGEList(counts=as.matrix(Exprs_table_new))


##Apply normalization
y <- calcNormFactors(y, method="TMM")

Sample_type <- as.vector(Pheno_data_table_new[,sample_type_colname])

design <- model.matrix(~0+Sample_type)

rownames(design) <- colnames(y);

v <- voom(y,design) 

fit <- lmFit(v,design)

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


assign(paste("TOPTABLE_",diff_exp_type, "_", data_type, sep="") ,TOPTABLE, envir = .GlobalEnv )

write.table(TOPTABLE,paste(outputFileFolder,"/", "TOPTABLE_",diff_exp_type, "_", data_type,".txt", sep=""),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)


png(paste(outputFileFolder, "/", "volcanoplot_", diff_exp_type, "_", data_type,".png",sep=""),width = 900, height = 600); 
volcanoplot(fit2,coef=1)
dev.off()

 

}