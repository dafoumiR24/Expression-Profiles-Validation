
##################################### R Scripts for RNA-seq/miRNA-seq/RNA microarray data analysis ##############################

#A. REQUIREMENTS
........................................................................................................................................

1. Windows 10
2. R version 3.2.1 (2015-06-18)
3. R packages: affy, simpleaffy, frma, RColorBrewer, affyPLM, inSilicoMerging, amap, genefilter, limma, hgu133plus2.db, annotate, gplots, devtools, easyGgplot2, edgeR, corrplot, reshape2, ggplot2, ggbiplot

#B. INSTALLATION AND SET UP
....................................................................................................................................... 1.Start R and install R/BioConductor dependencies  

download the BioC installation routines 
 source("http://bioconductor.org/biocLite.R")
install the core packages. It will take some time!!
 biocLite(c())
install the Bioconductor packages
 biocLite(c("affy", "simpleaffy", "frma", "affyPLM", "inSilicoMerging", "genefilter", "limma", "hgu133plus2.db", "annotate", "edgeR"))
install the CRAN packages
 install.packages(c("RColorBrewer", "amap", "gplots", "devtools", "corrplot", "reshape2", "ggplot2" ))
install the Github packages
 library(devtools)
 install_github(c("kassambara/easyGgplot2", "vqv/ggbiplot"))

2. Download files from 'scripts' directory to a preference directory. To set this directory as your working directory of R, type:
   WD <- "C:/Users/aaaa/bbbb/cccc/dddd"
   setwd(WD)

3. a. For GPL570 platform -Affymetrix Human Genome U133 Plus 2.0 microarray transcriptomic data- analyses you need to acquire the raw data directly
from GEO. Specifically for each GSE accession number (GSEXXXX) you need to download the .CEL files and put them into a new directory named "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/data". Also you need to capture the experimental information. This is just a text file which includes a data frame containing three columns. The first column (without colname) describes the chip names, and the last ones correspond the source of the biological samples hybridised to them and the GSE accession number of each chip respectively (with column names). You need to put this file into the "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/data" directory.

 b. For RNA and miRNA sequencing data analyses put the text files containing raw read count and biological phenotype data frames from TCGA or other sources into the working directory ("C:/Users/aaaa/bbbb/cccc/dddd"). The row names and the column names of raw read count data frames correspond to mRNAs/mature miRNA strands and Samples respectively. The row names and column names of biological phenotype data contain the samples and the biological characteristics of samples respectively.




#C. USAGE
........................................................................................................................................  
1. qc_microarray_data (GSE_acc, dir_download_GSES )
 
Description: This function performs some quality control checks in order to make sure that there are no issues with a specific GSE dataset. 

Input arguments 
GSE_acc: a character string indicating GSE number
dir_download_GSES: a character string indicating the pathway of directory of downloaded GSE dataset(s).

Output data files

This function generates a new directory "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/Quality_Control" containing .png files of boxplots and histograms of probe intensities before and after normalisation procedure (fRMA), pseudo-images of all microarray chips depicting the weights of probe level model fitting, NUSE and RLE boxplots along with their statistics in .txt files


2. remove_lq_samples (GSE_acc, cels_REMOVE_AFTER_QC,dir_download_GSES )
  
This function removes the low quality samples.

   
Input arguments 
    
GSE_acc: a character string indicating GSE number
cels_REMOVE_AFTER_QC: a character vector of names of .CEL files going to be removed
dir_download_GSES: a character string indicating the pathway of directory of downloaded GSE dataset(s).

Output data files

This function generates a new directory "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/Quality_Control" which contains the .CEL files of high quality microchips with the corresponding phenotype data frame in txt file and the frma normalized expression set in .RData file.



3. GSEs_merge_gene_filter (GSE_acc_4_merge, targetAnnot, batchAnnot,  transform_method= "COMBAT", dir_download_GSES )


This function merges multiple GSE datasets through batch effect removal, filter out uninformative data such as control probesets and other internal controls as well as removing genes with low variance
    
    
Input arguments 
GSE_acc_4_merge: a character vector of GSE accession numbers of different datasets to be merged.
targetAnnot: a character string of phenotype dataset column name corresponding the biological variable
batchAnnot: a character string of phenotype dataset column name corresponding the GSE accession number of each chip.
transform_method : Merging method aimed at removing batch effects. Possible options are: BMC, COMBAT, DWD, GENENORM, GENESHIFT, NONE and XPN
dir_download_GSES: a character string indicating the pathway of directory of downloaded GSE datasets

    
Output data files

GSEs_merge_gene_filter function generates four .png files including MDSs and dendrograms for visual inspection of the merged datasets with and without applying any transformation method, four .txt files containing merged expression and merged pheno data frames with and without applying tranformation, and two . RData files containing merged expression set with and without transformation. 




4. DEA_microarrays (Exprs_table, Pheno_data_table,sample_type_colname, diff_exp_type,adjust_pvalue_method, LogFC_cuttoff, Adj_p_value_cuttoff, annotation_microarray_probesets, outputFileFolder )

This function performs differential expression analysis for microarray data.

Input arguments
     
Exprs_table: Expression data frame containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples.
Pheno_data_table: Phenotype dataframe containing biological information of samples with rows corresponding to samples and columns to biological variables
sample_type_colname: a character string of phenotype dataset column name corresponding to a biological variable
diff_exp_type : a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")
adjust_pvalue_method: character string specifying p-value adjustment method. Possible options: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.
LogFC_cuttoff: log2 fold change threshold value
Adj_p_value_cuttoff: adjusted p-value threshold
annotation_microarray_probesets: logical value. If TRUE it will be performed annotation of the results with associated gene symbols
outputFileFolder: a character string indicating the output directory pathway

Output data files

DEA_microarrays function generates a .txt file containing the table of the top-ranked features from a linear model fit and a .png file including a volcano plot of log-fold changes versus log-odds of differential expression

    
5. cluster_heatmap (Exprs_table, Pheno_data_table, vars2study, varType, dist_method,clust_method,main_title, outputFileFolder)

This function generates an enhanced heatmap representation with a dendrogram added to the left side and to the top, according to cluster analysis.

Input arguments

Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
vars2study: a character vector of column names of Pheno_data_table object to be studied
varType : a character vector indicating the type of variables vars2study. Possible options: cat (categorical variable type), ord (ordinal variable type)
dist_method : the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"  
clust_method : the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
main_title : a character string indicating an overall title for the plot
outputFileFolder :  a character string indicating the output directory pathway


Output data files

Cluster_heatmap function generates a .png file containing the clustering heatmap plot





  6. Boxplots (Exprs_table, Pheno_data_table, var2study,group_names, genes_or_mirnas,addDot=TRUE, x_label, y_label, outputFileFolder)
  
This function produces boxplot(s) of the given grouped values.


Input arguments

Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
var2study: a character string of phenotype dataset column name corresponding to a biological variable to be studied
group_names: a character vector which contains the group names of biological variable to be studied
genes_or_mirnas: a character vector which contains feature (mRNA,mature miRNA) names including in expression data frame to be studied
addDot: logical value. If TRUE dots are added to boxplots
x_label: a character string indicating x-axis title
y_label: a character string indicating y-axis title
outputFileFolder: a character string indicating the output directory pathway


Output data files

For each feature (mRNA/mature miRNA), Boxplots function generates .png file containing boxplot of the given biologically grouped expression values.





7. prepr_common_miRNA_mRNA (RNA_seq_data, miRNA_seq_data,TCGA_source, outputFileFolder)

   
This function constructs miRNA- and RNA- seq raw read count data frames with common samples.

Input arguments

RNA_seq_data: a data frame containing RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples
miRNA_seq_data: a data frame containing miRNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples
TCGA_source: logical value. If TRUE, function will be able to manage TCGA data
outputFileFolder: a character string indicating the output directory pathway
  
Output data files
     
If TCGA_source=FALSE prepr_common_miRNA_mRNA generates two .txt files containing miRNA- and RNA- seq raw read count data frames with common samples respectively. By defining TCGA_source= TRUE it further creates two .txt files corresponding to phenotypic information of miRNA- and RNA- seq respectively. 
   


8. Prepr_DEA_mRNA_miRNA_seq (Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type,data_type, adjust_pvalue_method, LogFC_cuttoff, Adj_p_value_cuttoff, outputFileFolder)

  
This function filters lowly expressed features and applies TMM normalization, voom transformation as well as differential expression analysis of miRNA-, RNA- seq data.

Input arguments

Exprs_table: a data frame containing miRNA- or RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
sample_type_colname: a character string of phenotype dataset column name corresponding to a biological variable
diff_exp_type: a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")
data_type: a character string including the type of data. Possible options: mRNA, miRNA
adjust_pvalue_method: character string specifying p-value adjustment method. Possible options: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.
LogFC_cuttoff: log2 fold change threshold value
Adj_p_value_cuttoff: adjusted p-value threshold
outputFileFolder: a character string indicating the output directory pathway

Output data files

Prepr_DEA_mRNA_miRNA_seq function generates a .txt file containing the table of the top-ranked features (mature miRNAs/mRNAs) from a linear model fit and a .png file including a volcano plot of log-fold changes versus log-odds of differential expression



9. Feature_cor_net (NORM_EXPR_VALUES, feature_list,feature_type,cor_pval_cutoff, outputFileFolder)
  
This function performs correlation (co-expression) network analysis.


Input arguments

NORM_EXPR_VALUES: a data frame containing normalized expression values for a series of samples, with rows corresponding to features (mRNA, mature miRNAs) and columns to samples
feature_list: a character vector containing the features (mRNAs, mature miRNAs) to be studied.
feature_type: type of data. Possible options: mRNA, miRNA
cor_pval_cutoff: a correlation p-value threshold to define the most significant interactions
outputFileFolder: a character string indicating the output directory pathway


Output data files

Feature_cor_net function generates three .txt files. The first file contains a data frame with the correlated feature targets for each feature. The second file contains a data frame organizing the most statistically significant interactions in pairs with their corresponding statistics and the third one a data frame including the node attributes. The last two files can be opened in Cytoscape in order to visualize the correlation network.

  
10. cor_heatmap (NORM_EXPR_VALUES, features_of_interest, outputFileFolder)

For user defined feature, this function performs a graphical display of correlation matrices (correlation coefficients) of positive and and negative correlations of this feature with all other features respectively.

Input arguments
     
NORM_EXPR_VALUES: a data frame containing normalized expression values for a series of samples, with rows corresponding to features (mRNA, mature miRNAs) and columns to samples
features_of_interest: a charaster vector containing the features (mRNAs, mature miRNAs) to be studied
outputFileFolder: a character string indicating the output directory pathway

    
Output data files

For each feature of interest cor_heatmap function generates two .png files including the heatmap displays of negative and positive correlations of this feature with all others, as well as two .txt files which contain the data frames organizing the results in pairs with their corresponding statistics.
 


11. miRNA_mRNA_inters (NORM_EXPR_VALUES_MRNA,NORM_EXPR_VALUES_MIRNA,gene_list, miRNA_list,method= "pearson",adjust_pvalue_method= NULL, miRNA_annot, cor_pvalue_cutoff,outputFileFolder) 

This function performs miRNA-mRNA interaction network analysis combining expression profiles with target site information (targetScan and microCosm databases)

Input arguments

NORM_EXPR_VALUES_MRNA: a data frame containing normalized expression values for a series of samples, with rows corresponding to mRNAs and columns to samples
NORM_EXPR_VALUES_MIR: a data frame containing normalized expression values for a series of samples, with rows corresponding to mature miRNAs and columns to samples
gene_list: a character vector containing mRNAs to be included in the network
miRNA_list: a character vector containing miRNAs to be included in the network
method: a character string indicating which correlation coefficient is to be used for the test. One of pearson, kendall, or spearman, can be abbreviated
adjust_pvalue_method: character string specifying p-value adjustment method. Possible options: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.
miRNA_annot: a character string indicating the data frame of annotation of mature miRNAs. The first column contains the mature miRNA accession numbers (MIMAT SYMBOLS) and the second one the mature miRNA symbols.
cor_pvalue_cutoff: a correlation p-value threshold to define the most significant miRNA-mRNA interactions
outputFileFolder: a character string indicating the output directory pathway


Output data files

miRNA_mRNA_inters function generates four .txt files. The first file contains a data frame showing the targets for each miRNA. The second file contains a dataframe including the number of miRNAs that are targeting a specific mRNA. The third file contains a data frame organizing the most statistically significant miRNA-mRNA interactions in pairs with their corresponding statistics and the last one a data frame including the node attributes. The last two files can be opened in Cytoscape in order to visualize the miRNA-mRNA correlation network.



12. Dotchart (Exprs_table, Pheno_data_table, var2study,group_names, genes_or_mirnas,  x_label, y_label, outputFileFolder)

This function draws a Cleveland dot plot

    
Input arguments
Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
var2study: a character string of phenotype dataset column name corresponding to a biological variable to be studied
group_names: a character vector which contains the group names of biological variable to be studied
genes_or_mirnas: a character vector which contains feature (mRNA,mature miRNA) names including in expression data frame to be studied
x_label: a character string indicating x-axis title
y_label: a character string indicating y-axis title
outputFileFolder: a character string indicating the output directory pathway

Output data files

For each feature (mRNA/mature miRNA), Dotchart function generates .png file containing a Cleveland dot plot of the given biologically grouped expression values.




13. PCA (Exprs_table, Pheno_data_table, pheno_var2study,PCs, labels, circle, ellipse, var.axes, outputFileFolder)

This function performs Principal component analysis (PCA)

      
Input arguments

Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
pheno_var2study: a character string of phenotype dataset column name corresponding to a biological variable to be studied
PCs: a numeric vector containing the PCs to be ploted
labels: logical value. If TRUE it labels the observations
circle: logical value. If TRUE it draws a correlation circle
ellipse: logical value. If TRUE it draws a normal data ellipse for each group
var.axes: logical value. If TRUE it draws arrows for the variables (mRNAs, miRNAs).
outputFileFolder: a character string indicating the output directory pathway

Output data files

PCA function generates two .png files containing a plot of the variances (y-axis) associated with the PCs (x-axis) and a biplot of principal components.



14. obs_strat_by_expr (Exprs_table, Pheno_data_table, gene_mirnas, outputFileFolder)

According the expression level for each feature (mRNA, miRNA) of interest, this function divides the observations (samples) into 3 groups (high, intermediate, low).

Input arguments

Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
gene_mirnas: a character vector of features(miRNAs, mRNAs) to be used for splitting of observations
outputFileFolder: a character string indicating the output directory pathway

Output data files

obs_strat_by_expr function generates a .txt file containing a data frame with rows corresponding to samples and, for each feature of interest, a column including the expression groups of samples. 




15. Filter_normaliz_mRNA_miRNA_seq_data (Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type,data_type, outputFileFolder )
  
This function filters lowly expressed features (mRNAs/miRNAs), applies TMM normalization and voom transformation as well as calculates log cpm expression counts of miRNA-, RNA- seq data.

Input arguments

Exprs_table: a data frame containing miRNA- or RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples
Pheno_data_table: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables
sample_type_colname: a character string of phenotype dataset column name corresponding to a biological variable
diff_exp_type: a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")
data_type: a character string including the type of data. Possible options: mRNA, miRNA
outputFileFolder: a character string indicating the output directory pathway


Output data files
  
Filter_normaliz_mRNA_miRNA_seq_data function generates two .txt files which contain two data frames including log-CPM expression counts and biological phenotype data respectively. 
  



16.  km_plot (Exprs_table, time_event_data, genes_or_mirnas,xlabel, ylabel, outputFileFolder )

This function performs a univariate Kaplan-Meier survival analysis. Specifically according to the expression level of each feature of interest it divides the samples into two groups (low, high) and it performs survival analysis between these groups.


Input arguments
     
Exprs_table: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples
time_event_data: a data frame of survival data with rows corresponding to samples and columns to time (column name= time) and event (column name= event) data respectively
genes_or_mirnas: a character vector of features (mRNAs, miRNAs) to be studied
xlabel: a character string indicating x-axis title
ylabel: a character string indicating y-axis title
outputFileFolder: a character string indicating the output directory pathway

 
Output data files

For each feature (mRNA, miRNA) km_plot function generates a .png file containing Kaplan-Meier survival curves for high and low expression sample groups


