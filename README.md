

Tools for RNA-seq/miRNA-seq/RNA microarray data analysis
============================================================
#A. REQUIREMENTS


* 1. Windows 10
* 2. R version 3.2.1 (2015-06-18)
* 3. R packages: affy, simpleaffy, frma, RColorBrewer, affyPLM, inSilicoMerging, amap, genefilter, limma, hgu133plus2.db, annotate, gplots, devtools, easyGgplot2, edgeR, corrplot, reshape2, ggplot2, ggbiplot

#B. INSTALLATION AND SET UP<br/>

* 1. Start R and install R/BioConductor dependencies<br/>  

```{r}
#download the BioC installation routines
source("http://bioconductor.org/biocLite.R")
#install the core packages. It will take some time!!
biocLite(c())
#install the Bioconductor packages
biocLite(c("affy", "simpleaffy", "frma", "affyPLM", "inSilicoMerging", "genefilter", "limma", "hgu133plus2.db", "annotate", "edgeR"))
#install the CRAN packages
install.packages(c("RColorBrewer", "amap", "gplots", "devtools", "corrplot", "reshape2", "ggplot2" ))
#install the Github packages
library(devtools)
install_github(c("kassambara/easyGgplot2", "vqv/ggbiplot"))
```
* 2. Download files from 'scripts' and 'internal_data' directory to a preference directory. To set this directory as your working directory of R, type:<br/>

```{r}
WD <- "C:/Users/aaaa/bbbb/cccc/dddd"
setwd(WD)
```
To load the R functions (scripts) type:<br/>
```{r}
functions_list <- list.files(pattern = "script")
for(i in 1:length(functions_list)){
    source(functions_list[i])
  }
  ```
To extract files from zip archives type:<br/>
```{r}
unzip("MODIFIED_targetScan_v6.2_18.txt.zip")
unzip("MODIFIED_microCosm_v5_18.zip")
 ```

* 3.a For GPL570 platform -Affymetrix Human Genome U133 Plus 2.0 microarray transcriptomic data- analyses you need to acquire the raw data directly from GEO. Specifically for each GSE accession number (GSEXXXX) you need to download the .CEL files and put them into a new directory named "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/data". Also you need to capture the experimental information. This is just a text file which includes a data frame containing three columns. The first column (without colname) describes the chip names, and the last ones correspond the source of the biological samples hybridised to them and the GSE accession number of each chip respectively (with column names). You need to put this file into the "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/data" directory.<br/>

* b. For RNA and miRNA sequencing data analyses put the text files containing raw read count and biological phenotype data frames from TCGA or other sources into the working directory ("C:/Users/aaaa/bbbb/cccc/dddd"). The row names and the column names of raw read count data frames correspond to mRNAs/mature miRNA strands and Samples respectively. The row names and column names of biological phenotype data contain the samples and the biological characteristics of samples respectively.<br/>




#C. USAGE

 ```{r} 
qc_microarray_data (GSE_acc, dir_download_GSES)
 ```
__Description__

This function performs some quality control checks in order to make sure that there are no issues with a specific GSE dataset.<br/> 

__Input arguments__

_GSE_acc_: a character string indicating GSE number<br/>
_dir_download_GSES_: a character string indicating the pathway of directory of downloaded GSE dataset(s)<br/>

__Output data files__

This function generates a new directory "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/Quality_Control" containing .png files of boxplots and histograms of probe intensities before and after normalisation procedure (fRMA), pseudo-images of all microarray chips depicting the weights of probe level model fitting, NUSE and RLE boxplots along with their statistics in .txt files.<br/>
<br/>
<br/>
 ```{r} 
remove_lq_samples (GSE_acc, cels_REMOVE_AFTER_QC, dir_download_GSES)
  ```
__Description__

This function removes the low quality samples.

   
__Input arguments__ 
    
_GSE_acc_: a character string indicating GSE number<br/>
_cels_REMOVE_AFTER_QC_: a character vector of names of .CEL files going to be removed<br/>
_dir_download_GSES_: a character string indicating the pathway of directory of downloaded GSE dataset(s)<br/>

__Output data files__

This function generates a new directory "C:/Users/aaaa/bbbb/cccc/dddd/GSEXXXX/Quality_Control" which contains the .CEL files of high quality microchips with the corresponding phenotype data frame in txt file and the frma normalized expression set in .RData file.<br/>

<br/>
<br/>
```{r}
GSEs_merge_gene_filter (GSE_acc_4_merge, targetAnnot, batchAnnot, transform_method, dir_download_GSES)
```

__Description__

This function merges multiple GSE datasets through batch effect removal, filter out uninformative data such as control probesets and other internal controls as well as removing genes with low variance.<br/>
    
    
__Input arguments__

_GSE_acc_4_merge_: a character vector of GSE accession numbers of different datasets to be merged<br/>
_targetAnnot_: a character string of phenotype dataset column name corresponding the biological variable<br/>
_batchAnnot_: a character string of phenotype dataset column name corresponding the GSE accession number of each chip<br/>
_transform_method_: Merging method aimed at removing batch effects. Possible options are: "BMC", "COMBAT", "DWD", "GENENORM", "GENESHIFT", "NONE" and "XPN"<br/>
_dir_download_GSES_: a character string indicating the pathway of directory of downloaded GSE datasets<br/>

    
__Output data files__

GSEs_merge_gene_filter function generates four .png files including MDSs and dendrograms for visual inspection of the merged datasets with and without applying any transformation method, four .txt files containing merged expression and merged pheno data frames with and without applying tranformation, and two . RData files containing merged expression set with and without transformation.<br/>


<br/>
<br/>
```{r}
DEA_microarrays (Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type, adjust_pvalue_method, LogFC_cuttoff, Adj_p_value_cuttoff, annotation_microarray_probesets, outputFileFolder)
```
__Description__

This function performs differential expression analysis for microarray data.<br/>

__Input arguments__
     
_Exprs_table_: Expression data frame containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples<br/>
_Pheno_data_table_: Phenotype dataframe containing biological information of samples with rows corresponding to samples and columns to biological varables<br/>
_sample_type_colname_: a character string of phenotype dataset column name corresponding to a biological variable<br/>
_diff_exp_type_: a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")<br/>
_adjust_pvalue_method_: character string specifying p-value adjustment method. Possible options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"<br/>
_LogFC_cuttoff_: log2 fold change threshold value<br/>
_Adj_p_value_cuttoff_: adjusted p-value threshold<br/>
_annotation_microarray_probesets_: logical value. If TRUE it will be performed annotation of the results with associated gene symbols
_outputFileFolder_: a character string indicating the output directory pathway<br/>

__Output data files__

DEA_microarrays function generates a .txt file containing the table of the top-ranked features from a linear model fit and a .png file including a volcano plot of log-fold changes versus log-odds of differential expression.<br/>
<br/>
<br/>
```{r}    
cluster_heatmap (Exprs_table, Pheno_data_table, vars2study, varType, dist_method, clust_method, main_title, outputFileFolder)
```

__Description__

This function generates an enhanced heatmap representation with a dendrogram added to the left side and to the top, according to cluster analysis.<br/>

__Input arguments__

_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_vars2study_: a character vector of column names of Pheno_data_table object to be studied<br/>
_varType_ : a character vector indicating the type of variables vars2study. Possible options: "cat" (categorical variable type), "ord" (ordinal variable type)<br/>
_dist_method_ : the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"<br/>
_clust_method_: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)<br/>
_main_title_: a character string indicating an overall title for the plot<br/>
_outputFileFolder_:  a character string indicating the output directory pathway<br/>


__Output data files__

Cluster_heatmap function generates a .png file containing the clustering heatmap plot.<br/>



<br/>
<br/>
```{r}  
Boxplots (Exprs_table, Pheno_data_table, var2study, group_names, genes_or_mirnas, addDot, x_label, y_label, outputFileFolder)
  ```
__Description__

This function produces boxplot(s) of the given grouped values.<br/>


__Input arguments__

_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_var2study_: a character string of phenotype dataset column name corresponding to a biological variable to be studied<br/>
_group_names_: a character vector which contains the group names of biological variable to be studied<br/>
_genes_or_mirnas_: a character vector which contains feature (mRNA,mature miRNA) names including in expression data frame to be studied<br/>
_addDot_: logical value. If TRUE dots are added to boxplots<br/>
_x_label_: a character string indicating x-axis title<br/>
_y_label_: a character string indicating y-axis title<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>


__Output data files__

For each feature (mRNA/mature miRNA), Boxplots function generates .png file containing boxplot of the given biologically grouped expression values.<br/>


<br/>
<br/>

```{r}  
prepr_common_miRNA_mRNA (RNA_seq_data, miRNA_seq_data, TCGA_source, outputFileFolder)
```
   
__Description__

This function constructs miRNA- and RNA- seq raw read count data frames with common samples.<br/>

__Input arguments__

_RNA_seq_data_: a data frame containing RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples<br/>
_miRNA_seq_data_: a data frame containing miRNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples<br/>
_TCGA_source_: logical value. If TRUE, function will be able to manage TCGA data<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>
  
__Output data files__
     
If TCGA_source=FALSE prepr_common_miRNA_mRNA generates two .txt files containing miRNA- and RNA- seq raw read count data frames with common samples respectively. By defining TCGA_source= TRUE it further creates two .txt files corresponding to phenotypic information of miRNA- and RNA- seq respectively.<br/>
   
<br/>
<br/>
```{r}  
Prepr_DEA_mRNA_miRNA_seq (Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type,data_type, adjust_pvalue_method, LogFC_cuttoff, Adj_p_value_cuttoff, outputFileFolder)
```

  
__Description__

This function filters lowly expressed features and applies TMM normalization, voom transformation as well as differential expression analysis of miRNA-, RNA- seq data.<br/>

__Input arguments__

_Exprs_table_: a data frame containing miRNA- or RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_sample_type_colname_: a character string of phenotype dataset column name corresponding to a biological variable<br/>
_diff_exp_type_: a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")<br/>
_data_type_: a character string including the type of data. Possible options: "mRNA", "miRNA"<br/>
_adjust_pvalue_method_: character string specifying p-value adjustment method. Possible options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"<br/>
_LogFC_cuttoff_: log2 fold change threshold value<br/>
_Adj_p_value_cuttoff_: adjusted p-value threshold<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>

__Output data files__<br/>

Prepr_DEA_mRNA_miRNA_seq function generates a .txt file containing the table of the top-ranked features (mature miRNAs/mRNAs) from a linear model fit and a .png file including a volcano plot of log-fold changes versus log-odds of differential expression.<br/>

<br/>
<br/>
```{r}  
Feature_cor_net (NORM_EXPR_VALUES, feature_list, feature_type, cor_pval_cutoff, outputFileFolder)
```
  
__Description__

This function performs correlation (co-expression) network analysis.<br/>


__Input arguments__

_NORM_EXPR_VALUES_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features (mRNA, mature miRNAs) and columns to samples<br/>
_feature_list_: a character vector containing the features (mRNAs, mature miRNAs) to be studied<br/>
_feature_type_: type of data. Possible options: "mRNA", "miRNA"<br/>
_cor_pval_cutoff_: a correlation p-value threshold to define the most significant interactions<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>


__Output data files__<br/>

Feature_cor_net function generates three .txt files. The first file contains a data frame with the correlated feature targets for each feature. The second file contains a data frame organizing the most statistically significant interactions in pairs with their corresponding statistics and the third one a data frame including the node attributes. The last two files can be opened in Cytoscape in order to visualize the correlation network.<br/>
<br/>
<br/>
```{r}    
cor_heatmap (NORM_EXPR_VALUES, features_of_interest, outputFileFolder)
```

__Description__

For user defined feature, this function performs a graphical display of correlation matrices (correlation coefficients) of positive and and negative correlations of this feature with all other features respectively.<br/>

__Input arguments__
     
_NORM_EXPR_VALUES_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features (mRNA, mature miRNAs) and columns to samples<br/>
_features_of_interest_: a charaster vector containing the features (mRNAs, mature miRNAs) to be studied<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>

    
__Output data files__

For each feature of interest cor_heatmap function generates two .png files including the heatmap displays of negative and positive correlations of this feature with all others, as well as two .txt files which contain the data frames organizing the results in pairs with their corresponding statistics.<br/>
 
<br/>
<br/>
```{r}  
miRNA_mRNA_inters (NORM_EXPR_VALUES_MRNA, NORM_EXPR_VALUES_MIRNA, gene_list, miRNA_list, method, adjust_pvalue_method, miRNA_annot,  cor_pvalue_cutoff, outputFileFolder)
```

__Description__

This function performs miRNA-mRNA interaction network analysis combining expression profiles with target site information (targetScan and microCosm databases).<br/>

__Input arguments__

_NORM_EXPR_VALUES_MRNA_: a data frame containing normalized expression values for a series of samples, with rows corresponding to mRNAs and columns to samples<br/>
_NORM_EXPR_VALUES_MIR_: a data frame containing normalized expression values for a series of samples, with rows corresponding to mature miRNAs and columns to samples<br/>
_gene_list_: a character vector containing mRNAs to be included in the network<br/>
_miRNA_list_: a character vector containing miRNAs to be included in the network<br/>
_method_: a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated<br/>
_adjust_pvalue_method_: character string specifying p-value adjustment method. Possible options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"<br/>
_miRNA_annot_: a character string indicating the data frame of annotation of mature miRNAs. The first column contains the mature miRNA accession numbers (MIMAT SYMBOLS) and the second one the mature miRNA symbols<br/>
_cor_pvalue_cutoff_: a correlation p-value threshold to define the most significant miRNA-mRNA interactions<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>


__Output data files__

miRNA_mRNA_inters function generates four .txt files. The first file contains a data frame showing the targets for each miRNA. The second file contains a dataframe including the number of miRNAs that are targeting a specific mRNA. The third file contains a data frame organizing the most statistically significant miRNA-mRNA interactions in pairs with their corresponding statistics and the last one a data frame including the node attributes. The last two files can be opened in Cytoscape in order to visualize the miRNA-mRNA correlation network.<br/>

<br/>
<br/>
```{r}  
Dotchart (Exprs_table, Pheno_data_table, var2study, group_names, genes_or_mirnas, x_label, y_label, outputFileFolder)
```

__Description__

This function draws a Cleveland dot plot.<br/>

    
__Input arguments__

_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_var2study_: a character string of phenotype dataset column name corresponding to a biological variable to be studied<br/>
_group_names_: a character vector which contains the group names of biological variable to be studied<br/>
_genes_or_mirnas_: a character vector which contains feature (mRNA,mature miRNA) names including in expression data frame to be studied<br/>
_x_label_: a character string indicating x-axis title<br/>
_y_label_: a character string indicating y-axis title<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>

__Output data files__

For each feature (mRNA/mature miRNA), Dotchart function generates .png file containing a Cleveland dot plot of the given biologically grouped expression values.<br/>


<br/>
<br/>
```{r}  
PCA (Exprs_table, Pheno_data_table, pheno_var2study, PCs, labels, circle, ellipse, var.axes, outputFileFolder)
```

__Description__

This function performs Principal component analysis (PCA).<br/>

      
__Input arguments__

_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_pheno_var2study_: a character string of phenotype dataset column name corresponding to a biological variable to be studied<br/>
_PCs_: a numeric vector containing the PCs to be ploted<br/>
_labels_: logical value. If TRUE it labels the observations<br/>
_circle_: logical value. If TRUE it draws a correlation circle<br/>
_ellipse_: logical value. If TRUE it draws a normal data ellipse for each group<br/>
_var.axes_: logical value. If TRUE it draws arrows for the variables (mRNAs, miRNAs)<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>

__Output data files__

PCA function generates two .png files containing a plot of the variances (y-axis) associated with the PCs (x-axis) and a biplot of principal components.<br/>

<br/>
<br/>
```{r}  
obs_strat_by_expr (Exprs_table, Pheno_data_table, gene_mirnas, outputFileFolder)
```

__Description__

According the expression level for each feature (mRNA, miRNA) of interest, this function divides the observations (samples) into 3 groups (high, intermediate, low).<br/>

__Input arguments__

_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_gene_mirnas_: a character vector of features(miRNAs, mRNAs) to be used for splitting of observations<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>

__Output data files__

obs_strat_by_expr function generates a .txt file containing a data frame with rows corresponding to samples and, for each feature of interest, a column including the expression groups of samples.<br/>


<br/>
<br/>
```{r}  
Filter_normaliz_mRNA_miRNA_seq_data (Exprs_table, Pheno_data_table, sample_type_colname, diff_exp_type, data_type, outputFileFolder)
```
  
__Description__

This function filters lowly expressed features (mRNAs/miRNAs), applies TMM normalization and voom transformation as well as calculates log cpm expression counts of miRNA-, RNA- seq data.<br/>

__Input arguments__

_Exprs_table_: a data frame containing miRNA- or RNA-seq raw read counts for a series of samples, with rows corresponding to features and columns to samples<br/>
_Pheno_data_table_: a data frame containing biological information of samples with rows corresponding to samples and columns to biological variables<br/>
_sample_type_colname_: a character string of phenotype dataset column name corresponding to a biological variable<br/>
_diff_exp_type_: a character string including group names of samples to be studied separated by "_vs_" (e.g. "GroupA_vs_GroupB")<br/>
_data_type_: a character string including the type of data. Possible options: "mRNA", "miRNA"<br/>
_outputFileFolder_: a character string indicating the output directory pathway<br/>


__Output data files__
  
Filter_normaliz_mRNA_miRNA_seq_data function generates two .txt files which contain two data frames including log-CPM expression counts and biological phenotype data respectively.<br/>
  

<br/>
<br/>
```{r}  
km_plot (Exprs_table, time_event_data, genes_or_mirnas, xlabel, ylabel, outputFileFolder)
```

__Description__

This function performs a univariate Kaplan-Meier survival analysis. Specifically according to the expression level of each feature of interest it divides the samples into two groups (low, high) and it performs survival analysis between these groups.<br/>


__Input arguments__
     
_Exprs_table_: a data frame containing normalized expression values for a series of samples, with rows corresponding to features and columns to samples<br/>
_time_event_data_: a data frame of survival data with rows corresponding to samples and columns to time (column name= time) and event (column name= event) data respectively<br/>
_genes_or_mirnas_: a character vector of features (mRNAs, miRNAs) to be studied<br/>
_xlabel_: a character string indicating x-axis title<br/>
_ylabel_: a character string indicating y-axis title<br/>
outputFileFolder: a character string indicating the output directory pathway<br/>

 
__Output data files__

For each feature (mRNA, miRNA) km_plot function generates a .png file containing Kaplan-Meier survival curves for high and low expression sample groups.<br/>

<br/>
<br/>


Copyright (c) 2016 AUTH<br/>
Author: Panagiotis Mokos<br/>

