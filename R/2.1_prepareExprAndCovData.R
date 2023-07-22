


prepareExprAndCovData<-function(exprFile,covFile){
  # tissueType<-"Lung" #Sample size is 515.
  # numOfPCs<-44 #Chosen via BE.
  #
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
  #                 tissueType,"/",numOfPCs,"ExprPCs.txt")
  # #To comment out.

  dataGeneExpressionFP<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  colnames(dataGeneExpressionFP)[1:4]<-c("chr","start","end","gene_id") #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.

  dataCovariatesRaw<-readr::read_delim(covFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #52*516. The first column is the covariate ID. Import from text (base) does not work. Import from text (readr) works.
  if(!identical(colnames(dataGeneExpressionFP)[-(1:4)],colnames(dataCovariatesRaw)[-1])){
    stop("Sample names do not match between the expression file and the covariate file.\n")
  }
  dataCovariates<-t(dataCovariatesRaw[,-1]) #515*52. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.
  colnames(dataCovariates)<-unlist(dataCovariatesRaw[,1]) #This is unnecessary but may be helpful during code development.

  #Filter out constant covariates, e.g., sex in sex organs.
  # dataCovariates[,7]<-1 #Test print message.
  # dataCovariates[,8]<-1 #Test print message.
  vars<-apply(dataCovariates,2,var)
  if(any(vars==0)){
    cat("Index (or indices) of constant covariate(s): ")
    cat(which(vars==0),sep=", ")
    cat(". Constant covariates are filtered out for the analysis.\n")

    dataCovariates<-dataCovariates[,which(vars!=0)] #515*52. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.
  }

  toReturn<-list(dataGeneExpressionFP=dataGeneExpressionFP,dataCovariates=dataCovariates)
  return(toReturn)
}




