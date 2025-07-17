#' @title
#' callSigGeneSNPPairs
#'
#' @description
#' This function is used for identifying significant gene-SNP pairs.
#'
#' @details
#' Most arguments in this function (from \code{exprFile} to \code{MASamplesThreshold}) must be the same as those used in \code{ClipperQTL()}. However, since \code{outputDir} will be used in C++ in this function, it must be recognizable in C++ (see above).
#'
#' The main method parameters of \code{callSigGeneSNPPairs()} are \code{sigGeneSNPPairMethod} and \code{percent}. \code{sigGeneSNPPairMethod} must be \code{"FastQTL"}, \code{"topPercent"}, or \code{NULL}. \code{percent} is only relevant if \code{"topPercent"} is used.
#'
#' If \code{sigGeneSNPPairMethod="FastQTL"}, then significant gene-SNP pairs will be identified using the method in FastQTL (see Algorithm S2 of reference; inverse functions of cumulative distribution functions are replaced by quantile functions).
#'
#' If \code{sigGeneSNPPairMethod="topPercent"}, then the top 1\% local common SNPs (in terms of significance of association) will be identified as significant for each eGene (if \code{percent=1}, that is).
#'
#' If \code{sigGeneSNPPairMethod=NULL}, then if \code{approach="standard"}, \code{sigGeneSNPPairMethod} will be set to \code{"FastQTL"}; if \code{approach="Clipper"}, \code{sigGeneSNPPairMethod} will be set to \code{"topPercent"}.
#'
#' If \code{approach="standard"}, then \code{sigGeneSNPPairMethod} can be either \code{"FastQTL"} or \code{"topPercent"}. If \code{approach="Clipper"}, then \code{sigGeneSNPPairMethod} cannot be \code{"FastQTL"}.
#'
#' \code{callSigGeneSNPPairs()} outputs several files in the output directory. The most important one is named "_resultSigGeneSNPPairs_FastQTL.rds" or "_resultSigGeneSNPPairs_topPercent.rds" (depending on \code{sigGeneSNPPairMethod}), which can be read into R with \code{readRDS()}. Each row corresponds to a significant gene-SNP pair.
#'
#' @param exprFile Must be the same as \code{exprFile} in \code{ClipperQTL()}.
#' @param covFile Must be the same as \code{covFile} in \code{ClipperQTL()}.
#' @param genotypeFile Must be the same as \code{genotypeFile} in \code{ClipperQTL()}.
#' @param tabixProgram Must be the same as \code{tabixProgram} in \code{ClipperQTL()}.
#' @param outputDir Must be the same as \code{outputDir} in \code{ClipperQTL()}. Must be recognizable in C++. "~" may not be recognized as the home directory in C++; the user may need to spell out the directory name instead.
#'
#' @param approach Must be the same as \code{approach} in \code{ClipperQTL()} (\code{"standard"} or \code{"Clipper"}). Helps determine \code{sigGeneSNPPairMethod}.
#'
#' @param cisDistance Must be the same as \code{cisDistance} in \code{ClipperQTL()}.
#' @param MAFThreshold Must be the same as \code{MAFThreshold} in \code{ClipperQTL()}.
#' @param MASamplesThreshold Must be the same as \code{MASamplesThreshold} in \code{ClipperQTL()}.
#'
#' @param numOfCores The number of cores to be used. Default is \code{1}.
#'
#' @param FDR_eGene The target FDR for eGene identification. Default is \code{0.05}.
#' @param sigGeneSNPPairMethod Key argument. Must be \code{"FastQTL"}, \code{"topPercent"}, or \code{NULL}. See below.
#' @param percent Key argument. Only relevant if \code{"topPercent"} is used. Default is \code{1}, which means the top 1\% local common SNPs (in terms of significance of association) will be identified as significant for each eGene. See below.
#'
#' @references
#' Heather J. Zhou, Xinzhou Ge, and Jingyi Jessica Li. ClipperQTL: ultrafast and powerful eGene identification method. Genome Biology, 26:207, 2025.
#'
#' @export





# #For code development only:
# library(dplyr)
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.5_prepareDataGeneExpressionFPSub.R")
#
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/4.1_prepareSigGeneSNPPairsParameters.R")



callSigGeneSNPPairs<-function(exprFile,covFile,genotypeFile,tabixProgram,outputDir, #Must be the same as those used in ClipperQTL().
                              approach, #Must be the same as approach in ClipperQTL(). Helps determine sigGeneSNPPairMethod.
                              cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10, #Must be the same as those used in ClipperQTL().
                              numOfCores=1,
                              FDR_eGene=0.05,sigGeneSNPPairMethod=NULL,percent=1){ #sigGeneSNPPairMethod must be "FastQTL", "topPercent", or NULL. percent is only relevant if "topPercent" is used.
  # tissueType<-"Lung" #Sample size is 515.
  # numOfPCs<-44 #Chosen via BE.
  #
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
  #                 tissueType,"/",numOfPCs,"ExprPCs.txt")
  # genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # outputDir<-paste0("/home/heatherjzhou/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/") #Key.
  # rm(tissueType,numOfPCs)
  #
  # approach<-"standard"
  #
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # numOfCores<-1
  #
  # FDR_eGene<-0.05
  # sigGeneSNPPairMethod<-NULL
  # percent<-1
  # #To comment out.



  #Prepare parameters for calling significant gene-SNP pairs (do this here so that sigGeneSNPPairMethod can be used to name the log file. Stop messages are not diverted by sink() anyways).
  temp<-prepareSigGeneSNPPairsParameters(approach,sigGeneSNPPairMethod,percent)
  sigGeneSNPPairMethod<-temp$sigGeneSNPPairMethod #"FastQTL" or "topPercent".
  percent<-temp$percent
  rm(temp)

  filename<-paste0(outputDir,"_log_sigGeneSNPPairs_",sigGeneSNPPairMethod,".txt")
  sink(file=filename,type="output") #Save print messages to this file. Messages sent to stderr() (including those from message, warning, and stop) do not go here.

  #Prepare dataGeneExpressionFP and dataCovariates.
  cat("Reading in expression data and covariate data...\n")
  temp<-prepareExprAndCovData(exprFile,covFile)
  dataGeneExpressionFP<-temp$dataGeneExpressionFP #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  dataCovariates<-temp$dataCovariates #515*52. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.
  rm(temp)

  #Prepare sampleIndices, a vector representing which columns of the genotype file correspond to samples in the expression and covariate data (need to subtract 1 before using in Cpp).
  cat("Preparing sample indices...\n")
  sampleIndices<-prepareSampleIndices(genotypeFile,tabixProgram,outputDir,
                                      sampleNames=colnames(dataGeneExpressionFP)[-(1:4)]) #Vector of length 515.

  chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.


  #Get eGenes.
  resultGenes<-readRDS(paste0(outputDir,"_resultGenes.rds"))
  eGenes<-(resultGenes%>%filter(qValue<=FDR_eGene))$gene_id #Vector of length 13,629.

  if(sigGeneSNPPairMethod=="FastQTL"){
    #Calculate pt.
    resultGenes_yes<-resultGenes%>%filter(gene_id%in%eGenes) #13,629*1007.
    resultGenes_no<-resultGenes%>%filter(!gene_id%in%eGenes) #12,466*1007.

    p1<-max(resultGenes_yes$pValue)
    p2<-min(resultGenes_no$pValue)
    pt<-(p1+p2)/2 #0.08641359.

    #Calculate absCorThreshold (a variable in resultGenes_yes), the threshold for calling significant gene-SNP pairs.
    #The p-value of a gene-SNP pair being in the top X (pt*100) percent is equivalent to the absolute correlation of the gene-SNP pair being in the top X percent, where "top" means the most significant.
    maxAbsCors<-as.matrix(resultGenes_yes[,-c(1:5,(ncol(resultGenes_yes)-1):ncol(resultGenes_yes))]) #13,629*1000. The columns are: bg1, ..., bg1000.
    # dim(maxAbsCors)
    resultGenes_yes$absCorThreshold<-apply(maxAbsCors,1,function(x){ #1 means by row. This takes a moment.
      return(quantile(x,probs=1-pt)) #Use 1-pt because larger absolute correlations are more significant.
    }) #13,629*1008. absCorThreshold is 0.1724693, 0.1767774, etc.

    # hist(maxAbsCors[1,],breaks=50)
    # hist(resultGenes_yes$absCorThreshold)
  }else{
    #Do nothing.
  }

  #Run chunks. This creates resultChunk1.rds, resultChunk2.rds, etc.
  if(TRUE){
    results<-parallel::mclapply(chunkInfo$indexOfChunk,FUN=function(indexOfChunk){ #parallel is a base package, so it doesn't need to be imported in the description file of the package.
      # indexOfChunk<-5 #To comment out.

      cat("\nRunning Chunk",indexOfChunk,"out of",nrow(chunkInfo),"chunks...\n") #Include a new line at the beginning of this print message to emphasize it.

      #Prepare dataGeneExpressionFPSub (eGenes only).
      dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
      dataGeneExpressionFPSub<-dataGeneExpressionFPSub%>%filter(gene_id%in%eGenes) #127*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
      if(nrow(dataGeneExpressionFPSub)==0) return(0) #If there is no eGene in this chunk, then do nothing.

      if(sigGeneSNPPairMethod=="FastQTL"){
        absCorThresholds<-(resultGenes_yes%>%filter(gene_id%in%dataGeneExpressionFPSub$gene_id))$absCorThreshold #Vector of length 127. 0.1838889, 0.1879149, etc.
      }else{
        absCorThresholds<-0 #Use 0 as a placeholder.
      }

      runChunk_sigGeneSNPPairs(dataGeneExpressionFPSub,dataCovariates,
                               genotypeFile,tabixProgram,sampleIndices,
                               cisDistance,MAFThreshold,MASamplesThreshold,
                               indexOfChunk,outputDir,
                               sigGeneSNPPairMethod,absCorThresholds,percent) #absCorThresholds is only relevant when using "FastQTL". percent is only relevant when using "topPercent".
      gc()
      return(0)
    },mc.cores=numOfCores) #results is a list of length nrow(chunkInfo). Each core takes up to a few GB of memory.

    # chunkInfo<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/Lung/_chunkInfo.rds")
    # results<-rep(0,nrow(chunkInfo))
    # results<-as.list(results)
    cat("\n") #Include a new line here to emphasize the print message.
    print(unlist(results)) #0 means success.
  }

  #Save _resultSigGeneSNPPairs.rds.
  cat("\nCombining results...\n") #Include a new line at the beginning of this print message to emphasize it.
  combineChunks_sigGeneSNPPairs(outputDir,chunkInfo,
                                n=length(sampleIndices),K=ncol(dataCovariates),
                                sigGeneSNPPairMethod)

  cat("\ncallSigGeneSNPPairs() finished running.\n") #Include a new line at the beginning of this print message to emphasize it.

  sink() #Close the connection.
}




