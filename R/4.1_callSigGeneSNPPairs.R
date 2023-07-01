#' callSigGeneSNPPairs
#'
#' @export


# #For code development only:
# library(dplyr)
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")


# ClipperQTL<-function(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
#                      approach="standard",B=1000,
#                      cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
#                      numOfChunksTarget=100,seed=1,numOfCores=1){


#This function: save _resultSigGeneSNPPairs.rds.
callSigGeneSNPPairs<-function(exprFile,covFile,genotypeFile,tabixProgram,outputDir, #Need to be the same as those used in ClipperQTL().
                              approach="standard", #Need to be the same as that used ClipperQTL().
                              cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10, #Need to be the same as those used in ClipperQTL().
                              numOfCores=1,
                              sigGeneSNPPairMethod=NULL,percent=NULL){ #sigGeneSNPPairMethod is "FastQTL" or "topPercent".
  # tissueType<-"Lung" #Sample size is 515.
  # numOfPCs<-44 #Chosen via BE.
  #
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
  #                 tissueType,"/",numOfPCs,"ExprPCs.txt")
  # genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  # rm(tissueType,numOfPCs)
  #
  # approach<-"Clipper" #Use this to make testing faster.
  # B<-20 #Use this to make testing faster.
  #
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # numOfChunksTarget<-100 #The total number of chunks will be approximately this number.
  # seed<-1
  # numOfCores<-1
  # #To comment out.



  filename<-paste0(outputDir,"_log_callSigGeneSNPPairs.txt")
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

  chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds"))

  cat("\nCalling significant gene-SNP pairs...\n") #Include a new line at the beginning of this print message to emphasize it.

  #Run chunks.
  if(TRUE){
    RNGkind("L'Ecuyer-CMRG") #This is necessary to ensure reproducibility when using mclapply().
    results<-parallel::mclapply(chunkInfo$indexOfChunk,FUN=function(indexOfChunk){ #parallel is a base package, so it doesn't need to be imported in the description file of the package.
      # indexOfChunk<-1 #To comment out.

      cat("\nRunning Chunk",indexOfChunk,"out of",nrow(chunkInfo),"chunks...\n") #Include a new line at the beginning of this print message to emphasize it.

      dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.

      runChunk(dataGeneExpressionFPSub,dataCovariates,
               genotypeFile,tabixProgram,sampleIndices,
               approach,B,cisDistance,MAFThreshold,MASamplesThreshold,
               indexOfChunk,outputDir)
      gc()
      return(0)
    },mc.cores=numOfCores) #results is a list of length nrow(chunkInfo). Each core takes up to a few GB of memory.

    # chunkInfo<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/Lung/_chunkInfo.rds")
    # results<-rep(0,nrow(chunkInfo))
    # results<-as.list(results)
    cat("\n") #Include a new line here to emphasize the print message.
    print(unlist(results)) #0 means success.
  }

  sink() #Close the connection.

}
