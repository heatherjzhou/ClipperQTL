#' ClipperQTL
#'
#' @export


# #For code development only:
# library(dplyr)
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.3_prepareChunkInfo.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.4_prepareMethodParameters.R")
# source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.5_prepareDataGeneExpressionFPSub.R")



ClipperQTL<-function(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
                     approach="standard",B=1000, #approach is "standard" or "Clipper".
                     cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
                     numOfChunksTarget=100,seed=1,numOfCores=1){
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
  # approach<-"standard"
  # B<-1000
  #
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # numOfChunksTarget<-100 #The total number of chunks will be approximately this number.
  # seed<-1
  # numOfCores<-1
  # #To comment out.



  #Create output directory.
  if(!dir.exists(outputDir)) dir.create(outputDir)

  filename<-paste0(outputDir,"_log_ClipperQTL.txt")
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

  #Prepare chunkInfo.
  cat("Partitioning genes into approximately",numOfChunksTarget,"chunks...\n")
  chunkInfo<-prepareChunkInfo(dataGeneExpressionFP,numOfChunksTarget,outputDir) #103*4.

  #Prepare method parameters.
  temp<-prepareMethodParameters(approach,B,sampleSize=nrow(dataCovariates))
  approach<-temp$approach #"standard" or "Clipper".
  B<-temp$B
  rm(temp)

  #Run chunks. This creates resultChunk1.rds, resultChunk2.rds, etc.
  if(TRUE){
    RNGkind("L'Ecuyer-CMRG") #This is necessary to ensure reproducibility when using mclapply().
    set.seed(seed) #For permutations.
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

  #Save _resultCombined.rds (to be deleted in the next step).
  cat("\nCombining results...\n") #Include a new line at the beginning of this print message to emphasize it.
  resultCombined<-combineChunks(outputDir,chunkInfo)

  #Save _resultGenes.rds.
  cat("\nCalling eGenes...\n") #Include a new line at the beginning of this print message to emphasize it.
  resultGenes<-callEGenes(resultCombined,approach,B,
                          outputDir)

  cat("\nClipperQTL() finished running.\n") #Include a new line at the beginning of this print message to emphasize it.

  sink() #Close the connection.
}










