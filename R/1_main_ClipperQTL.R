#' @title
#' ClipperQTL
#'
#' @description
#' This function is used for identifying eGenes.
#'
#' @details
#' \code{ClipperQTL()} requires three main pieces of input data: expression data, covariate data, and genotype data (the file paths are specified by \code{exprFile}, \code{covFile}, and \code{genotypeFile}, respectively). All three data sets must have the same format as the data sets used in GTEx V8 analysis.
#'
#' Specifically, the expression file must be a .bed.gz file. The first four columns must be chr, start, end, and gene_id (the exact column names do not matter). Each remaining column must correspond to a sample. The third column, end, will be used as the transcription start site (following FastQTL; so make sure to put the desired transcription start site in this column). The second column, start, will not be used.
#'
#' The covariate file must be a .txt file. The first column must be the covariate ID (the exact column name does not matter). Each remaining column must correspond to a sample. All constant covariates will be filtered out before analysis.
#'
#' The genotype file must be a .vcf.gz file. A .tbi index file must be present in the same directory, which can be generated using HTSlib tabix. The non-empty genotype entries must start with "0|0", "0|1", "1|0", or "1|1" (trailing characters are ok). The missing genotype entries will be imputed as within-SNP averages. This file can contain more samples than the expression file and the covariate file (which must have the same samples).
#'
#' The main method parameters of \code{ClipperQTL()} are \code{approach} and \code{B}. If the sample size is less than or equal to 450, then we recommend setting \code{approach="standard"} and \code{B=1000}. If the sample size is greater than 450, then you may set \code{approach="standard"} and \code{B=1000}, or, you may set \code{approach="Clipper"} and \code{B=1} or \code{B} between 20 and 100 for faster computational speed.
#'
#' \code{ClipperQTL()} outputs several files in the output directory. The most important one is named "_resultGenes.rds", which can be read into R with \code{readRDS()}. The first four columns are identical to the first four columns in the expression file. The next \code{1+B} columns are the maximum absolute correlations from the experimental round and the permutation rounds. The second to last column is \code{pValue} (standard variant) or \code{contrastScore} (Clipper variant). The last column is \code{qValue}. Each row corresponds to a gene. The eGenes are those with \code{qValue} under the target FDR threshold, e.g., 0.05.
#'
#' @param exprFile The directory and filename of the expression file (.bed.gz file; see below).
#' @param covFile The directory and filename of the covariate file (.txt file; see below).
#' @param genotypeFile The directory and filename of the genotype file (.vcf.gz file; see below).
#' @param tabixProgram The directory and filename of the tabix executable file.
#' @param outputDir The output directory (ending with "/").
#'
#' @param approach Key argument. Must be \code{"standard"} or \code{"Clipper"}. See below.
#' @param B Key argument. The number of permutations. See below.
#'
#' @param cisDistance The maximum distance between a SNP and the transcription start site of a gene for the SNP to be considered local for the gene. Default is \code{1e6}.
#' @param MAFThreshold The threshold for the minor allele frequency (MAF) of a SNP. Default is \code{0.01}.
#' @param MASamplesThreshold The threshold for the number of samples with at least one copy of the minor allele (MA samples). Default is \code{10}.
#'
#' @param numOfChunksTarget The target number chunks that all genes in the expression file will be divided into. Default is \code{100}.
#' @param seed The seed. Default is \code{1}.
#' @param numOfCores The number of cores to be used. Default is \code{1}.
#'
#' @references
#' Heather J. Zhou, Xinzhou Ge, and Jingyi Jessica Li. ClipperQTL: ultrafast and powerful eGene identification method. Genome Biology, 26:207, 2025.
#'
#' GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science, 369(6509):1318–1330, 2020.
#'
#' Halit Ongen, Alfonso Buil, Andrew Anand Brown, Emmanouil T. Dermitzakis, and Olivier Delaneau. Fast and efficient QTL mapper for thousands of molecular phenotypes. Bioinformatics, 32(10):1479–1485, 2016.
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
                     approach,B, #approach must be "standard" or "Clipper".
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



  # cellTypeIndex<-0 #Sample size is 981.
  # numOfPCs<-19 #Chosen via elbow.
  #
  # exprFile<-paste0("~/_Data/2023.09.11_OneK1K/_2023.09.18_De_Novo/1.2_expr_fullyProcessed/withoutINT/cellType",cellTypeIndex,".bed.gz")
  # covFile<-paste0("~/_Data/2023.09.11_OneK1K/_2023.09.18_De_Novo/3.4_covariatesCombined/cellType",cellTypeIndex,"/withoutINT_",numOfPCs,"ExprPCs.txt")
  # genotypeFile<-"~/_Data/2023.09.11_OneK1K/_2023.09.12_Authors/filter_vcf_r08/_combined.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/cellType",cellTypeIndex,"/")
  # rm(cellTypeIndex,numOfPCs)
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




