



# #For code development only:
# library(dplyr)
# Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/src/3.2_getDataGenotypeCpp.cpp")
# Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/src/3.3_getTableMaxAbsCorsCpp.cpp")



runChunk<-function(dataGeneExpressionFPSub,dataCovariates,
                   genotypeFile,tabixProgram,sampleIndices,
                   approach,B,cisDistance,MAFThreshold,MASamplesThreshold,
                   indexOfChunk,outputDir){
  # #Load dataGeneExpressionFPSub and dataCovariates.
  # if(TRUE){
  #   tissueType<-"Lung" #Sample size is 515.
  #   numOfPCs<-44 #Chosen via BE.
  #
  #   exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                    tissueType,".v8.normalized_expression.bed.gz")
  #   covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
  #                   tissueType,"/",numOfPCs,"ExprPCs.txt")
  #
  #   source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")
  #   temp<-prepareExprAndCovData(exprFile,covFile)
  #   dataGeneExpressionFP<-temp$dataGeneExpressionFP #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #   dataCovariates<-temp$dataCovariates #515*52. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.
  #
  #   #Load dataGeneExpressionFPSub.
  #   indexOfChunk<-5
  #   outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  #   chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  #   source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.5_prepareDataGeneExpressionFPSub.R")
  #   dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  #   rm(list=setdiff(ls(),c("getDataGenotypeCpp","getTableMaxAbsCorsCpp",
  #                          "dataGeneExpressionFPSub","dataCovariates","indexOfChunk","outputDir")))
  # }
  #
  # genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
  # sampleIndices<-prepareSampleIndices(genotypeFile,tabixProgram,outputDir,
  #                                     sampleNames=colnames(dataGeneExpressionFPSub)[-(1:4)]) #Vector of length 515.
  # rm(prepareSampleIndices)
  #
  # approach<-"Clipper" #Use this to make testing faster.
  # B<-20 #Use this to make testing faster.
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # RNGkind("L'Ecuyer-CMRG")
  # set.seed(1)
  # #To comment out.



  # #Load dataGeneExpressionFPSub and dataCovariates.
  # if(TRUE){
  #   cellTypeIndex<-0 #Sample size is 981.
  #   numOfPCs<-19 #Chosen via elbow.
  #
  #   exprFile<-paste0("~/_Data/2023.09.11_OneK1K/_2023.09.18_De_Novo/1.2_expr_fullyProcessed/withoutINT/cellType",cellTypeIndex,".bed.gz")
  #   covFile<-paste0("~/_Data/2023.09.11_OneK1K/_2023.09.18_De_Novo/3.4_covariatesCombined/cellType",cellTypeIndex,"/withoutINT_",numOfPCs,"ExprPCs.txt")
  #
  #   source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")
  #   temp<-prepareExprAndCovData(exprFile,covFile)
  #   dataGeneExpressionFP<-temp$dataGeneExpressionFP #9290*985. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #   dataCovariates<-temp$dataCovariates #981*21. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.
  #
  #   #Load dataGeneExpressionFPSub.
  #   indexOfChunk<-5
  #   outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/cellType",cellTypeIndex,"/")
  #   chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  #   source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.5_prepareDataGeneExpressionFPSub.R")
  #   dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  #   rm(list=setdiff(ls(),c("getDataGenotypeCpp","getTableMaxAbsCorsCpp",
  #                          "dataGeneExpressionFPSub","dataCovariates","indexOfChunk","outputDir")))
  # }
  #
  # genotypeFile<-"~/_Data/2023.09.11_OneK1K/_2023.09.12_Authors/filter_vcf_r08/_combined.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
  # sampleIndices<-prepareSampleIndices(genotypeFile,tabixProgram,outputDir,
  #                                     sampleNames=colnames(dataGeneExpressionFPSub)[-(1:4)]) #Vector of length 515.
  # rm(prepareSampleIndices)
  #
  # approach<-"Clipper" #Use this to make testing faster.
  # B<-20 #Use this to make testing faster.
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # RNGkind("L'Ecuyer-CMRG")
  # set.seed(1)
  # #To comment out.



  geneTSSs<-dataGeneExpressionFPSub$end #Vector of length 257. Key variable.
  Y<-as.matrix(dataGeneExpressionFPSub[,-(1:4)]) #257*515. Key variable.
  # dim(Y)

  #Get dataGenotype (SNPInfo and X). This takes about 11 seconds.
  if(TRUE){
    timeStart<-Sys.time()

    cat("Obtaining genotype data for Chunk ",indexOfChunk,"...\n",sep="")

    dataGenotype<-getDataGenotypeCpp(chrCurr=as.character(dataGeneExpressionFPSub$chr[1]),geneTSSs,genotypeFile,tabixProgram,tempFileName=tempfile(),
                                     sampleIndices,cisDistance,MAFThreshold,MASamplesThreshold) #A list of two items: SNPInfo and X.

    SNPInfo<-dataGenotype$SNPInfo #131,682*3. The columns are: CHROM, POS, ID.
    SNPPositions<-as.integer(SNPInfo[,2]) #Vector of length 131,682. Key variable.
    X<-dataGenotype$X #131,682*515. Key variable.

    timeEnd<-Sys.time()
    cat("Obtaining genotype data for Chunk ",indexOfChunk," took: ",sep="")
    print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }

  # dim(SNPInfo) #131,682*3.
  # dim(X) #131,682*515.
  # sum(abs(SNPPositions-geneTSSs[1])<=1e6) #7261 local common SNPs for the first gene in Chunk 5 of Lung. Matches FastQTL.
  # sum(abs(SNPPositions-geneTSSs[2])<=1e6) #7278 local common SNPs for the second gene in Chunk 5 of Lung. Matches FastQTL.

  # dim(SNPInfo) #101,908*3.
  # dim(X) #101,908*981.
  # sum(abs(SNPPositions-geneTSSs[1])<=1e6) #3514 local common SNPs for the first gene in Chunk 5 of cell type 0. Matches FastQTL.
  # sum(abs(SNPPositions-geneTSSs[2])<=1e6) #3678 local common SNPs for the second gene in Chunk 5 of cell type 0. Matches FastQTL.

  #Get tableMaxAbsCors. This takes about 16 or 24 seconds, depending on conditions. Takes less time as of 2025/01/22.
  if(TRUE){
    timeStart<-Sys.time()

    cat("Calculating maximum absolute correlations for Chunk ",indexOfChunk,"...\n",sep="")

    tableMaxAbsCors<-getTableMaxAbsCorsCpp(Y, #257*515.
                                           X, #131,682*515.
                                           dataCovariates, #515*52.
                                           geneTSSs, #Vector of length 257.
                                           SNPPositions, #Vector of length 131,682.
                                           cisDistance,
                                           approach,
                                           B
                                           ) #257*21. Each row corresponds to a gene. The columns are: exp, bg1, ..., bg20.

    timeEnd<-Sys.time()
    cat("Calculating maximum absolute correlations for Chunk ",indexOfChunk," took: ",sep="")
    print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }

  #Create resultChunk.
  colnames(tableMaxAbsCors)<-c("exp",paste0("bg",c(1:B)))
  resultChunk<-cbind(dataGeneExpressionFPSub[,1:4],tableMaxAbsCors)

  #Save resultChunk.
  path<-paste0(outputDir,"resultChunk",indexOfChunk,".rds") #resultChunk5.rds.
  saveRDS(resultChunk,path)

  # #Compare resultChunk to resultChunk_freeze.
  # resultChunk_freeze<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/Lung/resultChunk5_freeze.rds")
  # identical(resultChunk_freeze[,1:4],resultChunk[,1:4]) #TRUE is good.
  # max(abs(resultChunk_freeze$exp-resultChunk$exp)) #9.992007e-16.
  # max(abs(resultChunk_freeze[,-(1:4)]-resultChunk[,-(1:4)])) #2.664535e-15.
}









