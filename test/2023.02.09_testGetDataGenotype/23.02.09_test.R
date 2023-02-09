



library(dplyr)
Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.02.09_testGetDataGenotype/3.3_getDataGenotypeCpp_V1.cpp")

#Prepare data. See runChunk.R.
if(TRUE){
  #Load dataGeneExpressionFPSub and dataCovariates.
  if(TRUE){
    tissueType<-"Lung" #Sample size is 515.
    numOfPCs<-44 #Chosen via BE.

    # tissueType<-"Brain_Cortex" #Sample size is 205.
    # numOfPCs<-23 #Chosen via BE.

    exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
                     tissueType,".v8.normalized_expression.bed.gz")
    covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
                    tissueType,"/",numOfPCs,"ExprPCs.txt")

    source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.1_prepareExprAndCovData.R")
    temp<-prepareExprAndCovData(exprFile,covFile)
    dataGeneExpressionFP<-temp$dataGeneExpressionFP #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
    dataCovariates<-temp$dataCovariates #515*52. The columns are the known and inferred covariates. No need to convert the matrix to a data frame.

    indexOfChunk<-5
    outputDir<-"~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/"
    chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
    source("~/2022.03.14_ClipperQTL/ClipperQTL/R/3.1_prepareDataGeneExpressionFPSub.R")
    dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.

    rm(list=setdiff(ls(),c("getDataGenotypeCpp","residualizeCpp","getTableMaxAbsCorsCpp",
                           "dataGeneExpressionFPSub","dataCovariates","indexOfChunk","outputDir")))
  }

  genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
  sampleIndices<-prepareSampleIndices(genotypeFile,tabixProgram,outputDir,
                                      sampleNames=colnames(dataGeneExpressionFPSub)[-(1:4)]) #Vector of length 515.
  rm(prepareSampleIndices)

  B<-20
  MAFThreshold<-0.01
  MASamplesThreshold<-10

  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  #To comment out.
}

geneTSSs<-dataGeneExpressionFPSub$end #Vector of length 257. Key variable.
Y<-as.matrix(dataGeneExpressionFPSub[,-(1:4)]) #257*515. Key variable.
# dim(Y)

#Get dataGenotype (SNPInfo and X). This takes about 11 seconds.
if(TRUE){
  cat("\nObtaining genotype data for Chunk ",indexOfChunk,"...\n",sep="")
  timeStart<-Sys.time()

  dataGenotype<-getDataGenotypeCpp(chrCurr=dataGeneExpressionFPSub$chr[1],geneTSSs,genotypeFile,tabixProgram,tempFileName=tempfile(),
                                   sampleIndices,MAFThreshold,MASamplesThreshold) #A list of two items: SNPInfo and X.

  SNPInfo<-dataGenotype$SNPInfo #131,682*3. The columns are: CHROM, POS, ID. 131,511. 136,253. 111,810.
  SNPPositions<-as.integer(SNPInfo[,2]) #Vector of length 131,682. Key variable.
  X<-dataGenotype$X #131,682*515. Key variable.

  cat(dim(SNPInfo))

  timeEnd<-Sys.time()
  print(timeEnd-timeStart) #print() is better than cat() for time difference.
}












