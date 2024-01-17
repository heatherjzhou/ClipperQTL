


# #For code development only:
# library(dplyr)
# Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/src/3.2_getDataGenotypeCpp.cpp")
# Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/src/4.3_saveTemporaryResultFilesCpp.cpp")



runChunk_sigGeneSNPPairs<-function(dataGeneExpressionFPSub,dataCovariates,
                                   genotypeFile,tabixProgram,sampleIndices,
                                   cisDistance,MAFThreshold,MASamplesThreshold,
                                   indexOfChunk,outputDir,
                                   sigGeneSNPPairMethod,absCorThresholds,percent){ #absCorThresholds is only relevant when using "FastQTL". percent is only relevant when using "topPercent".
  # #Load dataGeneExpressionFPSub, dataCovariates, and absCorThresholds.
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
  #   outputDir<-paste0("/home/heatherjzhou/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/") #Key.
  #   chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  #   source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.5_prepareDataGeneExpressionFPSub.R")
  #   dataGeneExpressionFPSub<-prepareDataGeneExpressionFPSub(dataGeneExpressionFP,indexOfChunk,chunkInfo) #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  #   #Filter dataGeneExpressionFPSub to only contain eGenes.
  #   FDR_eGene<-0.05
  #   resultGenes<-readRDS(paste0(outputDir,"_resultGenes.rds"))
  #   eGenes<-(resultGenes%>%filter(qValue<=FDR_eGene))$gene_id #Vector of length 13,629.
  #   dataGeneExpressionFPSub<-dataGeneExpressionFPSub%>%filter(gene_id%in%eGenes) #127*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  #   sigGeneSNPPairMethod<-"FastQTL" #Key.
  #   if(sigGeneSNPPairMethod=="FastQTL"){
  #     #Calculate pt.
  #     resultGenes_yes<-resultGenes%>%filter(gene_id%in%eGenes) #13,629*1007.
  #     resultGenes_no<-resultGenes%>%filter(!gene_id%in%eGenes) #12,466*1007.
  #     p1<-max(resultGenes_yes$pValue)
  #     p2<-min(resultGenes_no$pValue)
  #     pt<-(p1+p2)/2 #0.08641359.
  #
  #     #Calculate absCorThreshold (a variable in resultGenes_yes), the threshold for calling significant gene-SNP pairs.
  #     #The p-value of a gene-SNP pair being in the top X (pt*100) percent is equivalent to the absolute correlation of the gene-SNP pair being in the top X percent, where "top" means the most significant.
  #     maxAbsCors<-as.matrix(resultGenes_yes[,-c(1:5,(ncol(resultGenes_yes)-1):ncol(resultGenes_yes))]) #13,629*1000. The columns are: bg1, ..., bg1000.
  #     # dim(maxAbsCors)
  #     resultGenes_yes$absCorThreshold<-apply(maxAbsCors,1,function(x){ #1 means by row. This takes a moment.
  #       return(quantile(x,probs=1-pt)) #Use 1-pt because larger absolute correlations are more significant.
  #     }) #13,629*1008. absCorThreshold is 0.1724693, 0.1767774, etc.
  #
  #     absCorThresholds<-(resultGenes_yes%>%filter(gene_id%in%dataGeneExpressionFPSub$gene_id))$absCorThreshold #Vector of length 127. 0.1838889, 0.1879149, etc.
  #   }else{
  #     absCorThresholds<-0 #Use 0 as a placeholder.
  #   }
  #
  #   rm(list=setdiff(ls(),c("getDataGenotypeCpp","saveTemporaryResultFilesCpp",
  #                          "dataGeneExpressionFPSub","dataCovariates",
  #                          "indexOfChunk","outputDir",
  #                          "sigGeneSNPPairMethod","absCorThresholds")))
  # }
  #
  # genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # source("~/2022.03.14_ClipperQTL/ClipperQTL/R/2.2_prepareSampleIndices.R")
  # sampleIndices<-prepareSampleIndices(genotypeFile,tabixProgram,outputDir,
  #                                     sampleNames=colnames(dataGeneExpressionFPSub)[-(1:4)]) #Vector of length 515.
  # rm(prepareSampleIndices)
  #
  # cisDistance<-1e6
  # MAFThreshold<-0.01
  # MASamplesThreshold<-10
  #
  # percent<-1
  # #To comment out.

  geneTSSs<-dataGeneExpressionFPSub$end #Vector of length 127. Key variable.
  Y<-as.matrix(dataGeneExpressionFPSub[,-(1:4)]) #127*515. Key variable.
  # dim(Y)

  #Get dataGenotype (SNPInfo and X). This takes about 11 seconds.
  if(TRUE){
    timeStart<-Sys.time()

    cat("Obtaining genotype data for Chunk ",indexOfChunk,"...\n",sep="")

    dataGenotype<-getDataGenotypeCpp(chrCurr=as.character(dataGeneExpressionFPSub$chr[1]),geneTSSs,genotypeFile,tabixProgram,tempFileName=tempfile(),
                                     sampleIndices,cisDistance,MAFThreshold,MASamplesThreshold) #A list of two items: SNPInfo and X.

    SNPInfo<-dataGenotype$SNPInfo #131,649*3. The columns are: CHROM, POS, ID.
    SNPPositions<-as.integer(SNPInfo[,2]) #Vector of length 131,649. Key variable.
    X<-dataGenotype$X #131,649*515. Key variable.

    timeEnd<-Sys.time()
    cat("Obtaining genotype data for Chunk ",indexOfChunk," took: ",sep="")
    print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }

  # dim(SNPInfo) #131,649*3.
  # dim(X) #131,649*515.
  # sum(abs(SNPPositions-geneTSSs[1])<=1e6) #7278.
  # sum(abs(SNPPositions-geneTSSs[2])<=1e6) #7336.

  #Save temporary result files, one for each gene. This takes about 14 seconds.
  if(TRUE){
    timeStart<-Sys.time()

    cat("Creating temporary result files for Chunk ",indexOfChunk,"...\n",sep="")

    outputDirChunk<-paste0(outputDir,"temp_Chunk",indexOfChunk,"/")
    if(!dir.exists(outputDirChunk)) dir.create(outputDirChunk)

    success<-saveTemporaryResultFilesCpp(Y, #127*515.
                                         X, #131,649*515.
                                         dataCovariates, #515*52.
                                         geneTSSs, #Vector of length 127.
                                         SNPPositions, #Vector of length 131,649.
                                         cisDistance,
                                         sigGeneSNPPairMethod,
                                         absCorThresholds,
                                         percent,
                                         indexOfChunk,
                                         outputDirChunk
    ) #0 means success.

    timeEnd<-Sys.time()
    cat("Creating temporary result files for Chunk ",indexOfChunk," took: ",sep="")
    print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }

  # temp<-read_csv("ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/temp_Chunk5/Chunk5_Gene1.csv",col_names=FALSE) #883*2.
  # dim(temp)

  #Combine temporary result files into resultChunk5.txt, for example. This takes about 6 seconds.
  if(TRUE){
    timeStart<-Sys.time()

    cat("Combining temporary result files for Chunk ",indexOfChunk,"...\n",sep="")

    pathResultChunk<-paste0(outputDir,"resultChunk",indexOfChunk,".txt")
    for(indexOfGene in 1:nrow(Y)){
      # indexOfGene<-1

      resultGene<-readr::read_csv(paste0(outputDirChunk,"Chunk",indexOfChunk,"_Gene",indexOfGene,".csv"),col_names=FALSE,show_col_types=FALSE) #833*2.
      colnames(resultGene)<-c("indexOfSNP","partialCor")
      indicesOfSNPs<-resultGene$indexOfSNP
      # min(abs(resultGene$partialCor)) #0.1840818, slightly larger than absCorThresholds[1], which is 0.1838889.

      chr<-dataGeneExpressionFPSub$chr[indexOfGene]
      geneID<-dataGeneExpressionFPSub$gene_id[indexOfGene]
      toSave<-cbind(chr,geneID, #Chromosome, gene ID.
                    (SNPInfo[,2])[indicesOfSNPs],(SNPInfo[,3])[indicesOfSNPs], #SNP position, SNP ID.
                    resultGene$partialCor #Partial correlation.
      )
      # dim(toSave) #833*5.
      colnames(toSave)<-c("chr","gene_id",
                          "SNP_position","SNP_id",
                          "partialCor")

      if(indexOfGene==1){
        write.table(toSave,pathResultChunk,
                    append=FALSE,quote=FALSE,
                    sep="\t",row.names=FALSE,col.names=TRUE)
      }else{
        write.table(toSave,pathResultChunk,
                    append=TRUE,quote=FALSE,
                    sep="\t",row.names=FALSE,col.names=FALSE)
      }
    }

    timeEnd<-Sys.time()
    cat("Combining temporary result files for Chunk ",indexOfChunk," took: ",sep="")
    print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }

  # temp<-read.delim("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/resultChunk5.txt") #14,571*5. Import from text (base) works.
  # dim(temp)

  #Delete the entire folder of outputDirChunk.
  unlink(outputDirChunk,recursive=TRUE)
}








