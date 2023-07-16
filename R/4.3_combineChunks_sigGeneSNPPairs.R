

combineChunks_sigGeneSNPPairs<-function(outputDir,chunkInfo,n,K,sigGeneSNPPairMethod){
  # outputDir<-paste0("/home/heatherjzhou/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/")
  # chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  # n<-515 #Sample size.
  # K<-52 #Number of covariates (not including the intercept).
  # sigGeneSNPPairMethod<-"FastQTL"
  # #To comment out.

  #Calculate numOfRowsTotal. This takes about 6 seconds.
  if(TRUE){
    # timeStart<-Sys.time()

    # resultChunks<-vector(mode="list",length=nrow(chunkInfo))
    numOfRowsTotal<-0
    for(indexOfChunk in chunkInfo$indexOfChunk){
      # indexOfChunk<-1 #To comment out.

      pathResultChunk<-paste0(outputDir,"resultChunk",indexOfChunk,".txt") #resultChunk1.txt.

      resultChunk<-read.delim(pathResultChunk) #10,674*5. Import from text (base) works.

      # resultChunks[[indexOfChunk]]<-resultChunk
      numOfRowsTotal<-numOfRowsTotal+nrow(resultChunk)
    }

    # timeEnd<-Sys.time()
    # print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }


  #Create resultCombined (to be filled).
  resultChunk1<-read.delim(paste0(outputDir,"resultChunk1.txt")) #10,674*5. Import from text (base) works.
  resultCombined<-data.frame(matrix(nrow=numOfRowsTotal,ncol=ncol(resultChunk1))) #1,948,850*5.
  colnames(resultCombined)<-colnames(resultChunk1)

  #Fill resultCombined. This takes about 34 seconds.
  if(TRUE){
    # timeStart<-Sys.time()

    numOfRowsFilled<-0
    for(indexOfChunk in chunkInfo$indexOfChunk){
      # indexOfChunk<-1 #To comment out.

      pathResultChunk<-paste0(outputDir,"resultChunk",indexOfChunk,".txt") #resultChunk1.txt.
      resultChunk<-read.delim(pathResultChunk) #10,674*5. Import from text (base) works.

      resultCombined[(numOfRowsFilled+1):(numOfRowsFilled+nrow(resultChunk)),]<-resultChunk

      numOfRowsFilled<-numOfRowsFilled+nrow(resultChunk)

      unlink(pathResultChunk)
    }

    # timeEnd<-Sys.time()
    # print(timeEnd-timeStart) #print() is better than cat() for time difference.
  }


  #Calculate tStat and pVal.
  partialCors<-resultCombined$partialCor
  resultCombined$tStat<-partialCors*sqrt((n-2-K)/(1-partialCors^2)) #tStat.
  resultCombined$pVal<-2*pt(q=abs(resultCombined$tStat),df=n-2-K,lower.tail=FALSE,log.p=FALSE) #pVal.

  #Save _resultSigGeneSNPPairs.rds.
  pathResultCombined<-paste0(outputDir,"_resultSigGeneSNPPairs_",sigGeneSNPPairMethod,".rds")
  saveRDS(resultCombined,pathResultCombined)

  # return(resultCombined) #1,948,850*7. The columns are: chr, gene_id, SNP_position, SNP_id, paritalCor, tStat, pVal.
}
