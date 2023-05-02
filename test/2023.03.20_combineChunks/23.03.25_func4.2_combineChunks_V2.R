


combineChunks_V2<-function(outputDir,chunkInfo){
  # tissueType<-"Lung" #Sample size is 515.
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.03.20_combineChunks/_result/",tissueType,"/")
  # chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  # #To comment out.

  chunkInfo<-chunkInfo%>%mutate(numOfGenesInChunk=geneIndexWithinChrEnd-geneIndexWithinChrStart+1)
  numOfGenesTotal<-sum(chunkInfo$numOfGenesInChunk) #26,095.

  chunkInfo<-chunkInfo%>%mutate(indexOfRowEnd=cumsum(numOfGenesInChunk))
  chunkInfo<-chunkInfo%>%mutate(indexOfRowStart=indexOfRowEnd-numOfGenesInChunk+1)
  chunkInfo<-chunkInfo[,c(1:5,7,6)] #Switch indexOfRowStart and indexOfRowEnd. This is unnecessary but may be helpful during code development.
  # identical(chunkInfo$indexOfRowEnd-chunkInfo$indexOfRowStart+1,chunkInfo$numOfGenesInChunk) #TRUE is good.

  resultChunk1<-readRDS(paste0(outputDir,"resultChunk1.rds"))
  resultCombined<-data.frame(matrix(nrow=numOfGenesTotal,ncol=ncol(resultChunk1)))
  colnames(resultCombined)<-colnames(resultChunk1)

  for(indexOfChunk in chunkInfo$indexOfChunk){
    # indexOfChunk<-1 #To comment out.

    pathResultChunk<-paste0(outputDir,"resultChunk",indexOfChunk,".rds") #resultChunk1.rds.
    resultChunk<-readRDS(pathResultChunk)

    indexOfRowStartCurr<-chunkInfo$indexOfRowStart[indexOfChunk]
    indexOfRowEndCurr<-chunkInfo$indexOfRowEnd[indexOfChunk]

    resultCombined[indexOfRowStartCurr:indexOfRowEndCurr,]<-resultChunk
  }

  pathResultCombined<-paste0(outputDir,"_resultCombined.rds")
  saveRDS(resultCombined,pathResultCombined)
}













