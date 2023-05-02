


combineChunks_V1<-function(outputDir,chunkInfo){
  # tissueType<-"Lung" #Sample size is 515.
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  # chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.

  pathResultCombined<-paste0(outputDir,"_resultCombined.txt")
  for(indexOfChunk in chunkInfo$indexOfChunk){
    # indexOfChunk<-1 #To comment out.

    pathResultChunk<-paste0(outputDir,"resultChunk",indexOfChunk,".rds") #resultChunk1.rds.
    resultChunk<-readRDS(pathResultChunk)
    if(indexOfChunk==1){
      write.table(resultChunk,pathResultCombined,
                  append=FALSE,quote=FALSE,
                  sep="\t",row.names=FALSE,col.names=TRUE)
    }else{
      write.table(resultChunk,pathResultCombined,
                  append=TRUE,quote=FALSE,
                  sep="\t",row.names=FALSE,col.names=FALSE)
    }

    # unlink(pathResultChunk)
  }

  # resultCombined<-read.delim(pathResultCombined) #Import from text (base) works.
  # return(resultCombined) #26,095*25. The first four columns are chr, start, end, and gene_id. end is used as TSS.
}













