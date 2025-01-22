


prepareChunkInfo<-function(dataGeneExpressionFP,numOfChunksTarget,outputDir){
  # #Load dataGeneExpressionFP.
  # tissueType<-"Lung" #Sample size is 515.
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # dataGeneExpressionFP<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  # colnames(dataGeneExpressionFP)[1:4]<-c("chr","start","end","gene_id") #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  # numOfChunksTarget<-100 #The total number of chunks will be approximately this number.
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  # #To comment out.



  # cellTypeIndex<-0 #Sample size is 981.
  # exprFile<-paste0("~/_Data/2023.09.11_OneK1K/_2023.09.18_De_Novo/1.2_expr_fullyProcessed/withoutINT/cellType",cellTypeIndex,".bed.gz")
  # dataGeneExpressionFP<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  # colnames(dataGeneExpressionFP)[1:4]<-c("chr","start","end","gene_id") #9290*985. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  # numOfChunksTarget<-100 #The total number of chunks will be approximately this number.
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/cellType",cellTypeIndex,"/")
  # #To comment out.

  #Create dataChrs.
  dataChrs<-dataGeneExpressionFP%>%group_by(chr)%>%
    summarize(numOfGenes=n())
  if(substring(dataChrs$chr[1],1,3)=="chr"){
    dataChrs<-dataChrs%>%mutate(chrNumber=substring(chr,4,length(chr))) #If chromosome names start with "chr", then remove "chr" to create chrNumber.
  }else{
    dataChrs<-dataChrs%>%mutate(chrNumber=chr)
  }
  suppressWarnings(dataChrs$chrNumber<-as.numeric(dataChrs$chrNumber)) #Possible warning message: NAs introduced by coercion.
  dataChrs<-dataChrs%>%arrange(chrNumber,chr) #Arrange the chromosomes by chromosome number.
  dataChrs<-dataChrs%>%select(-chrNumber) #23*2.

  #Calculate numOfChunks and numOfGenesPerChunkApprox.
  dataChrs<-dataChrs%>%mutate(numOfChunks=round(numOfGenes/nrow(dataGeneExpressionFP)*numOfChunksTarget))
  dataChrs<-dataChrs%>%mutate(numOfGenesPerChunkApprox=round(numOfGenes/numOfChunks))
  # sum(dataChrs$numOfChunks) #103.

  #Create chunkInfo.
  chunkInfo<-data.frame(matrix(ncol=7,nrow=0))
  for(i in 1:nrow(dataChrs)){
    # i<-3

    dataChr<-dataChrs[i,]
    numOfChunksChr<-dataChr$numOfChunks #6.
    dataChr<-dataChr[rep(1,numOfChunksChr),]
    dataChr$indexOfChunkWithinChr<-1:nrow(dataChr)
    dataChr<-dataChr%>%mutate(geneIndexWithinChrStart=(indexOfChunkWithinChr-1)*numOfGenesPerChunkApprox+1)
    dataChr<-dataChr%>%mutate(geneIndexWithinChrEnd=indexOfChunkWithinChr*numOfGenesPerChunkApprox)
    dataChr$geneIndexWithinChrEnd[nrow(dataChr)]<-dataChr$numOfGenes[1] #Change the last geneIndexWithinChrEnd to the number of genes in the chromosome.

    chunkInfo<-rbind(chunkInfo,dataChr)
  }
  chunkInfo$indexOfChunk<-1:nrow(chunkInfo)
  # chunkInfo<-chunkInfo%>%select(indexOfChunk,everything())
  chunkInfo<-chunkInfo%>%select(indexOfChunk,chr,
                                geneIndexWithinChrStart,geneIndexWithinChrEnd)

  #Edge case (OneK1K data, cell type 0, 450 genes total).
  chunkInfo<-chunkInfo%>%filter(geneIndexWithinChrStart<=geneIndexWithinChrEnd)
  chunkInfo$indexOfChunk<-1:nrow(chunkInfo)

  # #Check chunkInfo.
  # chunkInfo<-chunkInfo%>%mutate(numOfGenesInChunk=geneIndexWithinChrEnd-geneIndexWithinChrStart+1)
  # sum(chunkInfo$numOfGenesInChunk) #26,095.

  saveRDS(chunkInfo,paste0(outputDir,"_chunkInfo.rds"))

  return(chunkInfo)
}









