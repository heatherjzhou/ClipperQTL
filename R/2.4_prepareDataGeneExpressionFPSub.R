


prepareDataGeneExpressionFPSub<-function(dataGeneExpressionFP,indexOfChunk,chunkInfo){
  # #Load dataGeneExpressionFP.
  # tissueType<-"Lung" #Sample size is 515.
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # dataGeneExpressionFP<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  # colnames(dataGeneExpressionFP)[1]<-"chr" #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  #
  # indexOfChunk<-5
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  # chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.
  # #To comment out.

  chrCurr<-chunkInfo$chr[indexOfChunk] #chr1.
  geneIndexWithinChrStart<-chunkInfo$geneIndexWithinChrStart[indexOfChunk] #1029.
  geneIndexWithinChrEnd<-chunkInfo$geneIndexWithinChrEnd[indexOfChunk] #1285.

  dataGeneExpressionFPSub<-dataGeneExpressionFP%>%filter(chr==chrCurr) #2571*519.
  if(is.unsorted(dataGeneExpressionFPSub$end)){ #If the genes are not already sorted by TSS, sort them.
    dataGeneExpressionFPSub<-dataGeneExpressionFPSub%>%arrange(end)
  }
  dataGeneExpressionFPSub<-dataGeneExpressionFPSub[geneIndexWithinChrStart:geneIndexWithinChrEnd,] #257*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.

  return(dataGeneExpressionFPSub)
}







