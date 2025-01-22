


callEGenes<-function(resultCombined,approach,B,
                     outputDir){

  if(approach=="standard"){
    # resultCombined<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/_resultGenes.rds")
    # approach<-"standard"
    # B<-1000
    # outputDir<-"~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=standard_B=1000/"
    # #To comment out.

    expScores<-resultCombined[,5] #Vector of length 26,095.
    bgScores<-resultCombined[,6:(4+1+B)]
    # dim(bgScores) #26,095*1000.

    temp<-bgScores-expScores
    # dim(temp) #26,095*1000.
    numsOfPermsMoreSig<-rowSums(bgScores-expScores>=0)

    resultCombined$pValue<-(numsOfPermsMoreSig+1)/(B+1) #pValue is calculated as (X+1)/(Y+1).
    resultCombined$qValue<-qvalue::qvalue(resultCombined$pValue)$qvalues
    # eGenes<-(resultCombined%>%filter(qValue<=0.05))$gene_id

  }else if(approach=="Clipper"){
    # resultCombined<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=Clipper_B=20/_resultGenes.rds")
    # approach<-"Clipper"
    # B<-20
    # outputDir<-"~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.06.27_runPackage/_result/Lung/approach=Clipper_B=20/"
    # #To comment out.

    if(B==1){
      procedure<-"BC"
      h<-NULL
    }else{
      procedure<-"GZ"
      h<-1
    }

    temp<-Clipper::Clipper(score.exp=as.matrix(resultCombined[,5]),
                           score.back=as.matrix(resultCombined[,6:(4+1+B)]), #Key argument.
                           analysis="enrichment",
                           procedure=procedure,
                           contrast.score="max",
                           n.permutation=h)
    resultCombined$contrastScore<-temp$contrast.score.value
    resultCombined$qValue<-temp$q
    # eGenes<-(resultCombined%>%filter(qValue<=0.05))$gene_id
  }

  pathResultGenes<-paste0(outputDir,"_resultGenes.rds")
  saveRDS(resultCombined,pathResultGenes)
  unlink(paste0(outputDir,"_resultCombined.rds"))

  return(resultCombined) #26,095*27. The first four columns are chr, start, end, and gene_id. end is used as TSS. The last column is qValue.
}

















