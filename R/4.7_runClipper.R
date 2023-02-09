






#There are different ways to run Clipper:
#Difference vs maximum constrast score.
#Number of nulls.
#Alpha.
#Using maxAbsCors, t-stats, or -log10(pValues).
runClipper<-function(n,k,maxNumOfNulls,resultCombined,outputDir){
  # n<-ncol(dataGeneExpressionFP)-4
  # k<-ncol(dataCovariates)
  # #To comment out.

  # #Convert maxAbsCors to t-stats. Optional.
  # tableMaxAbsCors<-as.matrix(resultCombined[,-(1:4)])
  # table_tStats<-tableMaxAbsCors*sqrt((n-2-k)/(1-tableMaxAbsCors^2)) #...*21. This is actually already the absolute value of the t-stats.
  # resultCombined[,-(1:4)]<-table_tStats

  #Convert maxAbsCors to -log10(pValues). Optional.
  tableMaxAbsCors<-as.matrix(resultCombined[,-(1:4)])
  table_tStats<-tableMaxAbsCors*sqrt((n-2-k)/(1-tableMaxAbsCors^2)) #...*21. This is actually already the absolute value of the t-stats.
  table_pValues<-2*pt(q=abs(table_tStats),df=n-2-k,lower.tail=FALSE,log.p=FALSE) #...*21.
  tableScores<--log10(table_pValues) #...*21. Each row corresponds to a gene. The columns are: exp, bg1, ..., bg20.
  resultCombined[,-(1:4)]<-tableScores

  for(numOfNulls in 1:maxNumOfNulls){
    # numOfNulls<-1

    temp<-Clipper::Clipper(score.exp=resultCombined$exp,
                           score.back=resultCombined[,(5+1):(5+numOfNulls)], #Key argument.
                           analysis="e",contrast.score="max")
    resultCombined$contrastScore<-temp$contrast.score.value
    resultCombined$qValue<-temp$q
    colnames(resultCombined)[(ncol(resultCombined)-1):ncol(resultCombined)]<-paste0(c("contrastScore","qValue"),"_",numOfNulls,"Nulls")
  }
  saveRDS(resultCombined,paste0(outputDir,"_resultCombined.rds"))

  return(resultCombined)



}













