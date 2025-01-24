library(dplyr,tidyr)
library(ggplot2)


#Lung.
#Method 1: FastQTL with B=1000 and beta approximation.
#Method 2: ClipperQTL with approach="standard" and sigGeneSNPPairMethod="FastQTL".


getOverlap<-function(vec1,vec2){
  overlap<-length(intersect(vec1,vec2))/min(length(vec1),length(vec2))
  return(overlap)
}


resultPairs1<-read.delim("~/2022.03.14_ClipperQTL/_Work_2023.02.11/2023.04.10_significantESNPs/1.1_resultFastQTL/Lung/B=1000/Lung.signifpairs.txt.gz") #2,279,107*12.
resultPairs2<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=standard_B=1000/_resultSigGeneSNPPairs_FastQTL.rds") #1,948,850*7.

#Only keep pairs that correspond to shared eGenes.
resultGenes1<-read.delim("~/2022.03.14_ClipperQTL/_Work_2023.02.11/2023.04.10_significantESNPs/1.1_resultFastQTL/Lung/B=1000/Lung.genes.txt.gz") #26,095*19.
resultGenes2<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=standard_B=1000/_resultGenes.rds") #26,095*1007.
eGenes1<-(resultGenes1%>%filter(qval<=0.05))$gene_id #13,578 eGenes.
eGenes2<-(resultGenes2%>%filter(qValue<=0.05))$gene_id #13,629 eGenes.
getOverlap(eGenes1,eGenes2) #0.9923406.
eGenes_shared<-intersect(eGenes1,eGenes2) #13,474 shared eGenes.

resultPairs1<-resultPairs1%>%filter(gene_id%in%eGenes_shared) #2,278,369*12.
resultPairs2<-resultPairs2%>%filter(gene_id%in%eGenes_shared) #1,948,567*7.

#Create pair_id.
resultPairs1<-resultPairs1%>%mutate(pair_id=paste0(gene_id,"_",variant_id)) #2,278,369*13.
resultPairs2<-resultPairs2%>%mutate(pair_id=paste0(gene_id,"_",SNP_id)) #1,948,567*8.

#Overlap is high.
getOverlap(resultPairs1$pair_id,resultPairs2$pair_id) #0.9999954.



#For most shared eGenes, FastQTL identifies more eSNPs than ClipperQTL.
if(TRUE){
  resultPairs1_summary<-resultPairs1%>%group_by(gene_id)%>%
    summarise(numOfESNPs=n()) #13,474*2.
  min(resultPairs1_summary$numOfESNPs) #1.
  sum(resultPairs1_summary$numOfESNPs==1) #550 out of 13,474 eGenes have only 1 eSNP identified.
  resultPairs2_summary<-resultPairs2%>%group_by(gene_id)%>%
    summarise(numOfESNPs=n()) #13,474*2.
  min(resultPairs2_summary$numOfESNPs) #1.
  sum(resultPairs2_summary$numOfESNPs==1) #1282 out of 13,474 eGenes have only 1 eSNP identified.

  summary<-resultPairs1_summary%>%left_join(resultPairs2_summary,by="gene_id") #13,474*3.
  summary<-summary%>%left_join(resultGenes1%>%select(gene_id,beta_shape1,beta_shape2),by="gene_id") #13,474*5.

  sum(summary$numOfESNPs.x>summary$numOfESNPs.y) #12,077. Most shared eGenes have more eSNPs identified by FastQTL than by ClipperQTL.
  sum(summary$beta_shape1<1) #2152. Whether a gene has more eSNPs identified by FastQTL cannot be explained by beta_shape1.

  summary%>%ggplot(aes(x=numOfESNPs.x,y=numOfESNPs.y))+
    geom_point()+
    geom_abline(slope=1,intercept=0)

  summary<-summary%>%mutate(diffPercentage=(numOfESNPs.y-numOfESNPs.x)/numOfESNPs.x) #13,474*6.
  hist(summary$diffPercentage,breaks=50)
}

#Use a shared eGene as an example.
#The beta distribution from FastQTL is more right-skewed than the histogram from ClipperQTL.
#Therefore, the threshold for significant gene-SNP pairs from FastQTL is less significant,
#which explains why FastQTL identifies more eSNPs.
if(TRUE){
  #FastQTL beta_shape1 is slightly below or above 1.
  hist(resultGenes1$beta_shape1)
  range(resultGenes1$beta_shape1) #0.880226, 1.241690.

  #FastQTL beta_shape2 is always greater than 1.
  hist(resultGenes1$beta_shape2,breaks=50)
  range(resultGenes1$beta_shape2) #75.0351, 5465.1100.

  n<-515
  K<-52
  table_maxAbsCors<-as.matrix(resultGenes2[,-c(1:5,1006:1007)]) #26,095*1000.
  dim(table_maxAbsCors)
  table_tStat<-table_maxAbsCors*sqrt((n-2-K)/(1-table_maxAbsCors^2)) #tStat. This is actually already the absolute values of the t-stats.
  table_pVal<-2*pt(q=abs(table_tStat),df=n-2-K,lower.tail=FALSE,log.p=FALSE) #pVal.

  gene_id<-"ENSG00000001630.15" #Shared eGene. Chosen from summary.
  indexOfGene<-match(gene_id,resultGenes1$gene_id) #9871.
  indexOfGeneAlt<-match(gene_id,resultGenes2$gene_id) #9871.
  identical(indexOfGene,indexOfGeneAlt) #TRUE is good.

  beta_shape1<-resultGenes1$beta_shape1[indexOfGene] #1.02625. beta_shape1 from FastQTL.
  beta_shape2<-resultGenes1$beta_shape2[indexOfGene] #452.953. beta_shape2 from FastQTL.

  pVals<-table_pVal[indexOfGene,]
  hist(pVals,probability=TRUE,breaks=50)

  x<-seq(min(pVals),max(pVals),length=1000)
  f<-dbeta(x,beta_shape1,beta_shape2)
  lines(x,f)

  #The beta distribution from FastQTL is more right-skewed than the histogram from ClipperQTL.
  #Therefore, the threshold for significant gene-SNP pairs from FastQTL is less significant,
  #which explains why FastQTL identifies more eSNPs.
}






















