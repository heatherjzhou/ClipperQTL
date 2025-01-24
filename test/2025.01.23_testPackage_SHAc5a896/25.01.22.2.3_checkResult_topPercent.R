library(dplyr,tidyr)
library(ggplot2)


#Lung.
#Method 1: ClipperQTL with approach="standard" and sigGeneSNPPairMethod="FastQTL".
#Method 2: ClipperQTL with approach="Clipper" and sigGeneSNPPairMethod="topPercent".


getOverlap<-function(vec1,vec2){
  overlap<-length(intersect(vec1,vec2))/min(length(vec1),length(vec2))
  return(overlap)
}


resultPairs1<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=standard_B=1000/_resultSigGeneSNPPairs_FastQTL.rds") #1,948,850*7.
resultPairs2<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=Clipper_B=1/_resultSigGeneSNPPairs_topPercent.rds") #1,075,778*7.

#Only keep pairs that correspond to shared eGenes.
resultGenes1<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=standard_B=1000/_resultGenes.rds") #26,095*1007.
resultGenes2<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=Clipper_B=1/_resultGenes.rds") #26,095*8.
eGenes1<-(resultGenes1%>%filter(qValue<=0.05))$gene_id #Vector of length 13,629.
eGenes2<-(resultGenes2%>%filter(qValue<=0.05))$gene_id #Vector of length 14,157.
getOverlap(eGenes1,eGenes2) #0.9824639.
eGenes_shared<-intersect(eGenes1,eGenes2) #13,390 shared eGenes.

resultPairs1<-resultPairs1%>%filter(gene_id%in%eGenes_shared) #1,947,626*7.
resultPairs2<-resultPairs2%>%filter(gene_id%in%eGenes_shared) #1,010,157*7.

#Create pair_id.
resultPairs1<-resultPairs1%>%mutate(pair_id=paste0(gene_id,"_",SNP_id)) #1,947,626*8.
resultPairs2<-resultPairs2%>%mutate(pair_id=paste0(gene_id,"_",SNP_id)) #1,010,157*8.



resultPairs1_summary<-resultPairs1%>%group_by(gene_id)%>%
  summarise(numOfESNPs=n()) #13,390*2.
hist(resultPairs1_summary$numOfESNPs,breaks=50) #Right-skewed.
min(resultPairs1_summary$numOfESNPs) #1.

resultPairs2_summary<-resultPairs2%>%group_by(gene_id)%>%
  summarise(numOfESNPs=n()) #13,390*2.
hist(resultPairs2_summary$numOfESNPs,breaks=50) #Normal with right tail.
min(resultPairs2_summary$numOfESNPs) #1.

summary<-resultPairs1_summary%>%left_join(resultPairs2_summary,by="gene_id")
summary%>%ggplot(aes(x=numOfESNPs.x,y=numOfESNPs.y))+
  geom_point()+
  geom_abline(slope=1,intercept=0) #topPercent usually identifies fewer eSNPs than FastQTL but sometimes identifies more.




#Check that for a given gene, the eSNPs identified by one method are a subset of the eSNPs identified by the other method. Success.
geneIDCurr<-"ENSG00000000003.14" #Shared eGene. Chosen from summary.
resultPairs1Sub<-resultPairs1%>%filter(gene_id==geneIDCurr) #162*8.
resultPairs2Sub<-resultPairs2%>%filter(gene_id==geneIDCurr) #35*8.

length(intersect(resultPairs1Sub$pair_id,resultPairs2Sub$pair_id)) #35 is good.




