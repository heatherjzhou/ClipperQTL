library(dplyr,tidyr)
library(ggplot2)


#This file: check ClipperQTL() results. Success.


getOverlap<-function(vec1,vec2){
  overlap<-length(intersect(vec1,vec2))/min(length(vec1),length(vec2))
  return(overlap)
}


result_FastQTL<-read.delim("~/2022.03.14_ClipperQTL/_Work_2023.02.11/2023.03.03_realDataComparison/2.1_resultFastQTL/Lung/B=1000/Lung.genes.txt.gz") #26,095*19.
result_ClipperQTL_standard<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=standard_B=1000/_resultGenes.rds") #26,095*1007.
result_ClipperQTL_Clipper<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/Lung/approach=Clipper_B=1/_resultGenes.rds") #26,095*8.

eGenes_FastQTL<-(result_FastQTL%>%filter(qval<=0.05))$gene_id #Vector of length 13,578.
eGenes_ClipperQTL_standard<-(result_ClipperQTL_standard%>%filter(qValue<=0.05))$gene_id #Vector of length 13,629.
eGenes_ClipperQTL_Clipper<-(result_ClipperQTL_Clipper%>%filter(qValue<=0.05))$gene_id #Vector of length 14,157.

getOverlap(eGenes_FastQTL,eGenes_ClipperQTL_standard) #0.9923406.
getOverlap(eGenes_ClipperQTL_standard,eGenes_ClipperQTL_Clipper) #0.9824639.




