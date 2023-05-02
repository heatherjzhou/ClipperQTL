library(dplyr)


source("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.03.20_combineChunks/23.03.25_func4.1_combineChunks_V1.R")
source("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.03.20_combineChunks/23.03.25_func4.2_combineChunks_V2.R")


tissueType<-"Lung" #Sample size is 515.
outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.03.20_combineChunks/_result/",tissueType,"/")
chunkInfo<-readRDS(paste0(outputDir,"_chunkInfo.rds")) #103*4.

#combineChunks_V1. write.table + read.delim.
if(FALSE){

  timeStart<-Sys.time()
  combineChunks_V1(outputDir,chunkInfo)

  timeEnd<-Sys.time()
  cat("Writing the result took: ")
  print(timeEnd-timeStart)
  #Time difference of 52.20569 secs on 2023/03/20 when FastQTL was using 20 cores.
  #Time difference of 32.60073 secs on 2023/03/21 when FastQTL was finished.


  timeStart<-Sys.time()

  pathToRead<-paste0(outputDir,"_resultCombined.txt")
  result<-read.delim(pathToRead)

  timeEnd<-Sys.time()
  cat("Reading in the result took: ")
  print(timeEnd-timeStart)
  #Time difference of 3.523671 mins on 2023/03/20 when FastQTL was using 20 cores.
  #Time difference of 1.836466 mins on 2023/03/21 when FastQTL was finished.
  #1.057613 mins on 2023/02/17.
}





#combineChunks_V2. saveRDS + readRDS.
if(TRUE){
  timeStart<-Sys.time()
  combineChunks_V2(outputDir,chunkInfo)

  timeEnd<-Sys.time()
  cat("Writing the result took: ")
  print(timeEnd-timeStart)
  #Time difference of 54.00183 secs on 2023/03/20 when FastQTL was using 20 cores.
  #Time difference of 28.69831 secs on 2023/03/20 when FastQTL was finished.


  timeStart<-Sys.time()

  pathToRead<-paste0(outputDir,"_resultCombined.rds")
  result<-readRDS(pathToRead)

  timeEnd<-Sys.time()
  cat("Reading in the result took: ")
  print(timeEnd-timeStart)
  #Time difference of 2.628547 secs on 2023/03/20 when FastQTL was using 20 cores.
  #Time difference of 1.46003 secs on 2023/03/20 when FastQTL was finished.


  # #Check result.
  # resultCombined_new<-readRDS("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.03.20_combineChunks/_result/Lung/_resultCombined.rds")
  # resultCombined_old<-readRDS("~/2022.03.14_ClipperQTL/_Work_2023.02.11/2023.02.16_ClipperQTLOnGTEx_archive/1.1_resultClipperQTL/Lung/BE_B=1000/_resultCombined.rds")
  # max(abs(resultCombined_new[,-(1:4)]-resultCombined_old[,-(1:4)])) #5.551115e-16.
}




