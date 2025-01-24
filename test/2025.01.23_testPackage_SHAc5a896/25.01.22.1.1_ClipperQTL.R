


tissueType<-"Lung"
numOfPCs<-44 #Chosen via BE.

exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
                 tissueType,".v8.normalized_expression.bed.gz")
covFile<-paste0("~/_Data/2020.09.21_GTEx_V8/3_De_Novo/2022.10.05_customCovariates_removeConstantCovariates/_Data_txt/",
                tissueType,"/",numOfPCs,"ExprPCs.txt")
genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"


#Run the standard variant of ClipperQTL().
if(TRUE){
  approach<-"standard" #Key.
  B<-1000 #Key.

  outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/",tissueType,"/")
  if(!dir.exists(outputDir)) dir.create(outputDir)
  outputDir<-paste0(outputDir,"approach=",approach,"_B=",B,"/")
  if(!dir.exists(outputDir)) dir.create(outputDir)

  library(ClipperQTL) #Loading ClipperQTL also loads dplyr. Loading dplyr is necessary for ClipperQTL() and callSigGeneSNPPairs() to run.
  ClipperQTL(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
             approach,B,
             cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
             numOfChunksTarget=100,seed=1,numOfCores=5)
}


#Run the Clipper variant of ClipperQTL().
if(TRUE){
  approach<-"Clipper" #Key.
  B<-1 #Key.

  outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/test/2025.01.23_testPackage_SHAc5a896/1.1_resultClipperQTL/",tissueType,"/")
  if(!dir.exists(outputDir)) dir.create(outputDir)
  outputDir<-paste0(outputDir,"approach=",approach,"_B=",B,"/")
  if(!dir.exists(outputDir)) dir.create(outputDir)

  library(ClipperQTL) #Loading ClipperQTL also loads dplyr. Loading dplyr is necessary for ClipperQTL() and callSigGeneSNPPairs() to run.
  ClipperQTL(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
             approach,B,
             cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
             numOfChunksTarget=100,seed=1,numOfCores=5)
}



# Rscript 25.01.22.1.1_ClipperQTL.R > 1.1_resultClipperQTL/_log.txt 2> 1.1_resultClipperQTL/_log2.txt




