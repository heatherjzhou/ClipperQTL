


prepareSampleIndices<-function(genotypeFile,tabixProgram,outputDir,
                               sampleNames){
  # tissueType<-"Lung" #Sample size is 515.
  #
  # genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
  # tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"
  # outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/R/_temp/",tissueType,"/")
  #
  # #Load sampleNames.
  # exprFile<-paste0("~/_Data/2020.09.21_GTEx_V8/1_Raw/Single_Tissue_cis_QTL_Data/GTEx_Analysis_v8_eQTL_expression_matrices/",
  #                  tissueType,".v8.normalized_expression.bed.gz")
  # dataGeneExpressionFP<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  # colnames(dataGeneExpressionFP)[1:4]<-c("chr","start","end","gene_id") #26,095*519. The first four columns are chr, start, end, and gene_id. end is used as TSS.
  # sampleNames<-colnames(dataGeneExpressionFP)[-(1:4)] #Vector of length 515.
  # #To comment out.

  filename<-paste0(outputDir,"_genotypeColumnNames.txt")
  command<-paste(tabixProgram,genotypeFile,"-H | tail -1 >",filename) #-H prints the header lines, i.e., the meta lines.
  system(command) #Invoke system command.

  temp<-readr::read_delim(filename,delim="\t",escape_double=FALSE,trim_ws=TRUE) #Import from text (base) does not work. Import from text (readr) works.
  genotypeColumnNames<-colnames(temp) #Vector of length 847.
  # genotypeColumnNames[1]<-"CHROM" #Remove the hashtag. This is unnecessary for the purpose of getting sampleIndices.
  # saveRDS(genotypeColumnNames,paste0(outputDir,"_genotypeColumnNames.rds"))

  unlink(filename)

  sampleIndices<-match(sampleNames,genotypeColumnNames) #Vector of length 515.
  if(anyNA(sampleIndices)) stop("At least one sample name present in the expression data is not present in the genotype data.\n")
  if(is.unsorted(sampleIndices)) stop("Sample names in the expression data must be sorted in the same way as sample names in the genotype data.\n")

  return(sampleIndices)
}




