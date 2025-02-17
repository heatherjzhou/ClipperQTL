ClipperQTL
================

We have shown that ClipperQTL performs as well as FastQTL and runs up to
500 times faster (2023). Without relying on GPUs, it is up to 30 times
more computationally efficient than tensorQTL, a GPU-based
implementation of FastQTL. Here we show how to install and use
ClipperQTL.

## 1. Installation

First, install qValue from Bioconductor and Clipper from GitHub.

``` r
#install.packages("BiocManager")
BiocManager::install("qvalue") #Needed for the standard variant of ClipperQTL.

#install.packages("devtools")
devtools::install_github("JSB-UCLA/Clipper") #Needed for the Clipper variant of ClipperQTL.
```

Next, install HTSlib following <http://www.htslib.org/download/>.

Lastly, install ClipperQTL. If prompted, make sure to update
RcppArmadillo to at least Version 0.12.0.1.0, otherwise the installation
of ClipperQTL will fail. Note that ClipperQTL is only maintained for use
in Linux.

``` r
devtools::install_github("heatherjzhou/ClipperQTL")
```

## 2. `ClipperQTL()` for identifying eGenes

The ClipperQTL package contains two main functions: `ClipperQTL()` for
identifying eGenes and `callSigGeneSNPPairs()` for identifying
significant gene-SNP pairs. Details of these functions can be found by
running `library(ClipperQTL)` followed by `?ClipperQTL` and
`?callSigGeneSNPPairs`.

`ClipperQTL()` requires three main pieces of input data: expression
data, covariate data, and genotype data (the file paths are specified by
`exprFile`, `covFile`, and `genotypeFile`, respectively; see example
code below). All three data sets must have the same format as the data
sets used in GTEx V8 analysis (2020)
(<https://gtexportal.org/home/datasets>).

Specifically, the expression file must be a .bed.gz file. The first four
columns must be chr, start, end, and gene_id (the exact column names do
not matter). Each remaining column must correspond to a sample. The
third column, end, will be used as the transcription start site
(following FastQTL; so make sure to put the desired transcription start
site in this column). The second column, start, will not be used. An
example expression file is provided in the folder named “example” in
this repository. This file is downloaded from
<https://gtexportal.org/home/datasets>.

``` r
#Take a look at the example expression file.

exprFile<-"~/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood.v8.normalized_expression.bed.gz"
temp1<-readr::read_delim(exprFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #20,315*674. Sample size is 670.
```

The covariate file must be a .txt file. The first column must be the
covariate ID (the exact column name does not matter). Each remaining
column must correspond to a sample. All constant covariates will be
filtered out before analysis. An example covariate file is provided in
the folder named “example” in this repository. The known covariates are
obtained from <https://gtexportal.org/home/datasets>. The inferred
covariates are the top expression PCs (2022)
(<https://github.com/heatherjzhou/PCAForQTL>).

``` r
#Take a look at the example covariate file.

covFile<-"~/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood.v8.covariates.txt"
temp2<-readr::read_delim(covFile,delim="\t",escape_double=FALSE,trim_ws=TRUE) #50*671.
```

The genotype file must be a .vcf.gz file. A .tbi index file must be
present in the same directory, which can be generated using HTSlib tabix
(<https://www.htslib.org/doc/tabix.html>). The non-empty genotype
entries must start with “0\|0”, “0\|1”, “1\|0”, or “1\|1” (trailing
characters are ok). The missing genotype entries will be imputed as
within-SNP averages. This file can contain more samples than the
expression file and the covariate file (which must have the same
samples). An example genotype file is not provided in this repository.
The GTEx V8 genotype file can be downloaded from the AnVIL repository
with an approved dbGaP application (see
<https://gtexportal.org/home/protectedDataAccess>).

The main method parameters of `ClipperQTL()` are `approach` and `B`. If
the sample size is less than or equal to 450, then we recommend setting
`approach="standard"` and `B=1000`. If the sample size is greater than
450, then you may set `approach="standard"` and `B=1000`, or, you may
set `approach="Clipper"` and `B=1` or `B` between 20 and 100 for faster
computational speed (2023). Regarding which variant of ClipperQTL should
be used when the sample size is large enough, we believe that if
computational efficiency is a priority, then the Clipper variant should
be used. However, if the majority of data sets in the study have smaller
sample sizes, then you may choose to use the standard variant on all
data sets for consistency.

`ClipperQTL()` outputs several files in the output directory. The most
important one is named “\_resultGenes.rds”, which can be read into R
with `readRDS()`. The first four columns are identical to the first four
columns in the expression file. The next `1+B` columns are the maximum
absolute correlations from the experimental round and the permutation
rounds. The second to last column is `pValue` (standard variant) or
`contrastScore` (Clipper variant). The last column is `qValue`. Each row
corresponds to a gene. The eGenes are those with `qValue` under the
target FDR threshold, e.g., 0.05.

``` r
#Example code.

exprFile<-"~/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood.v8.normalized_expression.bed.gz"
covFile<-"~/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood.v8.covariates.txt"
genotypeFile<-"~/_Data/2020.09.21_GTEx_V8/1_Raw/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
tabixProgram<-"~/_Applications/htslib-1.10.2/bin/tabix"

#Use the standard variant of ClipperQTL.
approach<-"standard"
B<-1000

# #Alternatively, use the Clipper variant of ClipperQTL for faster computational speed (only recommended for data sets with sample size greater than 450).
# approach<-"Clipper"
# B<-1 #B=1 or B between 20 and 100 is recommended.

outputDir<-paste0("~/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood_","approach=",approach,"_B=",B,"/")
if(!dir.exists(outputDir)) dir.create(outputDir)

library(ClipperQTL) #Loading ClipperQTL also loads dplyr. Loading dplyr is necessary for ClipperQTL() and callSigGeneSNPPairs() to run.
ClipperQTL(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
           approach,B,
           cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
           numOfChunksTarget=100,seed=1,numOfCores=5)
```

## 3. `callSigGeneSNPPairs()` for identifying significant gene-SNP pairs

Most arguments in `callSigGeneSNPPairs()` (from `exprFile` to
`MASamplesThreshold`) must be the same as those used in `ClipperQTL()`.
However, since `outputDir` will be used in C++ in this function, it must
be recognizable in C++ (see example code below).

The main method parameters of `callSigGeneSNPPairs()` are
`sigGeneSNPPairMethod` and `percent`. `sigGeneSNPPairMethod` must be
`"FastQTL"`, `"topPercent"`, or `NULL`. `percent` is only relevant if
`"topPercent"` is used.

If `sigGeneSNPPairMethod="FastQTL"`, then significant gene-SNP pairs
will be identified using the method in FastQTL (see Algorithm S2 of our
paper; inverse functions of cumulative distribution functions are
replaced by quantile functions).

If `sigGeneSNPPairMethod="topPercent"`, then the top 1% local common
SNPs (in terms of significance of association) will be identified as
significant for each eGene (if `percent=1`, that is).

If `sigGeneSNPPairMethod=NULL`, then if `approach="standard"`,
`sigGeneSNPPairMethod` will be set to `"FastQTL"`; if
`approach="Clipper"`, `sigGeneSNPPairMethod` will be set to
`"topPercent"`.

If `approach="standard"`, then `sigGeneSNPPairMethod` can be either
`"FastQTL"` or `"topPercent"`. If `approach="Clipper"`, then
`sigGeneSNPPairMethod` cannot be `"FastQTL"`.

`callSigGeneSNPPairs()` outputs several files in the output directory.
The most important one is named “\_resultSigGeneSNPPairs_FastQTL.rds” or
“\_resultSigGeneSNPPairs_topPercent.rds” (depending on
`sigGeneSNPPairMethod`), which can be read into R with `readRDS()`. Each
row corresponds to a significant gene-SNP pair.

``` r
#Example code.

#Follow up with the standard variant of ClipperQTL.
approach<-"standard"
B<-1000

# #Alternatively, follow up with the Clipper variant of ClipperQTL.
# approach<-"Clipper"
# B<-1

outputDir<-paste0("/home/heatherjzhou/2022.03.14_ClipperQTL/ClipperQTL/example/Whole_Blood_","approach=",approach,"_B=",B,"/") #"~" may not be recognized as the home directory in C++; you may need to spell out the directory name instead.

callSigGeneSNPPairs(exprFile,covFile,genotypeFile,tabixProgram,outputDir,
                    approach,
                    cisDistance=1e6,MAFThreshold=0.01,MASamplesThreshold=10,
                    numOfCores=1,
                    FDR_eGene=0.05,sigGeneSNPPairMethod=NULL,percent=1) #Let the function choose sigGeneSNPPairMethod automatically.
```

Please note that `callSigGeneSNPPairs()` only provides a rudimentary
analysis of significant gene-SNP pairs. For a more in-depth analysis,
you may try fine-mapping methods such as SuSiE (2020).

## 4. Memory usage

ClipperQTL should use no more than a few GB of memory per core on data
sets with under 1000 individuals and ~20,000 features (this is true for
both functions in ClipperQTL: `ClipperQTL()` and
`callSigGeneSNPPairs()`). The total memory usage is proportional to the
number of cores used. If memory is a concern, use fewer cores or the
default number of cores, which is 1 for both functions.

## Citation

To acknowledge this package or this tutorial, please cite our paper
(2023): <https://doi.org/10.1101/2023.08.28.555191>. For questions,
please email us at <lijy03@g.ucla.edu> or <heatherjzhou@ucla.edu>.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gtexconsortiumGTExConsortiumAtlas2020" class="csl-entry">

GTEx Consortium. 2020. “The GTEx Consortium Atlas of Genetic Regulatory
Effects Across Human Tissues.” *Science* 369 (6509): 1318–30.

</div>

<div id="ref-ongenFastEfficientQTL2016" class="csl-entry">

Ongen, Halit, Alfonso Buil, Andrew Anand Brown, Emmanouil T.
Dermitzakis, and Olivier Delaneau. 2016. “Fast and Efficient QTL Mapper
for Thousands of Molecular Phenotypes.” *Bioinformatics* 32 (10):
1479–85.

</div>

<div id="ref-taylor-weinerScalingComputationalGenomics2019"
class="csl-entry">

Taylor-Weiner, Amaro, François Aguet, Nicholas J. Haradhvala, Sager
Gosai, Shankara Anand, Jaegil Kim, Kristin Ardlie, Eliezer M. Van Allen,
and Gad Getz. 2019. “Scaling Computational Genomics to Millions of
Individuals with GPUs.” *Genome Biology* 20 (1): 228.

</div>

<div id="ref-wangSimpleNewApproach2020" class="csl-entry">

Wang, Gao, Abhishek Sarkar, Peter Carbonetto, and Matthew Stephens.
2020. “A Simple New Approach to Variable Selection in Regression, with
Application to Genetic Fine Mapping.” *Journal of the Royal Statistical
Society: Series B (Statistical Methodology)* 82 (5): 1273–1300.

</div>

<div id="ref-zhouClipperQTLUltrafastPowerful2023" class="csl-entry">

Zhou, Heather J., Xinzhou Ge, and Jingyi Jessica Li. 2023. “ClipperQTL:
Ultrafast and Powerful <span class="nocase">eGene</span> Identification
Method.” *bioRxiv*.

</div>

<div id="ref-zhouPCAOutperformsPopular2022b" class="csl-entry">

Zhou, Heather J., Lei Li, Yumei Li, Wei Li, and Jingyi Jessica Li. 2022.
“PCA Outperforms Popular Hidden Variable Inference Methods for Molecular
QTL Mapping.” *Genome Biology* 23 (1): 210.

</div>

</div>
