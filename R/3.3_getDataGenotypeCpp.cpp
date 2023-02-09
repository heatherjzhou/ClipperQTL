#include <iostream>
#include <fstream>
#include <sstream>
#include <RcppArmadillo.h>
using namespace arma;

//The dimension of the temp file is numOfLines*numOfColumns.
//The row index is indexOfLine, and the column index is indexOfColumn.
//sampleIndices represents which columns of the temp file correspond to the samples in the expression and covariate data (need to subtract 1 before using in Cpp).
//[[Rcpp::export()]]
Rcpp::List getDataGenotypeCpp(const std::string chrCurr,
                              const vec geneTSSs,
                              const std::string genotypeFile,
                              const std::string tabixProgram,
                              const std::string tempFileName,

                              const vec sampleIndices,
                              const double MAFThreshold,
                              const int MASamplesThreshold){

  // auto start=std::chrono::high_resolution_clock::now();

  //Run tabix.
  int posBeginning=min(geneTSSs)-1e6; //74,786,699.
  if(posBeginning<0){ //posBeginning should not be smaller than 0.
    posBeginning=0;
  }
  int posEnd=max(geneTSSs)+1e6; //112,200,771.

  std::string region=chrCurr+":"+std::to_string(posBeginning)+"-"+std::to_string(posEnd); //chr1:74786699-112200771.
  std::string command=tabixProgram+" "+genotypeFile+" "+region+" > "+tempFileName;
  int statusCode=system(command.c_str()); //.c_str() returns a char*. system() requires char* input. Save return value in statusCode to avoid warning message.

  //Get numOfLines, the total number of lines in the temp file. 136,353.
  int numOfLines=0;

  std::ifstream myFile(tempFileName); //Create an ifstream object for reading from file.
  std::string myLine; //Create a string to store each line.
  while(getline(myFile,myLine)){ //Read the file line by line.
    numOfLines=numOfLines+1;
  }
  myFile.close(); //Close the file.

  //Create SNPInfo and X (filter samples, recode genotype data).
  Rcpp::CharacterMatrix SNPInfo(numOfLines,3); //136,353*3. The columns are: CHROM, POS, ID.
  mat X(numOfLines,sampleIndices.n_elem); //136,353*515.

  vec numsOfSamples0(numOfLines); //Vector of length 136,353. Calculate it here rather than later to save about half a second to a second.
  vec numsOfSamples1(numOfLines); //Vector of length 136,353.
  vec numsOfSamples2(numOfLines); //Vector of length 136,353.

  std::ifstream myFile2(tempFileName); //Create an ifstream object for reading from file. Cannot reuse myFile even if it is not closed.
  std::string myLine2;
  int indexOfLine=0;
  while(std::getline(myFile2,myLine2)){
    std::stringstream myStringStream(myLine2); //Create a stringstream object for reading from myLine2.
    std::string myWord;
    int indexOfColumn=0;
    int numOfSamplesProcessed=0;
    while(std::getline(myStringStream,myWord,'\t')){ //Double quotes do not work.
      if(indexOfColumn<=2){ //Fill SNPInfo.
        SNPInfo(indexOfLine,indexOfColumn)=myWord;
      }else if(indexOfColumn==(sampleIndices(numOfSamplesProcessed)-1)){ //Fill X.
        if(myWord=="0|0"){
          X(indexOfLine,numOfSamplesProcessed)=0;
          numsOfSamples0(indexOfLine)+=1;
        }else if(myWord=="0|1" || myWord=="1|0"){
          X(indexOfLine,numOfSamplesProcessed)=1;
          numsOfSamples1(indexOfLine)+=1;
        }else if(myWord=="1|1"){
          X(indexOfLine,numOfSamplesProcessed)=2;
          numsOfSamples2(indexOfLine)+=1;
        }else{
          X(indexOfLine,numOfSamplesProcessed)=datum::nan;
        }

        numOfSamplesProcessed+=1;
        if(numOfSamplesProcessed==sampleIndices.n_elem){ //If the sample size is 515, and numOfSamplesProcessed has just been incremented to 515, then break because sampleIndices(numOfSamplesProcessed) will not work.
          break;
        }
      }

      indexOfColumn+=1;
    } //End of while(std::getline(myStringStream,myWord,'\t')).

    indexOfLine+=1;
  } //End of while(std::getline(myFile2,myLine2)).
  myFile2.close(); //Close the file.
  remove(tempFileName.c_str()); //Delete the temp file.

  //Filter SNPs. If MAFThreshold and MASamples are 0, then the following code doesn't change X.
  // vec numsOfSamples0=conv_to<vec>::from(sum(X==0,1)); //1 means to return the statistic for each row. NaN==0 returns false. Use conv_to because otherwise I get a uvec rather than vec, which leads to MAFs being rounded to integers.
  // vec numsOfSamples1=conv_to<vec>::from(sum(X==1,1));
  // vec numsOfSamples2=conv_to<vec>::from(sum(X==2,1));
  vec numsOfSamplesWithGenotype=numsOfSamples0+numsOfSamples1+numsOfSamples2; //Vector of length 136,353.

  vec MAFs_pseudo=(numsOfSamples1+2*numsOfSamples2)/(2*numsOfSamplesWithGenotype); //Vector of length 136,353. Alternative allele frequencies.
  vec MAFs=min(MAFs_pseudo,1-MAFs_pseudo); //Vector of length 136,353.
  vec MASamples=min(numsOfSamples0,numsOfSamples2)+numsOfSamples1; //Vector of length 136,353.

  uvec indicesOfSNPsToKeep=find(MAFs>=MAFThreshold && MASamples>=MASamplesThreshold); //Vector of length 131,682.
  X=X.rows(indicesOfSNPsToKeep); //131,682*515.

  //Impute missing genotype data.
  vec averages=MAFs_pseudo(indicesOfSNPsToKeep)*2; //Vector of length 131,682.
  numsOfSamplesWithGenotype=numsOfSamplesWithGenotype(indicesOfSNPsToKeep); //Vector of length 131,682.
    for(int i=0;i<X.n_rows;i++){
      if(numsOfSamplesWithGenotype(i)<sampleIndices.n_elem){
        X.row(i).replace(datum::nan,averages(i));
      }
    }

  //Create SNPInfoFiltered.
  Rcpp::CharacterMatrix SNPInfoFiltered(X.n_rows,3); //131,682*3. The columns are: CHROM, POS, ID.
  for(int i=0;i<X.n_rows;i++){
    SNPInfoFiltered(i,Rcpp::_)=SNPInfo(indicesOfSNPsToKeep(i),Rcpp::_);
  }

  Rcpp::List toReturn=Rcpp::List::create(Rcpp::_("SNPInfo")=SNPInfoFiltered,
                                         Rcpp::_["X"]=X);

  // auto stop=std::chrono::high_resolution_clock::now();
  // auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // cout<<"The C++ code took "<<duration.count()<<" milliseconds.\n";

  return toReturn;
}








